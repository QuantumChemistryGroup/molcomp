#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <omp.h>

#include "mc_types.h"
#include "mc_atominfo.h"
#include "mc_getconnect.h"
#include "mc_io.h"
#include "mc_core.h"

enum {
    MC_NAME_SZ = 256,
    MC_MODE_SZ = 32,
    MC_LINE_SZ = 512,
    MC_PATH_SZ = 512,
    MC_COPY_SZ = 8192,
    MC_DESC_SZ = 4,
    MC_OUT_SZ  = 3
};

typedef struct {
    char name[MC_NAME_SZ];
    long int orig_idx;
    double *xyz;              /* flat coordinates [3 * natoms] */
    double desc[MC_DESC_SZ];  /* E, NucRep, RMSD_ref1, RMSD_ref2 */
    int ndesc;
    bool removed;
} ConfItem;

typedef struct {
    int idx1;
    int idx2;
    double out[MC_OUT_SZ];
    double dE_kcal;
} DupRecord;

typedef struct {
    char mode[MC_MODE_SZ];
    McOptions mc_opts;
    double value;
    double rmsd_thld;
    double maxdist_thld;
    bool readtopo;
} DedupConfig;

/* =========================================================================
 * Compact open-addressing hash set for tested pairs (avoids O(N^2) bitset)
 * ========================================================================= */
typedef struct {
    uint64_t *keys;
    size_t cap;   /* power of 2 */
    size_t count;
} PairSet;

static uint64_t hash64(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

static void pairset_init(PairSet *ps, size_t init_cap) {
    size_t cap = 1024;
    while (cap < init_cap) cap <<= 1;
    ps->keys = calloc(cap, sizeof(uint64_t));
    ps->cap = cap;
    ps->count = 0;
}

static void pairset_grow(PairSet *ps) {
    size_t old_cap = ps->cap;
    uint64_t *old_keys = ps->keys;
    size_t new_cap = old_cap << 1;
    uint64_t *new_keys = calloc(new_cap, sizeof(uint64_t));
    size_t mask = new_cap - 1;

    for (size_t i = 0; i < old_cap; i++) {
        uint64_t k = old_keys[i];
        if (k != 0) {
            size_t slot = (size_t)hash64(k) & mask;
            while (new_keys[slot] != 0) {
                slot = (slot + 1) & mask;
            }
            new_keys[slot] = k;
        }
    }
    free(old_keys);
    ps->keys = new_keys;
    ps->cap = new_cap;
}

/* Returns true if pair (a, b) was already present; otherwise inserts it */
static bool pairset_test_and_insert(PairSet *ps, long int a, long int b) {
    if (ps->count * 10 >= ps->cap * 7) {
        pairset_grow(ps);
    }
    /* +1 ensures key is never 0 even for (0, 0) */
    uint64_t key = (((uint64_t)(uint32_t)a) << 32) | ((uint64_t)(uint32_t)b + 1ULL);
    size_t mask = ps->cap - 1;
    size_t slot = (size_t)hash64(key) & mask;

    while (ps->keys[slot] != 0) {
        if (ps->keys[slot] == key) return true;
        slot = (slot + 1) & mask;
    }
    ps->keys[slot] = key;
    ps->count++;
    return false;
}

static void pairset_free(PairSet *ps) {
    free(ps->keys);
    ps->keys = NULL;
    ps->cap = 0;
    ps->count = 0;
}

/* =========================================================================
 * Helpers: natural sort, energy parser, options, NucRep
 * ========================================================================= */
static int natural_cmp(const char *a, const char *b) {
    while (*a && *b) {
        if (isdigit((unsigned char)*a) && isdigit((unsigned char)*b)) {
            char *ea, *eb;
            unsigned long va = strtoul(a, &ea, 10);
            unsigned long vb = strtoul(b, &eb, 10);
            if (va != vb) return (va < vb) ? -1 : 1;
            a = ea;
            b = eb;
        } else {
            if (*a != *b) return ((unsigned char)*a < (unsigned char)*b) ? -1 : 1;
            a++;
            b++;
        }
    }
    return (*a) ? 1 : ((*b) ? -1 : 0);
}

static int cmp_str_natural(const void *a, const void *b) {
    const char *sa = *(const char * const *)a;
    const char *sb = *(const char * const *)b;
    return natural_cmp(sa, sb);
}

static bool parse_energy_line(const char *line, double *out_e, bool *is_kcal) {
    *is_kcal = false;
    if (!line) return false;

    char buf[MC_LINE_SZ];
    strncpy(buf, line, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
char *tok = strtok(buf, " \t\r\n");
if (!tok) return false;
{
char *endp = NULL;
double val = strtod(tok, &endp);
if (endp != tok && *endp == '\0' && isfinite(val)) {
*out_e = val;
if (strstr(line, "kcal") != NULL) *is_kcal = true;
return true;
}
}

    const char *p = strstr(line, "E =");
    if (p) {
        p += 3;
    } else {
        p = strstr(line, "E=");
        if (p) p += 2;
    }
    if (p) {
        while (*p == ' ' || *p == '\t') p++;
        char *endp = NULL;
        double val = strtod(p, &endp);
        if (endp != p && isfinite(val)) {
            *out_e = val;
            if (strstr(line, "kcal") != NULL) *is_kcal = true;
            return true;
        }
    }
    return false;
}

static void init_dedup_config(DedupConfig *cfg, const char *mode, bool readtopo) {
    memset(cfg, 0, sizeof(*cfg));
    snprintf(cfg->mode, sizeof(cfg->mode), "%s", mode);
    cfg->readtopo = readtopo;
    cfg->mc_opts = MC_DEFAULT_OPTIONS();
    cfg->mc_opts.STEREO = false;

    if (strcmp(mode, "rmsd-kabsch") == 0) {
        cfg->mc_opts.CONNECT = false;
        strcpy(cfg->mc_opts.ALG, "kabsch");
        cfg->value = 0.125;
    } else if (strcmp(mode, "rmsd-hungarian") == 0) {
        cfg->mc_opts.CONNECT = false;
        strcpy(cfg->mc_opts.ALG, "std");
        cfg->value = 0.125;
    } else if (strcmp(mode, "rmsd-graph") == 0) {
    cfg->mc_opts.CONNECT = true;
    strcpy(cfg->mc_opts.ALG, "std");
    cfg->value = 0.125;
    } else if (strcmp(mode, "angle") == 0) {
        cfg->mc_opts.CONNECT = true;
        strcpy(cfg->mc_opts.ALG, "std");
        cfg->value = 30.0;
        cfg->rmsd_thld = 0.5;
        cfg->maxdist_thld = 1.5;
    } else if (strcmp(mode, "dist") == 0) {
        printf("Dist is not implemented yet. Exit.\n");
        exit(0);
    } else {
        printf("Unknown mode '%s'. Use: rmsd-kabsch, rmsd-hungarian, rmsd-graph, angle\n", mode);
        exit(0);
    }
}

static bool is_duplicate(const DedupConfig *cfg, const McResult *res) {
    if (strcmp(cfg->mode, "rmsd-graph") == 0 ||
        strcmp(cfg->mode, "rmsd-hungarian") == 0 ||
        strcmp(cfg->mode, "rmsd-kabsch") == 0) {
        return (res->rmsd < cfg->value);
    } else if (strcmp(cfg->mode, "angle") == 0) {
        return (res->rmsd < cfg->rmsd_thld &&
                res->dmax < cfg->maxdist_thld &&
                res->amax < cfg->value);
    }
    return false;
}

static double calc_nucrep_flat(const double *xyz, const double *Z_arr, long int natoms) {
    double E = 0.0;
    const double ang_to_bohr = 1.8897261259077824;
    for (long int i = 0; i < natoms - 1; i++) {
        double q1 = Z_arr[i];
        double xi = xyz[3 * i + 0];
        double yi = xyz[3 * i + 1];
        double zi = xyz[3 * i + 2];
        for (long int j = i + 1; j < natoms; j++) {
            double dx = xi - xyz[3 * j + 0];
            double dy = yi - xyz[3 * j + 1];
            double dz = zi - xyz[3 * j + 2];
            double d = sqrt(dx * dx + dy * dy + dz * dz) * ang_to_bohr;
            if (d > 1e-12) {
                E += (q1 * Z_arr[j]) / d;
            }
        }
    }
    return E;
}

static ConfItem *read_multixyz(
    const char *filename,
    long int *out_nconfs,
    long int *out_natoms,
    char (**out_elements)[4],
    bool *out_mix,
    bool *out_kcal
) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Cannot open file %s\n", filename);
        exit(0);
    }

    long int cap = 256;
    long int count = 0;
    ConfItem *list = calloc((size_t)cap, sizeof(ConfItem));
    long int expected_natoms = -1;
    char (*elements)[ 4 ] = NULL;
    bool any_missing_e = false;
    bool detected_kcal = false;

    char line[MC_LINE_SZ];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0') continue;

        long int nat = strtol(p, NULL, 10);
        if (nat <= 0) break;

        if (expected_natoms < 0) {
            expected_natoms = nat;
            elements = calloc((size_t)nat, sizeof(*elements));
        } else if (nat != expected_natoms) {
            printf("Error: inconsistent NATOMS in %s (%ld vs %ld)\n", filename, nat, expected_natoms);
            exit(0);
        }

        if (!fgets(line, sizeof(line), fp)) break;
        double e_val = 0.0;
        bool is_kcal = false;
        bool has_e = parse_energy_line(line, &e_val, &is_kcal);
        if (!has_e) {
            any_missing_e = true;
        } else if (is_kcal) {
            detected_kcal = true;
        }

        double *xyz = malloc((size_t)(3 * nat) * sizeof(double));
        for (long int a = 0; a < nat; a++) {
            if (!fgets(line, sizeof(line), fp)) {
                printf("Unexpected EOF in %s at frame %ld\n", filename, count + 1);
                exit(0);
            }
            char el_buf[8] = "";
            if (sscanf(line, "%7s %lf %lf %lf",
                       el_buf,
                       &xyz[3 * a + 0],
                       &xyz[3 * a + 1],
                       &xyz[3 * a + 2]) != 4) {
                printf("Bad coordinate line in %s at frame %ld: %s\n", filename, count + 1, line);
                exit(0);
            }
            if (count == 0) {
                snprintf(elements[a], 4, "%s", el_buf);
            } else if (strncmp(elements[a], el_buf, 3) != 0) {
                printf("Error: atom order mismatch in %s at frame %ld atom %ld (%s vs %s)\n",
                       filename, count + 1, a + 1, el_buf, elements[a]);
                exit(1);
            }
        }

        if (count >= cap) {
            cap *= 2;
            list = realloc(list, (size_t)cap * sizeof(ConfItem));
        }
        memset(&list[count], 0, sizeof(ConfItem));
        snprintf(list[count].name, sizeof(list[count].name), "conf_%ld.xyz", count + 1);
        list[count].orig_idx = count;
        list[count].xyz = xyz;
        *(list[count].desc + 0) = has_e ? e_val : 0.0;
        list[count].ndesc = has_e ? 1 : 0;
        count++;
    }

    fclose(fp);
    *out_nconfs = count;
    *out_natoms = (expected_natoms > 0) ? expected_natoms : 0;
    *out_elements = elements;
    *out_mix = any_missing_e;
    *out_kcal = detected_kcal;
    return list;
}

static ConfItem *read_dir_xyz(
    long int *out_nconfs,
    long int *out_natoms,
    char (**out_elements)[4],
    bool *out_mix,
    bool *out_kcal
) {
    DIR *d = opendir(".");
    if (!d) {
        printf("Cannot open current directory.\n");
        exit(0);
    }

    char **fnames = NULL;
    long int nf = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        size_t len = strlen(ent->d_name);
        if (len > 4 && strcasecmp(ent->d_name + len - 4, ".xyz") == 0) {
            fnames = realloc(fnames, (size_t)(nf + 1) * sizeof(char *));
            fnames[nf] = strdup(ent->d_name);
            nf++;
        }
    }
    closedir(d);

    if (nf == 0) {
        printf("No .xyz files found in current directory.\n");
        exit(0);
    }

    qsort(fnames, (size_t)nf, sizeof(char *), cmp_str_natural);

    ConfItem *list = calloc((size_t)nf, sizeof(ConfItem));
    long int expected_natoms = -1;
    char (*elements)[4] = NULL;
    bool any_missing_e = false;
    bool detected_kcal = false;

    for (long int i = 0; i < nf; i++) {
        FILE *fp = fopen(fnames[i], "r");
        if (!fp) {
         printf("Cannot open file %s\n", fnames[i]);
         exit(1);
     }

        char line[MC_LINE_SZ];
        double e_val = 0.0;
        bool is_kcal = false;
        bool has_e = false;
        if (fgets(line, sizeof(line), fp) && fgets(line, sizeof(line), fp)) {
            has_e = parse_energy_line(line, &e_val, &is_kcal);
        }
        if (!has_e) {
            any_missing_e = true;
            printf("Warning: E= NaN in %s\n", fnames[i]);
        } else if (is_kcal) {
            detected_kcal = true;
        }

        long int nat = 0;
        McAtom *atoms = malloc(1 * sizeof(McAtom));
        MC_FREAD(fp, &atoms, &nat, false);
        fclose(fp);

        if (expected_natoms < 0) {
            expected_natoms = nat;
            elements = calloc((size_t)nat, sizeof(*elements));
            for (long int a = 0; a < nat; a++) {
                strncpy(elements[a], atoms[a].El, 3);
                elements[a][3] = '\0';
            }
        } else if (nat != expected_natoms) {
            printf("Error: NAT1 != NAT2 in %s (%ld vs %ld)\n", fnames[i], nat, expected_natoms);
            exit(0);
        } else {
            for (long int a = 0; a < nat; a++) {
                if (strncmp(elements[a], atoms[a].El, 3) != 0) {
                    printf("Error: atom order mismatch in %s at atom %ld (%s vs %s)\n",
                           fnames[i], a + 1, atoms[a].El, elements[a]);
                    exit(1);
                }
            }
        }

        double *xyz = malloc((size_t)(3 * nat) * sizeof(double));
        for (long int a = 0; a < nat; a++) {
            xyz[3 * a + 0] = *(atoms[a].xyz + 0);
            xyz[3 * a + 1] = *(atoms[a].xyz + 1);
            xyz[3 * a + 2] = *(atoms[a].xyz + 2);
        }
        free(atoms);

        strncpy(list[i].name, fnames[i], sizeof(list[i].name) - 1);
        list[i].orig_idx = i;
        list[i].xyz = xyz;
        *(list[i].desc + 0) = has_e ? e_val : 0.0;
        list[i].ndesc = has_e ? 1 : 0;

        free(fnames[i]);
    }
    free(fnames);

    *out_nconfs = nf;
    *out_natoms = expected_natoms;
    *out_elements = elements;
    *out_mix = any_missing_e;
    *out_kcal = detected_kcal;
    return list;
}

static int G_SORT_KEY = 0;
static bool G_SORT_REV = false;

static int cmp_conf_ptr(const void *a, const void *b) {
    const ConfItem *ca = *(const ConfItem * const *)a;
    const ConfItem *cb = *(const ConfItem * const *)b;
    double va = ca->desc[G_SORT_KEY];
    double vb = cb->desc[G_SORT_KEY];
    if (va < vb) return G_SORT_REV ? 1 : -1;
    if (va > vb) return G_SORT_REV ? -1 : 1;
    return natural_cmp(ca->name, cb->name);
}

static int cmp_duprec_canon(const void *a, const void *b) {
    const DupRecord *da = (const DupRecord *)a;
    const DupRecord *db = (const DupRecord *)b;
    if (da->idx1 != db->idx1) return (da->idx1 < db->idx1) ? -1 : 1;
    if (da->idx2 != db->idx2) return (da->idx2 < db->idx2) ? -1 : 1;
    return 0;
}
static int cmp_duprec(const void *a, const void *b) {
const DupRecord *da = (const DupRecord *)a;
const DupRecord *db = (const DupRecord *)b;
if (da->dE_kcal < db->dE_kcal) return -1;
if (da->dE_kcal > db->dE_kcal) return 1;
return cmp_duprec_canon(a, b);
}

static bool path_exists(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0);
}

static void copy_file_raw(const char *src, const char *dst) {
    FILE *fi = fopen(src, "rb");
    if (!fi) return;
    FILE *fo = fopen(dst, "wb");
    if (!fo) { fclose(fi); return; }
    char buf[MC_COPY_SZ];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), fi)) > 0) {
        fwrite(buf, 1, n, fo);
    }
    fclose(fi);
    fclose(fo);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s <conformers.xyz> <mode> [readtopo]\n", *(argv + 0));
        printf("  %s <mode> [readtopo]\n", *(argv + 0));
        printf("Modes: rmsd-kabsch, rmsd-hungarian, rmsd-graph, angle\n");
        return 0;
    }

    bool single_file = false;
    char input_file[MC_PATH_SZ] = "";
    char mode[MC_MODE_SZ] = "";
    bool readtopo = false;

    if (strstr(*(argv + 1), ".xyz") != NULL || strstr(*(argv + 1), ".XYZ") != NULL) {
        if (argc < 3) {
            printf("Please specify mode after xyz file.\n");
            return 0;
        }
        single_file = true;
        strncpy(input_file, *(argv + 1), sizeof(input_file) - 1);
        strncpy(mode, *(argv + 2), sizeof(mode) - 1);
    } else {
        single_file = false;
        strncpy(mode, *(argv + 1), sizeof(mode) - 1);
    }

    for (int i = 0; mode[i]; i++) mode[i] = (char)tolower((unsigned char)mode[i]);

    if (strcasecmp(*(argv + argc - 1), "readtopo") == 0) {
        readtopo = true;
    }

    if (readtopo && !path_exists("topo.txt")) {
        printf("Your topo.txt file is missing. Stop.\n");
        return 0;
    }

    char base_name[MC_PATH_SZ] = "";
    char out_unique_dir[MC_PATH_SZ] = "";
    char out_unique_xyz[MC_PATH_SZ] = "";
    char out_dup_txt[MC_PATH_SZ] = "";
    char out_dupxyz_txt[MC_PATH_SZ] = "";

    if (!single_file) {
        snprintf(out_unique_dir, sizeof(out_unique_dir), "unique-%s", mode);
        snprintf(out_dup_txt, sizeof(out_dup_txt), "dup-%s.txt", mode);
        snprintf(out_dupxyz_txt, sizeof(out_dupxyz_txt), "DUPXYZ-%s.txt", mode);
        if (path_exists(out_unique_dir)) {
            printf("Folder %s exists. Remove it first.\n", out_unique_dir);
            return 0;
        }
        if (path_exists(out_dup_txt)) {
            printf("File %s exists. Remove it first.\n", out_dup_txt);
            return 0;
        }
    } else {
        size_t len = strlen(input_file);
        size_t blen = (len > 4) ? (len - 4) : len;
        if (blen >= sizeof(base_name)) blen = sizeof(base_name) - 1;
        memcpy(base_name, input_file, blen);
        base_name[blen] = '\0';

        snprintf(out_unique_xyz, sizeof(out_unique_xyz), "%s-unique-%s.xyz", base_name, mode);
        snprintf(out_dup_txt, sizeof(out_dup_txt), "%s-dup-%s.txt", base_name, mode);
        snprintf(out_dupxyz_txt, sizeof(out_dupxyz_txt), "%s-DUPXYZ-%s.txt", base_name, mode);
        if (path_exists(out_unique_xyz)) {
            printf("File %s exists. Remove it first.\n", out_unique_xyz);
            return 0;
        }
        if (path_exists(out_dup_txt)) {
            printf("File %s exists. Remove it first.\n", out_dup_txt);
            return 0;
        }
    }

    DedupConfig cfg;
    init_dedup_config(&cfg, mode, readtopo);

    long int nconfs = 0;
    long int natoms = 0;
    char (*elements)[ 4 ] = NULL;
    bool mix_mode = false;
    bool is_kcal_units = false;

    ConfItem *confs = single_file
        ? read_multixyz(input_file, &nconfs, &natoms, &elements, &mix_mode, &is_kcal_units)
        : read_dir_xyz(&nconfs, &natoms, &elements, &mix_mode, &is_kcal_units);

    if (nconfs <= 0) {
        printf("No conformers loaded.\n");
        return 0;
    }

    double *Z_arr = malloc((size_t)natoms * sizeof(double));
    for (long int a = 0; a < natoms; a++) {
        Z_arr[a] = (double)MC_ATOMINFO(elements[a]).Z;
    }

    long int *nn_ref = NULL;
    long int (*nbrs_ref)[ 23 ] = NULL;
    if (cfg.mc_opts.CONNECT) {
        McAtom *ref_atoms = calloc((size_t)natoms, sizeof(McAtom));
        for (long int a = 0; a < natoms; a++) {
            strcpy(ref_atoms[a].El, elements[a]);
            *(ref_atoms[a].xyz + 0) = *(confs->xyz + 3 * a + 0);
            *(ref_atoms[a].xyz + 1) = *(confs->xyz + 3 * a + 1);
            *(ref_atoms[a].xyz + 2) = *(confs->xyz + 3 * a + 2);
        }
        if (readtopo) {
            FILE *ft = fopen("topo.txt", "r");
            MC_TXYZREAD(ft, &ref_atoms, natoms);
            fclose(ft);
        } else {
            MC_GETCONNECT(&ref_atoms, natoms);
        }
        MC_SYMMETRYCHECK(&ref_atoms, natoms);

        nn_ref = calloc((size_t)natoms, sizeof(long int));
        nbrs_ref = calloc((size_t)natoms, sizeof(*nbrs_ref));
        for (long int a = 0; a < natoms; a++) {
            nn_ref[a] = ref_atoms[a].nn;
            for (long int k = 0; k < ref_atoms[a].nn; k++) {
                nbrs_ref[a][k] = ref_atoms[a].n[k];
            }
        }
        free(ref_atoms);
    }
    McPrecomp *pc = MC_PRECOMP_CREATE(natoms, elements, (confs + 0)->xyz, nn_ref, nbrs_ref, &cfg.mc_opts);
    if (mix_mode) {
        for (long int c = 0; c < nconfs; c++) {
            confs[c].ndesc = 0;
        }
    }

    #pragma omp parallel for schedule(static)
    for (long int c = 0; c < nconfs; c++) {
        double nr = calc_nucrep_flat(confs[c].xyz, Z_arr, natoms);
        confs[c].desc[confs[c].ndesc] = nr;
    }
    for (long int c = 0; c < nconfs; c++) confs[c].ndesc++;

    McOptions rmsd_opts = MC_DEFAULT_OPTIONS();
    rmsd_opts.CONNECT = false;
    rmsd_opts.STEREO = false;
    strcpy(rmsd_opts.ALG, "std");

    #pragma omp parallel for schedule(dynamic, 4)
    for (long int c = 0; c < nconfs; c++) {
        McResult r = MC_COMPARE_FLAT(natoms, elements, NULL, NULL, confs[c].xyz, confs->xyz, &rmsd_opts);
        confs[c].desc[confs[c].ndesc] = (r.status == 0) ? r.rmsd : 0.0;
    }
    for (long int c = 0; c < nconfs; c++) confs[c].ndesc++;

    if (mix_mode) {
        #pragma omp parallel for schedule(dynamic, 4)
        for (long int c = 0; c < nconfs; c++) {
            McResult r = MC_COMPARE_FLAT(natoms, elements, NULL, NULL, confs[c].xyz, confs[nconfs - 1].xyz, &rmsd_opts);
            confs[c].desc[confs[c].ndesc] = (r.status == 0) ? r.rmsd : 0.0;
        }
        for (long int c = 0; c < nconfs; c++) confs[c].ndesc++;
    }
    int t_desc = confs->ndesc;

    ConfItem **active = malloc((size_t)nconfs * sizeof(ConfItem *));
    long int n_active = nconfs;
    for (long int c = 0; c < nconfs; c++) active[c] = &confs[c];

    G_SORT_KEY = 0;
    G_SORT_REV = false;
    qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);

    PairSet tested_pairs;
    pairset_init(&tested_pairs, (size_t)(nconfs > 1024 ? nconfs * 16 : 16384));

    DupRecord *dups = NULL;
    long int ndups = 0;
    long int dups_cap = 0;
    long int total_tested = 0;

    int chunk_base = 24;
    int chunk_sz = chunk_base;
    int iter_max = 3000;
    int perm_max = 8;
    int pm = 0;
    unsigned int rng_seed = 42;

    printf("Parallel starts\n\n");

    for (int it = 0; it < iter_max; it++) {
        long int prev_dups = ndups;
        printf("Iteration: %d  permutation: %d  (%d)\n", it, pm, chunk_sz);

        long int n_chunks = (n_active + chunk_sz - 1) / chunk_sz;
        DupRecord *iter_dups = NULL;
        long int n_iter_dups = 0;
        long int iter_cap = 0;

        #pragma omp parallel
        {
            DupRecord *loc_dups = NULL;
            long int loc_ndups = 0;
            long int loc_cap = 0;
            long int loc_tested = 0;

            #pragma omp for schedule(dynamic, 1)
            for (long int ch = 0; ch < n_chunks; ch++) {
                long int start = ch * chunk_sz;
                long int end = start + chunk_sz;
                if (end > n_active) end = n_active;
                long int len = end - start;
                if (len <= 1) continue;

                bool *loc_removed = calloc((size_t)len, sizeof(bool));

                for (long int i = 0; i < len - 1; i++) {
                    if (loc_removed[i]) continue;
                    ConfItem *ci = active[start + i];

                    for (long int m = i + 1; m < len; m++) {
                        if (loc_removed[m]) continue;
                        ConfItem *cm = active[start + m];

                        ConfItem *keep = ci;
                        ConfItem *cand = cm;
                        if (*(cm->desc + 0) < *(ci->desc + 0)) {
                            keep = cm;
                            cand = ci;
                        } else if (*(cm->desc + 0) == *(ci->desc + 0)) {
                            if (natural_cmp(cm->name, ci->name) < 0) {
                                keep = cm;
                                cand = ci;
                            }
                        }

                        long int id_a = (keep->orig_idx < cand->orig_idx) ? keep->orig_idx : cand->orig_idx;
                        long int id_b = (keep->orig_idx < cand->orig_idx) ? cand->orig_idx : keep->orig_idx;

                        bool already = false;
                        #pragma omp critical(pair_cache)
                        {
                            already = pairset_test_and_insert(&tested_pairs, id_a, id_b);
                        }
                        if (already) continue;

                        loc_tested++;
                        McResult res = MC_COMPARE_PRECOMP(pc, keep->xyz, cand->xyz, &cfg.mc_opts);
                        if (res.status == 0 && is_duplicate(&cfg, &res)) {
                            double dE = fabs(*(cand->desc + 0) - *(keep->desc + 0));
                            double dE_kcal = is_kcal_units ? dE : (dE * 627.51);

                            if (loc_ndups >= loc_cap) {
                                loc_cap = (loc_cap == 0) ? 16 : loc_cap * 2;
                                loc_dups = realloc(loc_dups, (size_t)loc_cap * sizeof(DupRecord));
                            }
                            loc_dups[loc_ndups].idx1 = (int)keep->orig_idx;
                            loc_dups[loc_ndups].idx2 = (int)cand->orig_idx;
                            *(loc_dups[loc_ndups].out + 0) = res.rmsd;
                            *(loc_dups[loc_ndups].out + 1) = res.dmax;
                            *(loc_dups[loc_ndups].out + 2) = res.amax;
                            loc_dups[loc_ndups].dE_kcal = dE_kcal;
                            loc_ndups++;

                            if (cand == cm) {
                                loc_removed[m] = true;
                            } else {
                                loc_removed[i] = true;
                                break;
                            }
                        }
                    }
                }
                free(loc_removed);
            }

            #pragma omp critical(merge_dups)
            {
                total_tested += loc_tested;
                for (long int k = 0; k < loc_ndups; k++) {
                    if (n_iter_dups >= iter_cap) {
                        iter_cap = (iter_cap == 0) ? 64 : iter_cap * 2;
                        iter_dups = realloc(iter_dups, (size_t)iter_cap * sizeof(DupRecord));
                    }
                    iter_dups[n_iter_dups++] = loc_dups[k];
                }
            }
            free(loc_dups);
        }

        if (n_iter_dups > 1) {
            qsort(iter_dups, (size_t)n_iter_dups, sizeof(DupRecord), cmp_duprec_canon);
        }
        for (long int k = 0; k < n_iter_dups; k++) {
            int keep_orig = iter_dups[k].idx1;
            int dup_orig  = iter_dups[k].idx2;
            if (!confs[keep_orig].removed && !confs[dup_orig].removed) {
                confs[dup_orig].removed = true;
                if (ndups >= dups_cap) {
                    dups_cap = (dups_cap == 0) ? 64 : dups_cap * 2;
                    dups = realloc(dups, (size_t)dups_cap * sizeof(DupRecord));
                }
                dups[ndups++] = iter_dups[k];
            }
        }
        free(iter_dups);

        printf("%ld\n", total_tested);
        printf("LEN: %ld\n", ndups);

        if (ndups == prev_dups && pm >= perm_max) {
            break;
        } else if (ndups == prev_dups && pm < perm_max) {
            int coin = 0;
            if (t_desc > 1) {
                coin = (int)(rand_r(&rng_seed) % (unsigned int)t_desc);
            }
            int step = (128 - chunk_base) / perm_max + coin;
            if (step < 1) step = 1;
            chunk_sz += step;
            G_SORT_KEY = coin;
            G_SORT_REV = (chunk_sz % 2 == 0);
            qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);
            pm++;
        } else {
            long int w = 0;
            for (long int r = 0; r < n_active; r++) {
                if (!active[r]->removed) {
                    active[w++] = active[r];
                }
            }
            n_active = w;
            pm = 0;
            chunk_sz = chunk_base;
        }
    }

    if (pm != perm_max) {
        printf("Warning: increase number of iter=%d\n", iter_max);
    }

    G_SORT_KEY = 0;
    G_SORT_REV = false;
    qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);


    if (pm != perm_max) {
        printf("Warning: increase number of iter=%d\n", iter_max);
    }

    /* =====================================================================
     * Final deterministic exact pass:
     *  - With energies: checks all pairs within dE <= 2.0 kcal/mol window
     *  - Without energies (MIX mode): checks all remaining pairs (no NucRep cutoff)
     * ===================================================================== */
    G_SORT_KEY = 0;
    G_SORT_REV = false;
    qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);

    if (n_active > 1) {
        printf("Final deterministic exact pass on %ld structures...\n", n_active);
        DupRecord *exact_dups = NULL;
        long int n_exact = 0;
        long int exact_cap = 0;
        long int exact_tested = 0;

        #pragma omp parallel
        {
            DupRecord *loc_dups = NULL;
            long int loc_ndups = 0;
            long int loc_cap = 0;
            long int loc_tested = 0;

            #pragma omp for schedule(dynamic, 4)
            for (long int i = 0; i < n_active - 1; i++) {
                ConfItem *ci = active[i];
                double d0_i = *(ci->desc + 0);

                for (long int m = i + 1; m < n_active; m++) {
                    ConfItem *cm = active[m];
                    double d0_m = *(cm->desc + 0);

                    if (!mix_mode) {
                        double dE_kcal = fabs(d0_m - d0_i) * (is_kcal_units ? 1.0 : 627.51);
                        if (dE_kcal > 2.0) break;
                    }

                    long int id_a = (ci->orig_idx < cm->orig_idx) ? ci->orig_idx : cm->orig_idx;
                    long int id_b = (ci->orig_idx < cm->orig_idx) ? cm->orig_idx : ci->orig_idx;

                    bool already = false;
                    #pragma omp critical(pair_cache)
                    {
                        already = pairset_test_and_insert(&tested_pairs, id_a, id_b);
                    }
                    if (already) continue;

                    loc_tested++;
                    McResult res = MC_COMPARE_PRECOMP(pc, ci->xyz, cm->xyz, &cfg.mc_opts);
                    if (res.status == 0 && is_duplicate(&cfg, &res)) {
                        double dE = fabs(d0_m - d0_i);
                        double dE_kcal = mix_mode ? 0.0 : (is_kcal_units ? dE : (dE * 627.51));

                        if (loc_ndups >= loc_cap) {
                            loc_cap = (loc_cap == 0) ? 16 : loc_cap * 2;
                            loc_dups = realloc(loc_dups, (size_t)loc_cap * sizeof(DupRecord));
                        }
                        loc_dups[loc_ndups].idx1 = (int)ci->orig_idx;
                        loc_dups[loc_ndups].idx2 = (int)cm->orig_idx;
                        *(loc_dups[loc_ndups].out + 0) = res.rmsd;
                        *(loc_dups[loc_ndups].out + 1) = res.dmax;
                        *(loc_dups[loc_ndups].out + 2) = res.amax;
                        loc_dups[loc_ndups].dE_kcal = dE_kcal;
                        loc_ndups++;
                    }
                }
            }

            #pragma omp critical(merge_exact)
            {
                exact_tested += loc_tested;
                for (long int k = 0; k < loc_ndups; k++) {
                    if (n_exact >= exact_cap) {
                        exact_cap = (exact_cap == 0) ? 64 : exact_cap * 2;
                        exact_dups = realloc(exact_dups, (size_t)exact_cap * sizeof(DupRecord));
                    }
                    exact_dups[n_exact++] = loc_dups[k];
                }
            }
            free(loc_dups);
        }
/* Deterministic merge in canonical order */
if (n_exact > 1) {
qsort(exact_dups, (size_t)n_exact, sizeof(DupRecord), cmp_duprec_canon);
}
for (long int k = 0; k < n_exact; k++) {
            int keep_orig = exact_dups[k].idx1;
            int dup_orig  = exact_dups[k].idx2;
            if (!confs[keep_orig].removed && !confs[dup_orig].removed) {
                confs[dup_orig].removed = true;
                if (ndups >= dups_cap) {
                    dups_cap = (dups_cap == 0) ? 64 : dups_cap * 2;
                    dups = realloc(dups, (size_t)dups_cap * sizeof(DupRecord));
                }
                dups[ndups++] = exact_dups[k];
            }
        }
        free(exact_dups);

        long int w = 0;
        for (long int r = 0; r < n_active; r++) {
            if (!active[r]->removed) {
                active[w++] = active[r];
            }
        }
        n_active = w;
        total_tested += exact_tested;
        printf("Exact pass tested %ld extra pairs, total DUP = %ld\n", exact_tested, ndups);
    }

    printf("UNIQUE = %ld\n", n_active);



    if (ndups > 1) {
        qsort(dups, (size_t)ndups, sizeof(DupRecord), cmp_duprec);
    }

    FILE *fd = fopen(out_dup_txt, "w");
    FILE *fdx = fopen(out_dupxyz_txt, "w");
    if (fd && fdx) {
        fprintf(fd, "DUP = %ld\n", ndups);
        fprintf(fdx, "DUP = %ld (2nd is the DUP)\n", ndups);
        for (long int i = 0; i < ndups; i++) {
            const ConfItem *c1 = &confs[dups[i].idx1];
            const ConfItem *c2 = &confs[dups[i].idx2];
            fprintf(fdx, "%s %s\n", c1->name, c2->name);
            fprintf(fd, "%s %s %.3f %.3f %.3f %.2f\n",
                    c1->name, c2->name,
                    *(dups[i].out + 0), *(dups[i].out + 1), *(dups[i].out + 2),
                    dups[i].dE_kcal);
        }
    }
    if (fd) fclose(fd);
    if (fdx) fclose(fdx);

    if (!single_file) {
        mkdir(out_unique_dir, 0755);
        for (long int i = 0; i < n_active; i++) {
            char dst[MC_PATH_SZ];
            snprintf(dst, sizeof(dst), "%s/%s", out_unique_dir, active[i]->name);
            copy_file_raw(active[i]->name, dst);
        }
    } else {
        FILE *fo = fopen(out_unique_xyz, "w");
        if (fo) {
            for (long int i = 0; i < n_active; i++) {
                const ConfItem *ci = active[i];
                fprintf(fo, "%ld\n", natoms);
                if (mix_mode) {
                    fprintf(fo, "%s\n", ci->name);
                } else if (is_kcal_units) {
                    fprintf(fo, "Conf %ld  E = %.6f kcal/mol\n", i + 1, *(ci->desc + 0));
                } else {
                    fprintf(fo, "%.10f\n", *(ci->desc + 0));
                }
                for (long int a = 0; a < natoms; a++) {
                    fprintf(fo, "%3s %14.8f %14.8f %14.8f\n",
                            elements[a],
                            *(ci->xyz + 3 * a + 0),
                            *(ci->xyz + 3 * a + 1),
                            *(ci->xyz + 3 * a + 2));
                }
            }
            fclose(fo);
        }
    }
    MC_PRECOMP_FREE(pc);
    for (long int c = 0; c < nconfs; c++) free(confs[c].xyz);
    free(confs);
    free(elements);
    free(Z_arr);
    free(nn_ref);
    free(nbrs_ref);
    free(active);
    pairset_free(&tested_pairs);
    free(dups);

    return 0;
}
