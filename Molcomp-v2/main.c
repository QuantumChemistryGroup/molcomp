#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include "mc_types.h"
#include "mc_io.h"
#include "mc_getconnect.h"
#include "mc_core.h"

static char *mc_strlwr(char *str) {
    unsigned char *p = (unsigned char *)str;
    while (*p) {
        *p = (unsigned char)tolower(*p);
        p++;
    }
    return str;
}

static bool parse_bool_flag(const char *arg, const char *prefix, bool *out_val) {
    size_t len = strlen(prefix);
    if (strncmp(arg, prefix, len) != 0) return false;
    const char *val = arg + len;
    if (strcmp(val, "true") == 0) {
        *out_val = true;
    } else if (strcmp(val, "false") == 0) {
        *out_val = false;
    } else {
        printf("Your %.*s=%s is not clear\n", (int)(len - 1), prefix, val);
        exit(0);
    }
    return true;
}

int main(int argc, char *argv[]) {
    if (argc < 3 || argv[1] == NULL || argv[2] == NULL) {
        printf("Files with coordinates are not specified.\n");
        exit(0);
    }

    McOptions opts = MC_DEFAULT_OPTIONS();
    bool PRINT = false;
    bool WRITETOPO = false;
    bool READTOPO = false;

    for (int i = 3; i < argc; i++) {
        char *a = mc_strlwr(argv[i]);
        if (parse_bool_flag(a, "--verbose=", &opts.VER)) continue;
        if (parse_bool_flag(a, "--print=", &PRINT)) continue;
        if (parse_bool_flag(a, "--keepenantio=", &opts.STEREO)) continue;
        if (parse_bool_flag(a, "--connect=", &opts.CONNECT)) continue;
        if (parse_bool_flag(a, "--writetopo=", &WRITETOPO)) continue;
        if (parse_bool_flag(a, "--readtopo=", &READTOPO)) continue;
        if (parse_bool_flag(a, "--mol2=", &opts.MOL2)) continue;
        if (parse_bool_flag(a, "--bindh2heavy=", &opts.BINDH2HEAVY)) continue;

        if (strncmp(a, "--linearthld=", 13) == 0) {
            opts.LINEARTHLD = atof(a + 13);
            continue;
        }
        if (strncmp(a, "--alg=", 6) == 0) {
            const char *val = a + 6;
            if (strcmp(val, "std") == 0 || strcmp(val, "adv") == 0 || strcmp(val, "kabsch") == 0) {
                strcpy(opts.ALG, val);
            } else {
                printf("Your --alg=%s is not clear\n", val);
                exit(0);
            }
            continue;
        }
        if (strncmp(a, "--compare=", 10) == 0) {
            const char *val = a + 10;
            if (strcmp(val, "angle") == 0) {
                opts.CONNECT = true;
            } else {
                printf("Your --compare=%s is not clear.\n", val);
                exit(0);
            }
            continue;
        }
    }

    char fileout[300];
    int base_len = (int)strlen(argv[1]) - 4;
    if (base_len < 0) base_len = (int)strlen(argv[1]);
    snprintf(fileout, sizeof(fileout), "%.*s-mod.xyz", base_len, argv[1]);

    long int NAT1 = 0;
    long int NAT2 = 0;
    McAtom *coord1 = malloc(1 * sizeof(McAtom));
    McAtom *coord2 = malloc(1 * sizeof(McAtom));

    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
        printf("File 1 with coordinates does not exist.\n");
        exit(0);
    }
    if (opts.VER) printf("Start reading first input file: %s \n", argv[1]);
    MC_FREAD(fp, &coord1, &NAT1, opts.VER);
    fclose(fp);

    fp = fopen(argv[2], "r");
    if (!fp) {
        printf("File 2 with coordinates does not exist.\n");
        exit(0);
    }
    if (opts.VER) printf("Start reading second input file: %s \n", argv[2]);
    MC_FREAD(fp, &coord2, &NAT2, opts.VER);
    fclose(fp);

    if (NAT1 != NAT2) {
        printf("NAT1 is not equal to NAT2. Stop. %li %li \n", NAT1, NAT2);
        exit(0);
    }
    long int NATOMS = NAT1;

    bool skip1 = false;
    bool skip2 = false;

    if (opts.CONNECT == true) {
        if (READTOPO == false) {
            MC_GETCONNECT(&coord1, NATOMS);
        } else {
            fp = fopen("topo.txt", "r");
            if (!fp) {
                printf("FILE WITH TOPOLOGY %s DOES NOT EXIST. STOP.\n", "topo.txt");
                exit(0);
            }
            MC_TXYZREAD(fp, &coord1, NATOMS);
            fclose(fp);
        }
        MC_SYMMETRYCHECK(&coord1, NATOMS);
        skip1 = true;

        if (WRITETOPO == true && opts.MOL2 == false) {
            fp = fopen("topo.txt", "w");
            fprintf(fp, "%3li MOLECULAR CONNECTIVITY \n", NATOMS);
            for (long int i = 0; i < NATOMS; i++) {
                fprintf(fp, "%5li %5s %10.6f %10.6f %10.6f %8s", i + 1, coord1[i].El, coord1[i].xyz[0], coord1[i].xyz[1], coord1[i].xyz[2], "UDF");
                for (long int j = 0; j < coord1[i].nn; j++) {
                    fprintf(fp, "%5li ", coord1[i].n[j]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }

        if (opts.MOL2 == true) {
            if (READTOPO == false) {
                MC_GETCONNECT(&coord2, NATOMS);
            } else {
                fp = fopen("topo.txt", "r");
                if (!fp) {
                    printf("FILE WITH TOPOLOGY %s DOES NOT EXIST. STOP.\n", "topo.txt");
                    exit(0);
                }
                MC_TXYZREAD(fp, &coord2, NATOMS);
                fclose(fp);
            }
            MC_SYMMETRYCHECK(&coord2, NATOMS);
            skip2 = true;

            if (WRITETOPO == true) {
                printf("As --mol2=true connectivity of second molecule is written.\n");
                fp = fopen("topo.txt", "w");
                fprintf(fp, "%3li MOLECULAR CONNECTIVITY \n", NATOMS);
                for (long int i = 0; i < NATOMS; i++) {
                    fprintf(fp, "%5li %5s %10.6f %10.6f %10.6f %8s", i + 1, coord2[i].El, coord2[i].xyz[0], coord2[i].xyz[1], coord2[i].xyz[2], "UDF");
                    for (long int j = 0; j < coord2[i].nn; j++) {
                        fprintf(fp, "%5li ", coord2[i].n[j]);
                    }
                    fprintf(fp, "\n");
                }
                fclose(fp);
            }
        }
    }

    McAtom *coord1best = malloc((size_t)NATOMS * sizeof(McAtom));
    McResult res = MC_COMPARE(coord1, coord2, NATOMS, &opts, skip1, skip2, coord1best);

    if (res.status == -3) {
        printf("Number of unique elements is not the same in both structures. Exit.\n");
        exit(0);
    } else if (res.status == -4) {
        printf("Mismatch in element counts between structures. Close.\n");
        exit(0);
    } else if (res.status != 0) {
        printf("Comparison error: %d\n", res.status);
        exit(0);
    }

    if (PRINT) {
        MC_FPRINT(fileout, coord1best, NATOMS, res.rmsd);
    }
    printf("%f %f %f\n", res.rmsd, res.dmax, res.amax);

    free(coord1);
    free(coord2);
    free(coord1best);
    return 0;
}
