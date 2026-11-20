#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "mc_types.h"
#include "mc_atominfo.h"
#include "mc_inertia.h"
#include "mc_kabshangle.h"
#include "mc_renumb_hung.h"
#include "mc_getconnect.h"
#include "mc_stdfunc.h"
#include "mc_core.h"

McOptions MC_DEFAULT_OPTIONS(void) {
    McOptions opts;
    opts.VER = false;
    opts.STEREO = false;
    opts.CONNECT = true;
    opts.BINDH2HEAVY = true;
    opts.MOL2 = false;
    opts.LINEARTHLD = 5.0;
    strcpy(opts.ALG, "std");
    return opts;
}

McResult MC_COMPARE(
    const McAtom *coord1_in,
    const McAtom *coord2_in,
    long int NATOMS,
    const McOptions *opts,
    bool skip_connect1,
    bool skip_connect2,
    McAtom *coord1best_out
) {
    McResult res = {0.0, 0.0, 0.0, 0};

    if (!coord1_in || !coord2_in || !opts || NATOMS <= 0) {
        res.status = -1;
        return res;
    }

    bool VER = opts->VER;
    bool STEREO = opts->STEREO;
    bool CONNECT = opts->CONNECT;
    bool BINDH2HEAVY = opts->BINDH2HEAVY;
    bool MOL2 = opts->MOL2;
    double LINEARTHLD = opts->LINEARTHLD;

    McAtom *coord1 = malloc((size_t)NATOMS * sizeof(McAtom));
    McAtom *coord2 = malloc((size_t)NATOMS * sizeof(McAtom));
    if (!coord1 || !coord2) {
        free(coord1);
        free(coord2);
        res.status = -2;
        return res;
    }
    memcpy(coord1, coord1_in, (size_t)NATOMS * sizeof(McAtom));
    memcpy(coord2, coord2_in, (size_t)NATOMS * sizeof(McAtom));

    long int (*ADJ1)[ NATOMS ] = malloc((size_t)NATOMS * sizeof(*ADJ1));
    long int (*ADJ2)[ NATOMS ] = malloc((size_t)NATOMS * sizeof(*ADJ2));

    if (CONNECT == true) {
        if (!skip_connect1) {
            MC_GETCONNECT(&coord1, NATOMS);
        }
        MC_SYMMETRYCHECK(&coord1, NATOMS);
        MC_POPADJ(NATOMS, ADJ1, coord1);

        if (MOL2 == true) {
            if (!skip_connect2) {
                MC_GETCONNECT(&coord2, NATOMS);
            }
            MC_SYMMETRYCHECK(&coord2, NATOMS);
            MC_POPADJ(NATOMS, ADJ2, coord2);
        }
    }

    /* Populate Z for both molecules */
    for (long int i = 0; i != NATOMS; i++) {
        coord1[ i ].Z = (double)MC_ATOMINFO(coord1[ i ].El).Z;
        coord2[ i ].Z = (double)MC_ATOMINFO(coord2[ i ].El).Z;
    }

    /* Identity matrix of atoms p and q */
    bool (*IDN)[ NATOMS ] = malloc((size_t)NATOMS * sizeof(*IDN));
    if (CONNECT == true) {
        for (long int i = 0; i != NATOMS; i++) {
            for (long int j = 0; j != NATOMS; j++) {
                bool EqEl = false;
                if (MOL2 == true) {
                    if (coord1[ i ].Z == coord2[ j ].Z) EqEl = true;
                } else {
                    if (coord1[ i ].Z == coord1[ j ].Z) EqEl = true;
                }
                if (EqEl == false) {
                    IDN[ i ][ j ] = false;
                } else {
                    IDN[ i ][ j ] = MC_RADIAL(i, j, coord1, coord2, NATOMS, MOL2, ADJ1, ADJ2);
                }
            }
        }
    }

    McAtom *coord1mod  = malloc((size_t)NATOMS * sizeof(McAtom));
    McAtom *coord1modk = malloc((size_t)NATOMS * sizeof(McAtom));
    McAtom *coord1best = malloc((size_t)NATOMS * sizeof(McAtom));
    memcpy(coord1best, coord1, (size_t)NATOMS * sizeof(McAtom));

    /* Principal moments of inertia and eigenvectors */
    double IP1[ 3 ] = {0, 0, 0};
    double IP2[ 3 ] = {0, 0, 0};
    double EF1[ 3 ][ 3 ] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double EF2[ 3 ][ 3 ] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    MC_INERTIA(&coord1, IP1, NATOMS, EF1, VER);
    MC_INERTIA(&coord2, IP2, NATOMS, EF2, VER);
    (void)IP1;
    (void)IP2;

    /* Count unique elements */
    char (*UnEl)[ 4 ]  = malloc(1 * sizeof(*UnEl));
    char (*UnEl0)[ 4 ] = malloc(1 * sizeof(*UnEl0));
    long int UE = 0;
    long int UE0 = 0;
    long int *ElNum  = malloc(1 * sizeof(long int));
    long int *ElNum0 = malloc(1 * sizeof(long int));

    for (long int i = 0; i != NATOMS; i++) {
        if (i == 0) {
            strcpy(UnEl[ 0 ], coord2[ i ].El);
            strcpy(UnEl0[ 0 ], coord1[ i ].El);
            UE = 1;
            UE0 = 1;
            UnEl  = realloc(UnEl,  (size_t)(UE + 1)  * sizeof(*UnEl));
            UnEl0 = realloc(UnEl0, (size_t)(UE0 + 1) * sizeof(*UnEl0));
            ElNum  = realloc(ElNum,  (size_t)(UE + 1)  * sizeof(long int));
            ElNum0 = realloc(ElNum0, (size_t)(UE0 + 1) * sizeof(long int));
            ElNum[ 0 ] = 1;
            ElNum0[ 0 ] = 1;
        } else {
            bool MATCH = false;
            for (long int j = 0; j != UE; j++) {
                if (strcmp(coord2[ i ].El, UnEl[ j ]) == 0) {
                    MATCH = true;
                    ElNum[ j ]++;
                    break;
                }
            }
            if (MATCH == false) {
                strcpy(UnEl[ UE ], coord2[ i ].El);
                UE++;
                UnEl  = realloc(UnEl,  (size_t)(UE + 1) * sizeof(*UnEl));
                ElNum = realloc(ElNum, (size_t)(UE + 1) * sizeof(long int));
                ElNum[ UE - 1 ] = 1;
            }

            MATCH = false;
            for (long int j = 0; j != UE0; j++) {
                if (strcmp(coord1[ i ].El, UnEl0[ j ]) == 0) {
                    MATCH = true;
                    ElNum0[ j ]++;
                    break;
                }
            }
            if (MATCH == false) {
                strcpy(UnEl0[ UE0 ], coord1[ i ].El);
                UE0++;
                UnEl0  = realloc(UnEl0,  (size_t)(UE0 + 1) * sizeof(*UnEl0));
                ElNum0 = realloc(ElNum0, (size_t)(UE0 + 1) * sizeof(long int));
                ElNum0[ UE0 - 1 ] = 1;
            }
        }
    }

    if (UE != UE0) {
        res.status = -3;
        goto cleanup;
    }
for (long int i = 0; i != UE; i++) {
bool found_el = false;
for (long int j = 0; j != UE0; j++) {
if (strcmp(UnEl[ i ], UnEl0[ j ]) == 0) {
found_el = true;
if (ElNum[ i ] != ElNum0[ j ]) {
res.status = -4;
goto cleanup;
}
}
}
if (!found_el) {
res.status = -4;
goto cleanup;
}
}

    /* Orient along principal axes */
    for (long int i = 0; i != NATOMS; i++) {
        double xyz[ 3 ] = {coord1[ i ].xyz[ 0 ], coord1[ i ].xyz[ 1 ], coord1[ i ].xyz[ 2 ]};
        coord1[ i ].xyz[ 0 ] = xyz[ 0 ] * EF1[ 0 ][ 0 ] + xyz[ 1 ] * EF1[ 1 ][ 0 ] + xyz[ 2 ] * EF1[ 2 ][ 0 ];
        coord1[ i ].xyz[ 1 ] = xyz[ 0 ] * EF1[ 0 ][ 1 ] + xyz[ 1 ] * EF1[ 1 ][ 1 ] + xyz[ 2 ] * EF1[ 2 ][ 1 ];
        coord1[ i ].xyz[ 2 ] = xyz[ 0 ] * EF1[ 0 ][ 2 ] + xyz[ 1 ] * EF1[ 1 ][ 2 ] + xyz[ 2 ] * EF1[ 2 ][ 2 ];
    }

    for (long int i = 0; i != NATOMS; i++) {
        double xyz[ 3 ] = {coord2[ i ].xyz[ 0 ], coord2[ i ].xyz[ 1 ], coord2[ i ].xyz[ 2 ]};
        coord2[ i ].xyz[ 0 ] = xyz[ 0 ] * EF2[ 0 ][ 0 ] + xyz[ 1 ] * EF2[ 1 ][ 0 ] + xyz[ 2 ] * EF2[ 2 ][ 0 ];
        coord2[ i ].xyz[ 1 ] = xyz[ 0 ] * EF2[ 0 ][ 1 ] + xyz[ 1 ] * EF2[ 1 ][ 1 ] + xyz[ 2 ] * EF2[ 2 ][ 1 ];
        coord2[ i ].xyz[ 2 ] = xyz[ 0 ] * EF2[ 0 ][ 2 ] + xyz[ 1 ] * EF2[ 1 ][ 2 ] + xyz[ 2 ] * EF2[ 2 ][ 2 ];
    }

    double Imag[ 8 ][ 3 ] = {
        { 1,  1,  1},
        { 1, -1,  1},
        {-1, -1,  1},
        {-1,  1,  1},
        {-1,  1, -1},
        { 1,  1, -1},
        { 1, -1, -1},
        {-1, -1, -1}
    };

    double RMSDmin = 1e10;

    if (strcmp(opts->ALG, "std") == 0) {
        McAtom *coord1tmp = calloc((size_t)NATOMS, sizeof(McAtom));
        for (long int i = 0; i != 8; i++) {
            if (STEREO == true) {
                long int REFL = (long int)(Imag[ i ][ 0 ] * Imag[ i ][ 1 ] * Imag[ i ][ 2 ]);
                if (REFL == -1) continue;
            }
            for (long int j = 0; j != NATOMS; j++) {
                strcpy(coord1tmp[ j ].El, coord1[ j ].El);
                for (long int m = 0; m != 3; m++) {
                    coord1tmp[ j ].xyz[ m ] = coord1[ j ].xyz[ m ] * Imag[ i ][ m ];
                }
                coord1tmp[ j ].Z = coord1[ j ].Z;
                if (CONNECT == true) {
                    coord1tmp[ j ].nn = coord1[ j ].nn;
                    for (long int m = 0; m != coord1tmp[ j ].nn; m++) {
                        coord1tmp[ j ].n[ m ] = coord1[ j ].n[ m ];
                    }
                }
            }

            if (CONNECT != true) {
                MC_RENUMBHUNG(&coord1tmp, coord2, NATOMS, UnEl, UE);
            } else {
                if (BINDH2HEAVY == true) {
                    MC_RENUMBHUNGCONNECT(&coord1tmp, coord2, NATOMS, UnEl, UE, IDN, MOL2);
                } else {
                    MC_RENUMBHUNGGRAPH(&coord1tmp, coord2, NATOMS, UnEl, UE, IDN);
                }
            }

            double RMSD = 0, DAV_tmp = 0, AMAX_tmp = 0;
            MC_KABSHANGLE(coord1tmp, coord2, &coord1mod, NATOMS, &RMSD, &DAV_tmp, &AMAX_tmp, MOL2, LINEARTHLD, VER);
            if (RMSD < RMSDmin) {
                RMSDmin = RMSD;
                for (long int j = 0; j != NATOMS; j++) {
                    coord1best[ j ] = coord1tmp[ j ];
                }
            }
        }
        free(coord1tmp);

        if (VER) printf("Standard Kabsch, permutation.\n");
        double RMSD = 0, AAV = 0, AMAX = 0;
        double RMSDST = 0, AAVST = 0, AMAXST = 0;

        MC_KABSHANGLE(coord1best, coord2, &coord1mod, NATOMS, &RMSD, &AAV, &AMAX, MOL2, LINEARTHLD, VER);
        if (VER) printf("Kabsch (renumbering) %f %f %f\n", RMSD, AAV, AMAX);

        MC_KABSHANGLE(coord1, coord2, &coord1modk, NATOMS, &RMSDST, &AAVST, &AMAXST, MOL2, LINEARTHLD, VER);
        if (VER) printf("Kabsch (no renumbering) %f %f %f\n", RMSDST, AAVST, AMAXST);

        if (RMSD < RMSDST) {
            res.rmsd = RMSD;
            res.dmax = AAV;
            res.amax = AMAX;
            if (coord1best_out) memcpy(coord1best_out, coord1mod, (size_t)NATOMS * sizeof(McAtom));
        } else {
            res.rmsd = RMSDST;
            res.dmax = AAVST;
res.amax = AMAXST;
if (coord1best_out) memcpy(coord1best_out, coord1modk, (size_t)NATOMS * sizeof(McAtom));
}
} else if (strcmp(opts->ALG, "adv") == 0) {
long int PERM[ 6 ][ 3 ] = {
{0, 1, 2},
{0, 2, 1},
{1, 0, 2},
{1, 2, 0},
{2, 0, 1},
{2, 1, 0}
};
McAtom *coord1tmp = calloc((size_t)NATOMS, sizeof(McAtom));
for (long int p = 0; p != 6; p++) {
for (long int i = 0; i != 8; i++) {
if (STEREO == true) {
long int perm_parity = (p == 0 || p == 3 || p == 4) ? 1 : -1;
long int REFL = (long int)(Imag[ i ][ 0 ] * Imag[ i ][ 1 ] * Imag[ i ][ 2 ]) * perm_parity;
if (REFL == -1) continue;
}
for (long int j = 0; j != NATOMS; j++) {
strcpy(coord1tmp[ j ].El, coord1[ j ].El);
for (long int m = 0; m != 3; m++) {
long int q = PERM[ p ][ m ];
coord1tmp[ j ].xyz[ m ] = coord1[ j ].xyz[ q ] * Imag[ i ][ m ];
}
coord1tmp[ j ].Z = coord1[ j ].Z;
if (CONNECT == true) {
coord1tmp[ j ].nn = coord1[ j ].nn;
for (long int m = 0; m != coord1tmp[ j ].nn; m++) {
coord1tmp[ j ].n[ m ] = coord1[ j ].n[ m ];
}
}
}
if (CONNECT != true) {
MC_RENUMBHUNG(&coord1tmp, coord2, NATOMS, UnEl, UE);
} else {
if (BINDH2HEAVY == true) {
MC_RENUMBHUNGCONNECT(&coord1tmp, coord2, NATOMS, UnEl, UE, IDN, MOL2);
} else {
MC_RENUMBHUNGGRAPH(&coord1tmp, coord2, NATOMS, UnEl, UE, IDN);
}
}
double RMSD = 0, DAV_tmp = 0, AMAX_tmp = 0;
MC_KABSHANGLE(coord1tmp, coord2, &coord1mod, NATOMS, &RMSD, &DAV_tmp, &AMAX_tmp, MOL2, LINEARTHLD, VER);
if (RMSD < RMSDmin) {
RMSDmin = RMSD;
for (long int j = 0; j != NATOMS; j++) {
coord1best[ j ] = coord1tmp[ j ];
}
}
}
}
free(coord1tmp);
if (VER) printf("Advanced Kabsch, permutation.\n");
double RMSD = 0, AAV = 0, AMAX = 0;
double RMSDST = 0, AAVST = 0, AMAXST = 0;
MC_KABSHANGLE(coord1best, coord2, &coord1mod, NATOMS, &RMSD, &AAV, &AMAX, MOL2, LINEARTHLD, VER);
if (VER) printf("Kabsch (renumbering) %f %f %f\n", RMSD, AAV, AMAX);
MC_KABSHANGLE(coord1, coord2, &coord1modk, NATOMS, &RMSDST, &AAVST, &AMAXST, MOL2, LINEARTHLD, VER);
if (VER) printf("Kabsch (no renumbering) %f %f %f\n", RMSDST, AAVST, AMAXST);
if (RMSD < RMSDST) {
res.rmsd = RMSD;
res.dmax = AAV;
res.amax = AMAX;
if (coord1best_out) memcpy(coord1best_out, coord1mod, (size_t)NATOMS * sizeof(McAtom));
} else {
res.rmsd = RMSDST;
res.dmax = AAVST;
res.amax = AMAXST;
if (coord1best_out) memcpy(coord1best_out, coord1modk, (size_t)NATOMS * sizeof(McAtom));
}
} else if (strcmp(opts->ALG, "kabsch") == 0) {
if (VER) printf("Standard Kabsch, no permutation.\n");
double RMSD = 0, AAV = 0, AMAX = 0;
MC_KABSHANGLE(coord1, coord2, &coord1mod, NATOMS, &RMSD, &AAV, &AMAX, MOL2, LINEARTHLD, VER);
res.rmsd = RMSD;
res.dmax = AAV;
res.amax = AMAX;
if (coord1best_out) memcpy(coord1best_out, coord1mod, (size_t)NATOMS * sizeof(McAtom));
} else {
res.status = -5;
}
cleanup:
free(coord1);
free(coord2);
free(ADJ1);
free(ADJ2);
free(IDN);
free(coord1mod);
free(coord1modk);
free(coord1best);
free(UnEl);
free(UnEl0);
free(ElNum);
free(ElNum0);
return res;
}
McResult MC_COMPARE_FLAT(
long int NATOMS,
const char (*elements)[ 4 ],
const long int *nn,
const long int (*nbrs)[ 23 ],
const double *xyz1,
const double *xyz2,
const McOptions *opts
) {
McResult res = {0.0, 0.0, 0.0, 0};
if (NATOMS <= 0 || !elements || !xyz1 || !xyz2 || !opts) {
res.status = -1;
return res;
}
McAtom *c1 = calloc((size_t)NATOMS, sizeof(McAtom));
McAtom *c2 = calloc((size_t)NATOMS, sizeof(McAtom));
if (!c1 || !c2) {
free(c1);
free(c2);
res.status = -2;
return res;
}
bool has_conn = (nn != NULL && nbrs != NULL);
for (long int i = 0; i != NATOMS; i++) {
strncpy(c1[ i ].El, elements[ i ], 3);
c1[ i ].El[ 3 ] = '\0';
c1[ i ].xyz[ 0 ] = xyz1[ 3 * i + 0 ];
c1[ i ].xyz[ 1 ] = xyz1[ 3 * i + 1 ];
c1[ i ].xyz[ 2 ] = xyz1[ 3 * i + 2 ];
strncpy(c2[ i ].El, elements[ i ], 3);
c2[ i ].El[ 3 ] = '\0';
c2[ i ].xyz[ 0 ] = xyz2[ 3 * i + 0 ];
c2[ i ].xyz[ 1 ] = xyz2[ 3 * i + 1 ];
c2[ i ].xyz[ 2 ] = xyz2[ 3 * i + 2 ];
if (has_conn) {
long int nn_i = nn[ i ];
if (nn_i < 0) nn_i = 0;
if (nn_i > 23) nn_i = 23;
long int accepted = 0;
for (long int k = 0; k < nn_i; k++) {
long int nb = nbrs[ i ][ k ];
if (nb >= 1 && nb <= NATOMS) {
c1[ i ].n[ accepted ] = nb;
c2[ i ].n[ accepted ] = nb;
accepted++;
}
}
c1[ i ].nn = accepted;
c2[ i ].nn = accepted;
}
}
res = MC_COMPARE(c1, c2, NATOMS, opts, has_conn, has_conn, NULL);
free(c1);
free(c2);
return res;
}

struct McPrecomp {
    long int natoms;
    char (*elements)[ 4 ];
    double *Z;
    long int *nn;
    long int (*nbrs)[ 23 ];
    bool *idn;              /* flat [natoms * natoms] */
    char (*un_el)[ 4 ];
    long int ue;
    McMM *tors;
    long int ntors;
};

McPrecomp *MC_PRECOMP_CREATE(
    long int NATOMS,
    const char (*elements)[ 4 ],
    const double *xyz_ref,
    const long int *nn_in,
    const long int (*nbrs_in)[ 23 ],
    const McOptions *opts
) {
    if (NATOMS <= 0 || !elements || !opts) return NULL;

    McPrecomp *pc = calloc(1, sizeof(McPrecomp));
    if (!pc) return NULL;
    pc->natoms = NATOMS;
    pc->elements = calloc((size_t)NATOMS, sizeof(*pc->elements));
    pc->Z = calloc((size_t)NATOMS, sizeof(double));
    pc->nn = calloc((size_t)NATOMS, sizeof(long int));
    pc->nbrs = calloc((size_t)NATOMS, sizeof(*pc->nbrs));

    McAtom *ref = calloc((size_t)NATOMS, sizeof(McAtom));
    for (long int i = 0; i < NATOMS; i++) {
        strncpy(pc->elements[ i ], elements[ i ], 3);
        pc->elements[ i ][ 3 ] = '\0';
        pc->Z[ i ] = (double)MC_ATOMINFO(pc->elements[ i ]).Z;

        strcpy(ref[ i ].El, pc->elements[ i ]);
        ref[ i ].Z = pc->Z[ i ];
        if (xyz_ref) {
            *(ref[ i ].xyz + 0) = *(xyz_ref + 3 * i + 0);
            *(ref[ i ].xyz + 1) = *(xyz_ref + 3 * i + 1);
            *(ref[ i ].xyz + 2) = *(xyz_ref + 3 * i + 2);
        }
        if (nn_in && nbrs_in) {
        long int nn_i = nn_in[ i ];
        if (nn_i < 0) nn_i = 0;
        if (nn_i > 23) nn_i = 23;
        long int accepted = 0;
        for (long int k = 0; k < nn_i; k++) {
        long int nb = nbrs_in[ i ][ k ];
        if (nb >= 1 && nb <= NATOMS) {
        ref[ i ].n[ accepted ] = nb;
        accepted++;
        }
        }
        ref[ i ].nn = accepted;
        }
    }

    if (opts->CONNECT) {
        if (!nn_in || !nbrs_in) {
            MC_GETCONNECT(&ref, NATOMS);
        }
        MC_SYMMETRYCHECK(&ref, NATOMS);

        for (long int i = 0; i < NATOMS; i++) {
            pc->nn[ i ] = ref[ i ].nn;
            for (long int k = 0; k < ref[ i ].nn; k++) {
                pc->nbrs[ i ][ k ] = ref[ i ].n[ k ];
            }
        }

        long int (*ADJ)[ NATOMS ] = malloc((size_t)NATOMS * sizeof(*ADJ));
        MC_POPADJ(NATOMS, ADJ, ref);

        pc->idn = calloc((size_t)NATOMS * (size_t)NATOMS, sizeof(bool));
        for (long int i = 0; i < NATOMS; i++) {
            for (long int j = 0; j < NATOMS; j++) {
                if (pc->Z[ i ] == pc->Z[ j ]) {
                    pc->idn[ i * NATOMS + j ] = MC_RADIAL(i, j, ref, ref, NATOMS, false, ADJ, ADJ);
                }
            }
        }
        free(ADJ);

        /* Precompute torsions once */
        long int cap = 64;
        pc->tors = malloc((size_t)cap * sizeof(McMM));
        long int l = 0;
        for (long int i = 0; i < NATOMS; i++) {
            for (long int j = 0; j < ref[ i ].nn; j++) {
                long int z = ref[ i ].n[ j ];
                for (long int m = 0; m < ref[ i ].nn; m++) {
                    for (long int q = 0; q < ref[ z - 1 ].nn; q++) {
                        if (ref[ i ].n[ m ] == z || ref[ z - 1 ].n[ q ] == i + 1) continue;
                        if (l >= cap) {
                            cap *= 2;
                            pc->tors = realloc(pc->tors, (size_t)cap * sizeof(McMM));
                        }
                        *(pc->tors[ l ].A + 0) = ref[ i ].n[ m ];
                        *(pc->tors[ l ].A + 1) = i + 1;
                        *(pc->tors[ l ].A + 2) = z;
                        *(pc->tors[ l ].A + 3) = ref[ z - 1 ].n[ q ];
                        l++;
                    }
                }
            }
        }
        pc->ntors = l;
    }

    /* Unique elements list */
    pc->un_el = malloc((size_t)(NATOMS + 1) * sizeof(*pc->un_el));
    pc->ue = 0;
    for (long int i = 0; i < NATOMS; i++) {
        bool found = false;
        for (long int j = 0; j < pc->ue; j++) {
            if (strcmp(pc->elements[ i ], pc->un_el[ j ]) == 0) {
                found = true;
                break;
            }
        }
        if (!found) {
            strcpy(pc->un_el[ pc->ue ], pc->elements[ i ]);
            pc->ue++;
        }
    }

    free(ref);
    return pc;
}

void MC_PRECOMP_FREE(McPrecomp *pc) {
    if (!pc) return;
    free(pc->elements);
    free(pc->Z);
    free(pc->nn);
    free(pc->nbrs);
    free(pc->idn);
    free(pc->un_el);
    free(pc->tors);
    free(pc);
}

McResult MC_COMPARE_PRECOMP(
    const McPrecomp *pc,
    const double *xyz1_in,
    const double *xyz2_in,
    const McOptions *opts
) {
    McResult res = {0.0, 0.0, 0.0, 0};
    if (!pc || !xyz1_in || !xyz2_in || !opts) {
        res.status = -1;
        return res;
    }

    long int NATOMS = pc->natoms;
    bool VER = opts->VER;
    bool STEREO = opts->STEREO;
    bool CONNECT = opts->CONNECT;
    bool BINDH2HEAVY = opts->BINDH2HEAVY;
    double LINEARTHLD = opts->LINEARTHLD;

    McAtom *coord1 = calloc((size_t)NATOMS, sizeof(McAtom));
    McAtom *coord2 = calloc((size_t)NATOMS, sizeof(McAtom));
    for (long int i = 0; i < NATOMS; i++) {
        strcpy(coord1[ i ].El, pc->elements[ i ]);
        strcpy(coord2[ i ].El, pc->elements[ i ]);
        coord1[ i ].Z = pc->Z[ i ];
        coord2[ i ].Z = pc->Z[ i ];
        *(coord1[ i ].xyz + 0) = *(xyz1_in + 3 * i + 0);
        *(coord1[ i ].xyz + 1) = *(xyz1_in + 3 * i + 1);
        *(coord1[ i ].xyz + 2) = *(xyz1_in + 3 * i + 2);
        *(coord2[ i ].xyz + 0) = *(xyz2_in + 3 * i + 0);
        *(coord2[ i ].xyz + 1) = *(xyz2_in + 3 * i + 1);
        *(coord2[ i ].xyz + 2) = *(xyz2_in + 3 * i + 2);
        if (CONNECT) {
            coord1[ i ].nn = pc->nn[ i ];
            coord2[ i ].nn = pc->nn[ i ];
            for (long int k = 0; k < pc->nn[ i ]; k++) {
                coord1[ i ].n[ k ] = pc->nbrs[ i ][ k ];
                coord2[ i ].n[ k ] = pc->nbrs[ i ][ k ];
            }
        }
    }

    double IP1[ 3 ] = {0, 0, 0};
    double IP2[ 3 ] = {0, 0, 0};
    double EF1[ 3 ][ 3 ] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double EF2[ 3 ][ 3 ] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    MC_INERTIA(&coord1, IP1, NATOMS, EF1, VER);
    MC_INERTIA(&coord2, IP2, NATOMS, EF2, VER);
    (void)IP1; (void)IP2;

    double *xyz1_pa = malloc((size_t)(3 * NATOMS) * sizeof(double));
    double *xyz2_pa = malloc((size_t)(3 * NATOMS) * sizeof(double));
    for (long int i = 0; i < NATOMS; i++) {
        double x = *(coord1[ i ].xyz + 0);
        double y = *(coord1[ i ].xyz + 1);
        double z = *(coord1[ i ].xyz + 2);
        *(xyz1_pa + 3 * i + 0) = x * EF1[ 0 ][ 0 ] + y * EF1[ 1 ][ 0 ] + z * EF1[ 2 ][ 0 ];
        *(xyz1_pa + 3 * i + 1) = x * EF1[ 0 ][ 1 ] + y * EF1[ 1 ][ 1 ] + z * EF1[ 2 ][ 1 ];
        *(xyz1_pa + 3 * i + 2) = x * EF1[ 0 ][ 2 ] + y * EF1[ 1 ][ 2 ] + z * EF1[ 2 ][ 2 ];

        x = *(coord2[ i ].xyz + 0);
        y = *(coord2[ i ].xyz + 1);
        z = *(coord2[ i ].xyz + 2);
        *(xyz2_pa + 3 * i + 0) = x * EF2[ 0 ][ 0 ] + y * EF2[ 1 ][ 0 ] + z * EF2[ 2 ][ 0 ];
        *(xyz2_pa + 3 * i + 1) = x * EF2[ 0 ][ 1 ] + y * EF2[ 1 ][ 1 ] + z * EF2[ 2 ][ 1 ];
        *(xyz2_pa + 3 * i + 2) = x * EF2[ 0 ][ 2 ] + y * EF2[ 1 ][ 2 ] + z * EF2[ 2 ][ 2 ];
    }
    /* Copy rotated xyz2_pa back into coord2 for correct Hungarian matching */
    for (long int i = 0; i < NATOMS; i++) {
        *(coord2[ i ].xyz + 0) = *(xyz2_pa + 3 * i + 0);
        *(coord2[ i ].xyz + 1) = *(xyz2_pa + 3 * i + 1);
        *(coord2[ i ].xyz + 2) = *(xyz2_pa + 3 * i + 2);
    }

    if (strcmp(opts->ALG, "kabsch") == 0) {
        MC_KABSHANGLE_FAST(xyz1_pa, xyz2_pa, NULL, NATOMS, pc->tors, pc->ntors, true,
                           &res.rmsd, &res.dmax, &res.amax, LINEARTHLD, VER);
        free(coord1); free(coord2); free(xyz1_pa); free(xyz2_pa);
        return res;
    }

    bool (*IDN)[ NATOMS ] = (bool (*)[ NATOMS ])pc->idn;
    double Imag[ 8 ][ 3 ] = {
        { 1,  1,  1}, { 1, -1,  1}, {-1, -1,  1}, {-1,  1,  1},
        {-1,  1, -1}, { 1,  1, -1}, { 1, -1, -1}, {-1, -1, -1}
    };
    long int PERM[ 6 ][ 3 ] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };
    int n_perm = (strcmp(opts->ALG, "adv") == 0) ? 6 : 1;
McAtom *coord1tmp = calloc((size_t)NATOMS, sizeof(McAtom));
double *xyz_tmp   = malloc((size_t)(3 * NATOMS) * sizeof(double));
double *xyz_best  = malloc((size_t)(3 * NATOMS) * sizeof(double));
long int *best_old_of_new = malloc((size_t)NATOMS * sizeof(long int));
long int *cur_old_of_new = malloc((size_t)NATOMS * sizeof(long int));
if (!xyz_tmp || !xyz_best || !best_old_of_new || !cur_old_of_new) exit(1);
for (long int j = 0; j < NATOMS; j++) best_old_of_new[j] = j;
memcpy(xyz_best, xyz1_pa, (size_t)(3 * NATOMS) * sizeof(double));

    double RMSDmin = 1e10;
    for (int p = 0; p < n_perm; p++) {
    for (int i = 0; i < 8; i++) {
    if (STEREO) {
    long int perm_parity = (p == 0 || p == 3 || p == 4) ? 1 : -1;
    long int REFL = (long int)(Imag[ i ][ 0 ] * Imag[ i ][ 1 ] * Imag[ i ][ 2 ]) * perm_parity;
    if (REFL == -1) continue;
    }
            for (long int j = 0; j < NATOMS; j++) {
                strcpy(coord1tmp[ j ].El, pc->elements[ j ]);
                for (int m = 0; m < 3; m++) {
                    long int q = PERM[ p ][ m ];
                    *(coord1tmp[ j ].xyz + m) = *(xyz1_pa + 3 * j + q) * Imag[ i ][ m ];
                }
                coord1tmp[ j ].Z = pc->Z[ j ];
                if (CONNECT) {
                    coord1tmp[ j ].nn = pc->nn[ j ];
                    for (long int m = 0; m < pc->nn[ j ]; m++) {
                        coord1tmp[ j ].n[ m ] = pc->nbrs[ j ][ m ];
                    }
                }
            }

            if (!CONNECT) {
                MC_RENUMBHUNG(&coord1tmp, coord2, NATOMS, pc->un_el, pc->ue);
            } else if (BINDH2HEAVY) {
                MC_RENUMBHUNGCONNECT(&coord1tmp, coord2, NATOMS, pc->un_el, pc->ue, IDN, false);
            } else {
            MC_RENUMBHUNGGRAPH(&coord1tmp, coord2, NATOMS, pc->un_el, pc->ue, IDN);
            }
            for (long int j = 0; j < NATOMS; j++) {
            cur_old_of_new[j] = coord1tmp[j].N;
            if (cur_old_of_new[j] < 0 || cur_old_of_new[j] >= NATOMS) cur_old_of_new[j] = j;
            }
            for (long int j = 0; j < NATOMS; j++) {
            *(xyz_tmp + 3 * j + 0) = *(coord1tmp[ j ].xyz + 0);
                *(xyz_tmp + 3 * j + 1) = *(coord1tmp[ j ].xyz + 1);
                *(xyz_tmp + 3 * j + 2) = *(coord1tmp[ j ].xyz + 2);
            }

            double RMSD = 0.0;
            /* compute_angles = false inside orientation search! */
            MC_KABSHANGLE_FAST(xyz_tmp, xyz2_pa, NULL, NATOMS, NULL, 0, false,
                               &RMSD, NULL, NULL, LINEARTHLD, false);
            if (RMSD < RMSDmin) {
            RMSDmin = RMSD;
            memcpy(xyz_best, xyz_tmp, (size_t)(3 * NATOMS) * sizeof(double));
            memcpy(best_old_of_new, cur_old_of_new, (size_t)NATOMS * sizeof(long int));
            }
}
}
double RMSD = 0, DMAX = 0, AMAX = 0;
double RMSDST = 0, DMAXST = 0, AMAXST = 0;
long int *new_of_old = malloc((size_t)NATOMS * sizeof(long int));
if (!new_of_old) exit(1);
for (long int j = 0; j < NATOMS; j++) new_of_old[j] = j;
for (long int j = 0; j < NATOMS; j++) {
long int old_idx = best_old_of_new[j];
if (old_idx >= 0 && old_idx < NATOMS) new_of_old[old_idx] = j;
}
McMM *tors_best = NULL;
if (pc->ntors > 0) {
tors_best = malloc((size_t)pc->ntors * sizeof(McMM));
if (!tors_best) exit(1);
for (long int t = 0; t < pc->ntors; t++) {
for (int a = 0; a < 4; a++) {
long int old_a = pc->tors[t].A[a] - 1;
if (old_a >= 0 && old_a < NATOMS) {
tors_best[t].A[a] = new_of_old[old_a] + 1;
} else {
tors_best[t].A[a] = pc->tors[t].A[a];
}
}
tors_best[t].V = pc->tors[t].V;
}
}
MC_KABSHANGLE_FAST(xyz_best, xyz2_pa, NULL, NATOMS, tors_best ? tors_best : pc->tors, pc->ntors, true,
&RMSD, &DMAX, &AMAX, LINEARTHLD, VER);
MC_KABSHANGLE_FAST(xyz1_pa, xyz2_pa, NULL, NATOMS, pc->tors, pc->ntors, true,
&RMSDST, &DMAXST, &AMAXST, LINEARTHLD, VER);
free(tors_best);
free(new_of_old);
if (RMSD < RMSDST) {
res.rmsd = RMSD;
res.dmax = DMAX;
res.amax = AMAX;
} else {
res.rmsd = RMSDST;
res.dmax = DMAXST;
res.amax = AMAXST;
}
free(coord1);
free(coord2);
free(xyz1_pa);
free(xyz2_pa);
free(coord1tmp);
free(xyz_tmp);
free(xyz_best);
free(best_old_of_new);
free(cur_old_of_new);
return res;
}
