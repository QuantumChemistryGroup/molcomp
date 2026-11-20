#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "mc_types.h"
#include "mc_hungarian.h"
#include "mc_renumb_hung.h"

void MC_RENUMBHUNGCONNECT(McAtom **coord1tmp, McAtom *coord2, long int NATOMS, char (*UnEl)[4], long int UE, bool IDN[][NATOMS], bool MOL2) {
    long int i = 0;
    long int j = 0;
    // temporary array for renumbering
    McAtom *coord = malloc(NATOMS * sizeof(McAtom));

    long int UnElN[UE];
    for (i = 0; i < UE; i++) {
        UnElN[i] = 0;
    }
    for (i = 0; i < NATOMS; i++) {
        (*coord1tmp)[i].N = -1;
        for (j = 0; j < UE; j++) {
            if (strcmp(coord2[i].El, UnEl[j]) == 0) {
                UnElN[j] = UnElN[j] + 1;
            }
        }
        coord2[i].N = -1;
    }

    // check if any unique element is polyvalent
    bool POLYVAL[UE];
    for (i = 0; i < UE; i++) {
        POLYVAL[i] = false;
        for (j = 0; j < NATOMS; j++) {
            if (strcmp(UnEl[i], (*coord1tmp)[j].El) == 0) {
                if ((*coord1tmp)[j].nn != 1) {
                    POLYVAL[i] = true;
                    break;
                }
            }
        }
    }

    // for each element construct the cost matrix (polyvalent atoms first)
    for (i = 0; i < UE; i++) {
        if (POLYVAL[i] != false) {
            long long int IND = UnElN[i];
            long long int (*COST)[IND] = malloc(IND * sizeof(*COST));
            long int p = 0, q = 0, m = 0, l = 0;

            for (p = 0; p < NATOMS; p++) {
                if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                    for (q = 0; q < NATOMS; q++) {
                        if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                            double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
                            double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
                            double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
                            long long int d = (long long int)(sqrt(x * x + y * y + z * z) * 1e6);
                            if (IDN[p][q] == true) {
                                d = d / 1000;
                            }
                            COST[m][l] = d;
                            l = l + 1;
                        }
                    }
                    l = 0;
                    m = m + 1;
                }
            }

            long int *MATCH = malloc(IND * sizeof(long int));
            MC_HUNGMAIN(IND, IND, COST, &MATCH);

            m = 0;
            l = 0;
            for (p = 0; p < NATOMS; p++) {
                if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                    for (q = 0; q < NATOMS; q++) {
                        if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                            if (MATCH[m] == l) {
                                (*coord1tmp)[p].N = q;
                            }
                            l = l + 1;
                        }
                    }
                    l = 0;
                    m = m + 1;
                }
            }
            free(MATCH);
            free(COST);
        }
    }

    // deal with univalent atoms
    for (i = 0; i < UE; i++) {
        if (POLYVAL[i] == false) {
            long long int IND = UnElN[i];
            long long int (*COST)[IND] = malloc(IND * sizeof(*COST));
            long int p = 0, q = 0, m = 0, l = 0;

            for (p = 0; p < NATOMS; p++) {
                if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                    for (q = 0; q < NATOMS; q++) {
                        if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                            double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
                            double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
                            double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
                            long long int d = (long long int)(sqrt(x * x + y * y + z * z) * 1e6);
long int ino = ((*coord1tmp)[p].nn > 0 ? (*coord1tmp)[p].n[0] - 1 : -1);
long int inw = -1;
if (ino >= 0 && ino < NATOMS) inw = (*coord1tmp)[ino].N;
if (MOL2 == false) {
if (ino >= 0 && ino < NATOMS && inw >= 0 && inw < NATOMS &&
(*coord1tmp)[q].nn > 0 && (*coord1tmp)[q].n[0] - 1 == inw &&
IDN[ino][inw] == true) {
d = d / 1000;
}
} else {
if (ino >= 0 && ino < NATOMS && inw >= 0 && inw < NATOMS &&
coord2[q].nn > 0 && coord2[q].n[0] - 1 == inw &&
IDN[ino][inw] == true) {
d = d / 1000;
}
}
                            COST[m][l] = d;
                            l = l + 1;
                        }
                    }
                    l = 0;
                    m = m + 1;
                }
            }

            long int *MATCH = malloc(IND * sizeof(long int));
            MC_HUNGMAIN(IND, IND, COST, &MATCH);

            m = 0;
            l = 0;
            for (p = 0; p < NATOMS; p++) {
                if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                    for (q = 0; q < NATOMS; q++) {
                        if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                            if (MATCH[m] == l) {
                                (*coord1tmp)[p].N = q;
                            }
                            l = l + 1;
                        }
                    }
                    l = 0;
                    m = m + 1;
                }
            }
            free(MATCH);
            free(COST);
        }
    }

    // populate coord (temporary array) according to renumbering
    long int *old_to_new = malloc((size_t)NATOMS * sizeof(long int));
    if (!old_to_new) exit(1);
    for (i = 0; i < NATOMS; i++) {
    long int it = (*coord1tmp)[i].N;
    if (it < 0 || it >= NATOMS) {
    fprintf(stderr, "MC_RENUMB: unmatched atom during renumbering\n");
    exit(1);
    }
    old_to_new[i] = it;
    coord[it] = (*coord1tmp)[i];
    coord[it].N = i;
    }
    for (i = 0; i < NATOMS; i++) {
    for (j = 0; j < coord[i].nn && j < 23; j++) {
    long int old_nb = coord[i].n[j] - 1;
    if (old_nb >= 0 && old_nb < NATOMS) {
    coord[i].n[j] = old_to_new[old_nb] + 1;
    }
    }
    }
    free(old_to_new);
    // overwrite coord1tmp
    for (i = 0; i < NATOMS; i++) {
        (*coord1tmp)[i] = coord[i];
        (*coord1tmp)[i].N = coord[i].N;
    }

    free(coord);
}

void MC_RENUMBHUNG(McAtom **coord1tmp, McAtom *coord2, long int NATOMS, char (*UnEl)[4], long int UE) {
    long int i = 0;
    long int j = 0;
    McAtom *coord = malloc(NATOMS * sizeof(McAtom));

    long int UnElN[UE];
    for (i = 0; i < UE; i++) {
        UnElN[i] = 0;
    }
    for (i = 0; i < NATOMS; i++) {
        (*coord1tmp)[i].N = -1;
        for (j = 0; j < UE; j++) {
            if (strcmp(coord2[i].El, UnEl[j]) == 0) {
                UnElN[j] = UnElN[j] + 1;
            }
        }
        coord2[i].N = -1;
    }

    for (i = 0; i < UE; i++) {
        long long int IND = UnElN[i];
        long long int (*COST)[IND] = malloc(IND * sizeof(*COST));
        long int p = 0, q = 0, m = 0, l = 0;

        for (p = 0; p < NATOMS; p++) {
            if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                for (q = 0; q < NATOMS; q++) {
                    if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                        double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
                        double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
                        double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
                        long long int d = (long long int)(sqrt(x * x + y * y + z * z) * 1e6);
                        COST[m][l] = d;
                        l = l + 1;
                    }
                }
                l = 0;
                m = m + 1;
            }
        }

        long int *MATCH = malloc(IND * sizeof(long int));
        MC_HUNGMAIN(IND, IND, COST, &MATCH);

        m = 0;
        l = 0;
        for (p = 0; p < NATOMS; p++) {
            if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                for (q = 0; q < NATOMS; q++) {
                    if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                        if (MATCH[m] == l) {
                            (*coord1tmp)[p].N = q;
                        }
                        l = l + 1;
                    }
                }
                l = 0;
                m = m + 1;
            }
        }
        free(MATCH);
        free(COST);
    }
long int *old_to_new = malloc((size_t)NATOMS * sizeof(long int));
if (!old_to_new) exit(1);
for (i = 0; i < NATOMS; i++) {
long int it = (*coord1tmp)[i].N;
if (it < 0 || it >= NATOMS) {
fprintf(stderr, "MC_RENUMB: unmatched atom during renumbering\n");
exit(1);
}
old_to_new[i] = it;
coord[it] = (*coord1tmp)[i];
coord[it].N = i;
}
for (i = 0; i < NATOMS; i++) {
for (j = 0; j < coord[i].nn && j < 23; j++) {
long int old_nb = coord[i].n[j] - 1;
if (old_nb >= 0 && old_nb < NATOMS) {
coord[i].n[j] = old_to_new[old_nb] + 1;
}
}
}
free(old_to_new);
    for (i = 0; i < NATOMS; i++) {
        (*coord1tmp)[i] = coord[i];
        (*coord1tmp)[i].N = coord[i].N;
    }

    free(coord);
}

void MC_RENUMBHUNGGRAPH(McAtom **coord1tmp, McAtom *coord2, long int NATOMS, char (*UnEl)[4], long int UE, bool IDN[][NATOMS]) {
    long int i = 0;
    long int j = 0;
    McAtom *coord = malloc(NATOMS * sizeof(McAtom));

    long int UnElN[UE];
    for (i = 0; i < UE; i++) {
        UnElN[i] = 0;
    }
    for (i = 0; i < NATOMS; i++) {
        (*coord1tmp)[i].N = -1;
        for (j = 0; j < UE; j++) {
            if (strcmp(coord2[i].El, UnEl[j]) == 0) {
                UnElN[j] = UnElN[j] + 1;
            }
        }
        coord2[i].N = -1;
    }

    for (i = 0; i < UE; i++) {
        long long int IND = UnElN[i];
        long long int (*COST)[IND] = malloc(IND * sizeof(*COST));
        long int p = 0, q = 0, m = 0, l = 0;

        for (p = 0; p < NATOMS; p++) {
            if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
                for (q = 0; q < NATOMS; q++) {
                    if (strcmp(coord2[q].El, UnEl[i]) == 0) {
                        double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
                        double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
                        double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
                        long long int d = (long long int)(sqrt(x * x + y * y + z * z) * 1e6);
                        if (IDN[p][q] == true) {
                            d = d / 1000;
                        }
                        COST[m][l] = d;
                        l = l + 1;
                    }
                }
                l = 0;
                m = m + 1;
            }
        }

long int *MATCH = malloc(IND * sizeof(long int));
MC_HUNGMAIN(IND, IND, COST, &MATCH);
m = 0;
l = 0;
for (p = 0; p < NATOMS; p++) {
if (strcmp((*coord1tmp)[p].El, UnEl[i]) == 0) {
for (q = 0; q < NATOMS; q++) {
if (strcmp(coord2[q].El, UnEl[i]) == 0) {
if (MATCH[m] == l) {
(*coord1tmp)[p].N = q;
}
l = l + 1;
}
}
l = 0;
m = m + 1;
}
}
free(MATCH);
free(COST);
}
long int *old_to_new = malloc((size_t)NATOMS * sizeof(long int));
if (!old_to_new) exit(1);
for (i = 0; i < NATOMS; i++) {
long int it = (*coord1tmp)[i].N;
if (it < 0 || it >= NATOMS) {
fprintf(stderr, "MC_RENUMB: unmatched atom during renumbering\n");
exit(1);
}
old_to_new[i] = it;
coord[it] = (*coord1tmp)[i];
coord[it].N = i;
}
for (i = 0; i < NATOMS; i++) {
for (j = 0; j < coord[i].nn && j < 23; j++) {
long int old_nb = coord[i].n[j] - 1;
if (old_nb >= 0 && old_nb < NATOMS) {
coord[i].n[j] = old_to_new[old_nb] + 1;
}
}
}
free(old_to_new);
for (i = 0; i < NATOMS; i++) {
(*coord1tmp)[i] = coord[i];
(*coord1tmp)[i].N = coord[i].N;
}
free(coord);
}
