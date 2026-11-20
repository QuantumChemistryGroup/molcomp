#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mc_types.h"
#include "mc_stdfunc.h"

static inline int comp(const void *a, const void *b) {
    long int va = *(const long int *)a;
    long int vb = *(const long int *)b;
    return (va > vb) - (va < vb);
}

bool MC_BONDEQUAL(long int p, long int q, McAtom *coord1, McAtom *coord2, long int NATOMS, bool MOL2) {
    (void)NATOMS;
    bool EQ = true;

    if (MOL2 == false) {
        if (coord1[p].nn != coord1[q].nn) {
            return false;
        }
    } else {
        if (coord1[p].nn != coord2[q].nn) {
            return false;
        }
    }
long int N = coord1[p].nn;
if (N <= 0) return EQ;
if (N > 23) return false;
long int array1[23];
long int array2[23];
    for (long int i = 0; i < N; i++) {
        long int seq1 = coord1[p].n[i] - 1;
        array1[i] = (long int)coord1[seq1].Z;
        if (MOL2 == false) {
            long int seq2 = coord1[q].n[i] - 1;
            array2[i] = (long int)coord1[seq2].Z;
        } else {
            long int seq2 = coord2[q].n[i] - 1;
            array2[i] = (long int)coord2[seq2].Z;
        }
    }

    qsort(array1, N, sizeof(long int), comp);
    qsort(array2, N, sizeof(long int), comp);

    for (long int i = 0; i < N; i++) {
        if (array1[i] != array2[i]) {
            EQ = false;
            break;
        }
    }

    return EQ;
}

void MC_POPADJ(long int NATOMS, long int ADJ[][NATOMS], McAtom *coord) {
    for (long int j = 0; j < NATOMS; j++) {
        for (long int l = 0; l < NATOMS; l++) {
            ADJ[j][l] = 0;
        }
    }

    for (long int j = 0; j < NATOMS; j++) {
        for (long int l = 0; l < coord[j].nn; l++) {
            long int target_row = coord[j].n[l] - 1;
            ADJ[target_row][j] = 1;
        }
    }
}

bool MC_RADIAL(long int p, long int q, McAtom *coord1, McAtom *coord2, long int NATOMS, bool MOL2, long int ADJ1[][NATOMS], long int ADJ2[][NATOMS]) {
    (void)ADJ1;
    (void)ADJ2;
    bool EQ = true;

    /* FIRST QUICK CHECK - immediate surrounding coord number */
    if (MOL2 == false) {
        if (coord1[p].nn != coord1[q].nn) {
            return false;
        }
    } else {
        if (coord1[p].nn != coord2[q].nn) {
            return false;
        }
    }

    /* SECOND QUICK CHECK - immediate surrounding elements */
    long int N = coord1[p].nn;
    if (N <= 0) return EQ;
    if (N > 23) return false;
    long int array1[23];
    long int array2[23];
    for (long int i = 0; i < N; i++) {
        long int seq1 = coord1[p].n[i] - 1;
        array1[i] = (long int)coord1[seq1].Z;
        if (MOL2 == false) {
            long int seq2 = coord1[q].n[i] - 1;
            array2[i] = (long int)coord1[seq2].Z;
        } else {
            long int seq2 = coord2[q].n[i] - 1;
            array2[i] = (long int)coord2[seq2].Z;
        }
    }

    qsort(array1, N, sizeof(long int), comp);
    qsort(array2, N, sizeof(long int), comp);

    for (long int i = 0; i < N; i++) {
        if (array1[i] != array2[i]) {
            EQ = false;
            break;
        }
    }

    /* BFS search */
    if (EQ != false) {
        long int visited1[NATOMS];
        long int visited2[NATOMS];

        for (long int ind = 0; ind < NATOMS; ind++) {
            visited1[ind] = 0;
            visited2[ind] = 0;
        }

        long int queue1[NATOMS];
        long int front1 = 0;
        long int rear1 = 0;

        long int queue2[NATOMS];
        long int front2 = 0;
        long int rear2 = 0;

        visited1[p] = 1;
        queue1[rear1++] = p;

        visited2[q] = 1;
        queue2[rear2++] = q;

        long int level1[NATOMS];
        long int level2[NATOMS];
        for (long int i = 0; i < NATOMS; i++) {
            level1[i] = 0;
            level2[i] = 0;
        }

        long int pl1 = 0;
        long int pl2 = 0;
        long int pl1arr[NATOMS];
        long int pl2arr[NATOMS];
        long int plind = 0;
        for (long int i = 0; i < NATOMS; i++) {
            pl1arr[i] = 0;
            pl2arr[i] = 0;
        }

        while ((front1 < rear1) || (front2 < rear2)) {
            if (!((front1 < rear1) && (front2 < rear2))) {
                return false;
            }
            long int u1 = queue1[front1++];
            long int u2 = queue2[front2++];

            if ((pl1 != pl2) || (level1[u1] != level2[u2])) {
                return false;
            }

            if (pl1 != level1[u1]) {
                for (long int j = 0; j < plind; j++) {
                    long int k = pl1arr[j];
                    long int m = pl2arr[j];
                    pl1arr[j] = (long int)coord1[k].Z;
                    if (MOL2 == true) {
                        pl2arr[j] = (long int)coord2[m].Z;
                    } else {
                        pl2arr[j] = (long int)coord1[m].Z;
                    }
                }
                qsort(pl1arr, plind, sizeof(long int), comp);
                qsort(pl2arr, plind, sizeof(long int), comp);

                for (long int j = 0; j < plind; j++) {
                    if (pl1arr[j] != pl2arr[j]) {
                        return false;
                    }
                }

                for (long int j = 0; j < plind; j++) {
                    pl1arr[j] = 0;
                    pl2arr[j] = 0;
                }
                plind = 0;
                pl1arr[plind] = u1;
                pl2arr[plind] = u2;
                plind++;
            } else {
                pl1arr[plind] = u1;
                pl2arr[plind] = u2;
                plind++;
            }

            // neighbors of u1
            for (long int i = 0; i < coord1[u1].nn; i++) {
                long int idx = coord1[u1].n[i] - 1;
                if (visited1[idx] == 0) {
                    visited1[idx] = 1;
                    level1[idx] = level1[u1] + 1;
                    queue1[rear1++] = idx;
                }
            }
            pl1 = level1[u1];

            // neighbors of u2
            if (MOL2 == false) {
                for (long int q_i = 0; q_i < coord1[u2].nn; q_i++) {
                    long int idx = coord1[u2].n[q_i] - 1;
                    if (visited2[idx] == 0) {
                        visited2[idx] = 1;
                        level2[idx] = level2[u2] + 1;
                        queue2[rear2++] = idx;
                    }
                }
            } else {
                for (long int q_i = 0; q_i < coord2[u2].nn; q_i++) {
                    long int idx = coord2[u2].n[q_i] - 1;
                    if (visited2[idx] == 0) {
                        visited2[idx] = 1;
                        level2[idx] = level2[u2] + 1;
                        queue2[rear2++] = idx;
                    }
                }
            }
            pl2 = level2[u2];
            }
            if (plind > 0) {
            for (long int j = 0; j < plind; j++) {
            long int k = pl1arr[j];
            long int m = pl2arr[j];
            pl1arr[j] = (long int)coord1[k].Z;
            if (MOL2 == true) {
            pl2arr[j] = (long int)coord2[m].Z;
            } else {
            pl2arr[j] = (long int)coord1[m].Z;
            }
            }
            qsort(pl1arr, plind, sizeof(long int), comp);
            qsort(pl2arr, plind, sizeof(long int), comp);
            for (long int j = 0; j < plind; j++) {
            if (pl1arr[j] != pl2arr[j]) {
            return false;
            }
            }
            }
            }
            return EQ;
            }

