#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "mc_types.h"
#include "mc_atominfo.h"
#include "mc_getconnect.h"

static bool mc_isvalueinarray(long int val, long int *arr, long int size, long int start) {
    bool Belong = false;
    long int i;
    for (i = start; i < size; i++) {
        if (arr[i] == val) {
            Belong = true;
        }
    }
    return Belong;
}

static double MC_CONN_ANGLE(double xyz1[3], double xyz2[3], double xyz3[3]) {
    double angle = 0;

    // get the vectors
    double v1[3] = {xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2]};
    double v2[3] = {xyz3[0] - xyz2[0], xyz3[1] - xyz2[1], xyz3[2] - xyz2[2]};
double v1mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
double v2mag = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
if (v1mag < 1e-12 || v2mag < 1e-12) return 0.0;
double v1norm[3] = {v1[0]/v1mag, v1[1]/v1mag, v1[2]/v1mag};
double v2norm[3] = {v2[0]/v2mag, v2[1]/v2mag, v2[2]/v2mag};

    // calculate the dot product
    double res = v1norm[0]*v2norm[0] + v1norm[1]*v2norm[1] + v1norm[2]*v2norm[2];
    if (res > 1.0) res = 1.0;
    if (res < -1.0) res = -1.0;

    angle = acos(res) * (180.0 / M_PI);
    return angle;
}

void MC_GETCONNECT(McAtom **coordinates, long int NATOMS) {
    long int i = 0, j = 0, i1 = 0;
    long int maxbond = 0;
    long int *bonds;
    double *dist;
    double d;

    for (i = 0; i < NATOMS; i++) {
        maxbond = 0;
        bonds = malloc(1 * sizeof(long int));
        dist = malloc(1 * sizeof(double));

        for (j = 0; j < NATOMS; j++) {
            // skip comparison with itself
            if (i != j) {
                d = sqrt(pow(((*coordinates)[i].xyz[0] - (*coordinates)[j].xyz[0]), 2)
                       + pow(((*coordinates)[i].xyz[1] - (*coordinates)[j].xyz[1]), 2)
                       + pow(((*coordinates)[i].xyz[2] - (*coordinates)[j].xyz[2]), 2));
                // check whether the distance is less than sum of covalent radii + 0.45A
                if (d < MC_ATOMINFO((*coordinates)[i].El).R + MC_ATOMINFO((*coordinates)[j].El).R + 0.45) {
                    bonds[maxbond] = j + 1;
                    dist[maxbond] = d;
                    maxbond++;
                    bonds = (long int *)realloc(bonds, (maxbond + 1) * sizeof(long int));
                    dist = (double *)realloc(dist, (maxbond + 1) * sizeof(double));
                }
            }
        }

        // check with maximum valency
        if (maxbond <= MC_ATOMINFO((*coordinates)[i].El).V) {
            (*coordinates)[i].nn = maxbond;
            for (j = 0; j < maxbond; j++) {
                (*coordinates)[i].n[j] = bonds[j];
            }
        } else {
            // bubble sorting to remove the largest distances
            for (j = 0; j < maxbond; j++) {
                for (i1 = 0; i1 < (maxbond - j - 1); i1++) {
                    if (dist[i1] > dist[i1 + 1]) {
                        long int N;
                        double dist1;
                        dist1 = dist[i1];
                        N = bonds[i1];
                        dist[i1] = dist[i1 + 1];
                        bonds[i1] = bonds[i1 + 1];
                        dist[i1 + 1] = dist1;
                        bonds[i1 + 1] = N;
                    }
                }
            }

            (*coordinates)[i].nn = MC_ATOMINFO((*coordinates)[i].El).V;
            for (j = 0; j < (*coordinates)[i].nn; j++) {
                (*coordinates)[i].n[j] = bonds[j];
            }
        }

        free(bonds);
        free(dist);
    }

    // get rid of some bonds which give the angle less than 50
    for (i = 0; i < NATOMS; i++) {
        long int pairs[1000][2];
        double angle;
        long int count = 0;
        for (i1 = 0; i1 < (*coordinates)[i].nn; i1++) {
            for (j = i1 + 1; j < (*coordinates)[i].nn; j++) {
                angle = MC_CONN_ANGLE((*coordinates)[(*coordinates)[i].n[i1] - 1].xyz,
                                      (*coordinates)[i].xyz,
                                      (*coordinates)[(*coordinates)[i].n[j] - 1].xyz);
                if (angle < 50.0) {
                    pairs[count][0] = (*coordinates)[i].n[i1];
                    pairs[count][1] = (*coordinates)[i].n[j];
                    long int z = pairs[count][0];
                    long int z1 = pairs[count][1];
                    d = sqrt(pow(((*coordinates)[i].xyz[0] - (*coordinates)[z - 1].xyz[0]), 2)
                           + pow(((*coordinates)[i].xyz[1] - (*coordinates)[z - 1].xyz[1]), 2)
                           + pow(((*coordinates)[i].xyz[2] - (*coordinates)[z - 1].xyz[2]), 2));
                    double cr = MC_ATOMINFO((*coordinates)[i].El).R + MC_ATOMINFO((*coordinates)[z - 1].El).R + 0.45;

                    double d1 = sqrt(pow(((*coordinates)[i].xyz[0] - (*coordinates)[z1 - 1].xyz[0]), 2)
                                   + pow(((*coordinates)[i].xyz[1] - (*coordinates)[z1 - 1].xyz[1]), 2)
                                   + pow(((*coordinates)[i].xyz[2] - (*coordinates)[z1 - 1].xyz[2]), 2));
                    double cr1 = MC_ATOMINFO((*coordinates)[i].El).R + MC_ATOMINFO((*coordinates)[z1 - 1].El).R + 0.45;

                    if ((cr - d) < (cr1 - d1)) {
                        pairs[count][0] = z;
                    } else {
                        pairs[count][0] = z1;
                    }
                    count++;
                }
            }
        }

        if ((*coordinates)[i].nn <= 1) continue;
            long int atoms[ 23 ];
        long int count1 = 0;
        for (i1 = 0; i1 < (*coordinates)[i].nn; i1++) {
            bool IN = false;
            for (j = 0; j < count; j++) {
                if ((*coordinates)[i].n[i1] == pairs[j][0]) {
                    IN = true;
                    break;
                }
            }
            if (IN == false) {
                atoms[count1] = (*coordinates)[i].n[i1];
                count1++;
            }
        }

        (*coordinates)[i].nn = count1;
        for (i1 = 0; i1 < (*coordinates)[i].nn; i1++) {
            (*coordinates)[i].n[i1] = atoms[i1];
        }

        // disconnect the atoms which were connected to these elements
        for (i1 = 0; i1 < count; i1++) {
            long int z = pairs[i1][0];
            if (mc_isvalueinarray(i + 1, (*coordinates)[z - 1].n, (*coordinates)[z - 1].nn, 0) == true) {
                long int atoms1[(*coordinates)[z - 1].nn - 1];
                long int count_rem = 0;
                for (j = 0; j < (*coordinates)[z - 1].nn; j++) {
                    if ((*coordinates)[z - 1].n[j] != i + 1) {
                        atoms1[count_rem] = (*coordinates)[z - 1].n[j];
                        count_rem++;
                    }
                }
                (*coordinates)[z - 1].nn = (*coordinates)[z - 1].nn - 1;
                for (j = 0; j < (*coordinates)[z - 1].nn; j++) {
                    (*coordinates)[z - 1].n[j] = atoms1[j];
                }
            }
        }
    }

    MC_SYMMETRYCHECK(coordinates, NATOMS);
}

void MC_SYMMETRYCHECK(McAtom **coordinates, long int NATOMS) {
    // check symmetry: if 1 is connected to 2, then 2 should be connected to 1
    long int i = 0, j = 0, i1 = 0;
    for (i = 0; i < NATOMS; i++) {
        long int N = 0;
        long int DEL[23];
        for (j = 0; j < (*coordinates)[i].nn; j++) {
            long int ind = (*coordinates)[i].n[j] - 1;
            bool F = false;
            for (i1 = 0; i1 < (*coordinates)[ind].nn; i1++) {
                if ((i + 1) == (*coordinates)[ind].n[i1]) {
                    F = true;
                    break;
                }
            }
            if (F == false) {
                DEL[N] = ind + 1;
                N = N + 1;
            }
        }
        if (N > 0) {
            long int nn = (*coordinates)[i].nn - N;
            if (nn < 0) nn = 0;
            if (nn > 23) nn = 23;
            long int tmp[23];
            long int x = 0;
            for (j = 0; j < (*coordinates)[i].nn; j++) {
                bool Found = mc_isvalueinarray((*coordinates)[i].n[j], DEL, N, 0);
                if (Found == false) {
                    tmp[x] = (*coordinates)[i].n[j];
                    x = x + 1;
                }
            }
            (*coordinates)[i].nn = (*coordinates)[i].nn - N;
            for (j = 0; j < x; j++) {
                (*coordinates)[i].n[j] = tmp[j];
            }
        }
    }
}

void MC_SYMMETRYCHECKDEF(McAtom *coorddef, McAtom **coord1, long int NATOMS) {
    long int i = 0;
    long int j = 0;
    long int k = 0;
    for (i = 0; i < NATOMS; i++) {
        for (j = 0; j < (*coord1)[i].nn; j++) {
            bool F = mc_isvalueinarray((*coord1)[i].n[j], coorddef[i].n, coorddef[i].nn, 0);
            if (F == false) {
                long int m = (*coord1)[i].n[j] - 1;
                bool F1 = mc_isvalueinarray(i + 1, (*coord1)[m].n, (*coord1)[m].nn, 0);
                if (F1 == false) {
                    (*coord1)[m].nn = (*coord1)[m].nn + 1;
                    long int n = (*coord1)[m].nn;
                    (*coord1)[m].n[n - 1] = i + 1;
                }
            }
        }
        for (j = 0; j < (coorddef)[i].nn; j++) {
            bool F = mc_isvalueinarray((coorddef)[i].n[j], (*coord1)[i].n, (*coord1)[i].nn, 0);
            if (F == false) {
                long int m = (coorddef)[i].n[j] - 1;
                bool F1 = mc_isvalueinarray(i + 1, (*coord1)[m].n, (*coord1)[m].nn, 0);
                if (F1 == true) {
                    (*coord1)[m].nn = (*coord1)[m].nn - 1;
                    long int n = (*coord1)[m].nn;
                    long int new_arr[n];
                    long int count = 0;
                    for (k = 0; k < (*coord1)[m].nn + 1; k++) {
                        if ((*coord1)[m].n[k] != i + 1) {
                            new_arr[count] = (*coord1)[m].n[k];
                            count = count + 1;
                        }
                    }

for (k = 0; k < (*coord1)[m].nn; k++) {
(*coord1)[m].n[k] = new_arr[k];
}
}
}
}
}
}
