#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "mc_types.h"
#include "mc_kabshangle.h"

static double MC_DIHEDRAL(double p1[3], double p2[3], double p3[3], double p4[3]);
static double MC_ANGLEB(double xyz1[3], double xyz2[3], double xyz3[3]);
static void MC_GETVEC(double atom2[3], double atom1[3], double v[3]);
static void MC_GETUNITVEC(double v[3], double u[3]);
static double mc_scalar_product(double v1[3], double v2[3]);
static void mc_vector_product(double u[3], double w[3], double v[3]);

// Jacobi eigensolver for real symmetric 4x4 matrix (Horn's quaternion Kabsch method)
static void jacobi_sym4(double A[4][4], double eval[4], double evec[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            evec[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < 60; sweep++) {
        double off = 0.0;
        for (int p = 0; p < 3; p++) {
            for (int q = p + 1; q < 4; q++) {
                off += fabs(A[p][q]);
            }
        }
        if (off < 1e-15) break;

        for (int p = 0; p < 3; p++) {
            for (int q = p + 1; q < 4; q++) {
                double apq = A[p][q];
                if (fabs(apq) < 1e-16) continue;

                double app = A[p][p];
                double aqq = A[q][q];
                double tau = (aqq - app) / (2.0 * apq);
                double t = (tau >= 0.0)
                    ? (1.0 / (tau + sqrt(1.0 + tau * tau)))
                    : (-1.0 / (-tau + sqrt(1.0 + tau * tau)));
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;

                A[p][p] = app - t * apq;
                A[q][q] = aqq + t * apq;
                A[p][q] = A[q][p] = 0.0;

                for (int r = 0; r < 4; r++) {
                    if (r == p || r == q) continue;
                    double arp = A[r][p];
                    double arq = A[r][q];
                    A[r][p] = A[p][r] = c * arp - s * arq;
                    A[r][q] = A[q][r] = s * arp + c * arq;
                }

                for (int k = 0; k < 4; k++) {
                    double vkp = evec[k][p];
                    double vkq = evec[k][q];
                    evec[k][p] = c * vkp - s * vkq;
                    evec[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        eval[i] = A[i][i];
    }
}

void MC_KABSHANGLE(
    McAtom *coord1,
    McAtom *coord2,
    McAtom **coord1mod,
    long int NATOMS,
    double *RMSD,
    double *AAV,
    double *AMAX,
    bool MOL2,
    double LINEARTHLD,
    bool VER
) {
    long int i = 0;
    long int j = 0;

    // keep unperturbed structures metadata
    for (i = 0; i < NATOMS; i++) {
        strcpy((*coord1mod)[i].El, coord1[i].El);
        for (j = 0; j < 3; j++) {
            (*coord1mod)[i].xyz[j] = coord1[i].xyz[j];
        }
        (*coord1mod)[i].N = coord1[i].N;
    }

    // calculate centroids
    double cent1[3] = {0, 0, 0};
    double cent2[3] = {0, 0, 0};

    for (i = 0; i < NATOMS; i++) {
        for (j = 0; j < 3; j++) {
            cent1[j] = coord1[i].xyz[j] + cent1[j];
            cent2[j] = coord2[i].xyz[j] + cent2[j];
        }
    }

    for (i = 0; i < 3; i++) {
        cent1[i] = cent1[i] / NATOMS;
        cent2[i] = cent2[i] / NATOMS;
    }

    // calculate covariance matrix H = (P^T)Q (3x3)
    double H[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (i = 0; i < NATOMS; i++) {
        double p0 = coord1[i].xyz[0] - cent1[0];
        double p1 = coord1[i].xyz[1] - cent1[1];
        double p2 = coord1[i].xyz[2] - cent1[2];

        double q0 = coord2[i].xyz[0] - cent2[0];
        double q1 = coord2[i].xyz[1] - cent2[1];
        double q2 = coord2[i].xyz[2] - cent2[2];

        H[0][0] += p0 * q0; H[0][1] += p0 * q1; H[0][2] += p0 * q2;
        H[1][0] += p1 * q0; H[1][1] += p1 * q1; H[1][2] += p1 * q2;
        H[2][0] += p2 * q0; H[2][1] += p2 * q1; H[2][2] += p2 * q2;
    }

    // Horn's 4x4 symmetric matrix N for optimal proper rotation (det(R) = +1)
    double Sxx = H[0][0], Sxy = H[0][1], Sxz = H[0][2];
    double Syx = H[1][0], Syy = H[1][1], Syz = H[1][2];
    double Szx = H[2][0], Szy = H[2][1], Szz = H[2][2];

    double Nmat[4][4] = {
        { Sxx + Syy + Szz,  Syz - Szy,        Szx - Sxz,        Sxy - Syx       },
        { Syz - Szy,        Sxx - Syy - Szz,  Sxy + Syx,        Szx + Sxz       },
        { Szx - Sxz,        Sxy + Syx,       -Sxx + Syy - Szz,  Syz + Szy       },
        { Sxy - Syx,        Szx + Sxz,        Syz + Szy,       -Sxx - Syy + Szz }
    };

    double eval[4];
    double evec[4][4];
    jacobi_sym4(Nmat, eval, evec);

    int max_idx = 0;
    for (int k = 1; k < 4; k++) {
        if (eval[k] > eval[max_idx]) max_idx = k;
    }

    double q0 = evec[0][max_idx];
    double q1 = evec[1][max_idx];
    double q2 = evec[2][max_idx];
    double q3 = evec[3][max_idx];
    double qnorm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    if (qnorm > 1e-15) {
        q0 /= qnorm; q1 /= qnorm; q2 /= qnorm; q3 /= qnorm;
    } else {
        q0 = 1.0; q1 = q2 = q3 = 0.0;
    }

    // Rotation matrix R (row-vector convention: PROT = P * R)
    double R[3][3];
    R[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R[0][1] = 2.0 * (q1*q2 + q0*q3);
    R[0][2] = 2.0 * (q1*q3 - q0*q2);

    R[1][0] = 2.0 * (q1*q2 - q0*q3);
    R[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    R[1][2] = 2.0 * (q2*q3 + q0*q1);

    R[2][0] = 2.0 * (q1*q3 + q0*q2);
    R[2][1] = 2.0 * (q2*q3 - q0*q1);
    R[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    // get all dihedrals
    McMM *tors = malloc(1 * sizeof(*tors));
    McAtom *coord = (MOL2 == true) ? coord2 : coord1;

    long int l = 0;
    long int q = 0;
    long int m = 0;
    for (i = 0; i < NATOMS; i++) {
        for (j = 0; j < coord[i].nn; j++) {
            long int z = coord[i].n[j];
            for (m = 0; m < coord[i].nn; m++) {
                for (q = 0; q < coord[z - 1].nn; q++) {
                    if (coord[i].n[m] == z || coord[z - 1].n[q] == i + 1) {
                    } else {
                        tors[l].A[0] = coord[i].n[m];
                        tors[l].A[1] = i + 1;
                        tors[l].A[2] = z;
                        tors[l].A[3] = coord[z - 1].n[q];
                        l++;
                        tors = realloc(tors, (l + 1) * sizeof(*tors));
                    }
                }
            }
        }
    }

    (*RMSD) = 0;
    double Aav = 0;  // average angle
    double Amax = 0; // maximum angle difference
    double Dmax = 0; // maximum distance between aligned atoms

    for (i = 0; i < NATOMS; i++) {
        double p0 = coord1[i].xyz[0] - cent1[0];
        double p1 = coord1[i].xyz[1] - cent1[1];
        double p2 = coord1[i].xyz[2] - cent1[2];

        double prot[3];
        prot[0] = p0 * R[0][0] + p1 * R[1][0] + p2 * R[2][0];
        prot[1] = p0 * R[0][1] + p1 * R[1][1] + p2 * R[2][1];
        prot[2] = p0 * R[0][2] + p1 * R[1][2] + p2 * R[2][2];

        double d = 0; // distance between atom pair
        for (j = 0; j < 3; j++) {
            (*coord1mod)[i].xyz[j] = prot[j];
            double x0 = prot[j];
            double x1 = coord2[i].xyz[j] - cent2[j];
            double dx2 = (x0 - x1) * (x0 - x1);
            (*RMSD) = dx2 + (*RMSD);
            d = dx2 + d;
        }
        d = sqrt(d);
        if (d > Dmax) {
            Dmax = d;
        }
    }

    long int count = 0;
    // compare angles
    for (i = 0; i < l; i++) {
        double A1 = MC_DIHEDRAL(coord1[tors[i].A[0] - 1].xyz, coord1[tors[i].A[1] - 1].xyz,
                                coord1[tors[i].A[2] - 1].xyz, coord1[tors[i].A[3] - 1].xyz);
        double A2 = MC_DIHEDRAL(coord2[tors[i].A[0] - 1].xyz, coord2[tors[i].A[1] - 1].xyz,
                                coord2[tors[i].A[2] - 1].xyz, coord2[tors[i].A[3] - 1].xyz);

        double A11 = MC_ANGLEB(coord1[tors[i].A[0] - 1].xyz, coord1[tors[i].A[1] - 1].xyz, coord1[tors[i].A[2] - 1].xyz);
        double A12 = MC_ANGLEB(coord1[tors[i].A[1] - 1].xyz, coord1[tors[i].A[2] - 1].xyz, coord1[tors[i].A[3] - 1].xyz);
        double A21 = MC_ANGLEB(coord2[tors[i].A[0] - 1].xyz, coord2[tors[i].A[1] - 1].xyz, coord2[tors[i].A[2] - 1].xyz);
        double A22 = MC_ANGLEB(coord2[tors[i].A[1] - 1].xyz, coord2[tors[i].A[2] - 1].xyz, coord2[tors[i].A[3] - 1].xyz);

        // threshold for linearity
        double thld = LINEARTHLD;
        if ((fabs(180 - A11) < thld) || (fabs(180 - A12) < thld) ||
            (fabs(180 - A21) < thld) || (fabs(180 - A22) < thld)) {
            continue;
        }

        double DIFF = (A2 - A1);
        if (DIFF > 180.0) {
            DIFF = DIFF - 360;
        } else if (DIFF <= -180.0) {
            DIFF = DIFF + 360;
        }
        DIFF = fabs(DIFF);
        if ((isnan(A1) == true) || (isnan(A2) == true)) {
            continue;
        }
        Aav = Aav + DIFF;
        if (DIFF > Amax) {
            Amax = DIFF;
            if (VER) {
                printf("i:%li   j:%li   k:%li   l:%li   A1:%f  A2:%f DIFF: %f count=%li\n",
                       tors[i].A[0], tors[i].A[1], tors[i].A[2], tors[i].A[3], A1, A2, DIFF, count);
            }
        }
        count = count + 1;
    }

    if (count != 0) {
        Aav = Aav / count;
    } else {
        Aav = 0;
    }

    (*RMSD) = sqrt((1.0 / NATOMS) * (*RMSD));
    if (VER) printf("RMSD = %f\n", (*RMSD));
    (*AAV) = Dmax; // middle value is the largest distance
    (*AMAX) = Amax;

    free(tors);
}

static double MC_DIHEDRAL(double p1[3], double p2[3], double p3[3], double p4[3]) {
    double b1[3], b2[3], b3[3];
    MC_GETVEC(p2, p1, b1);
    MC_GETVEC(p3, p2, b2);
    MC_GETVEC(p4, p3, b3);

    double b1xb2[3];
    double b2xb3[3];
    mc_vector_product(b1, b2, b1xb2);
    mc_vector_product(b2, b3, b2xb3);

    double n1[3], n2[3], b2x[3], m1[3];
    MC_GETUNITVEC(b1xb2, n1);

MC_GETUNITVEC(b2xb3, n2);
MC_GETUNITVEC(b2, b2x);
mc_vector_product(n1, b2x, m1);
double x = mc_scalar_product(n1, n2);
double y = mc_scalar_product(m1, n2);
return atan2(y, x) * 180.0 / M_PI;
}
static void MC_GETVEC(double atom2[3], double atom1[3], double v[3]) {
for (int i = 0; i < 3; i++) {
v[i] = atom2[i] - atom1[i];
}
}
static void MC_GETUNITVEC(double v[3], double u[3]) {
double lenv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
if (lenv < 1e-12) {
u[0] = 0.0; u[1] = 0.0; u[2] = 0.0;
return;
}
for (int i = 0; i < 3; i++) {
u[i] = v[i] / lenv;
}
}
static double mc_scalar_product(double v1[3], double v2[3]) {
return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}
static void mc_vector_product(double u[3], double w[3], double v[3]) {
v[0] = u[1]*w[2] - u[2]*w[1];
v[1] = u[2]*w[0] - u[0]*w[2];
v[2] = u[0]*w[1] - u[1]*w[0];
}
static double MC_ANGLEB(double xyz1[3], double xyz2[3], double xyz3[3]) {
double v1[3] = {xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2]};
double v2[3] = {xyz3[0] - xyz2[0], xyz3[1] - xyz2[1], xyz3[2] - xyz2[2]};
double v1mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
double v2mag = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
if (v1mag < 1e-12 || v2mag < 1e-12) return 0.0;
double res = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) / (v1mag * v2mag);
if (res > 1.0) res = 1.0;
if (res < -1.0) res = -1.0;
return acos(res) * (180.0 / M_PI);
}

void MC_KABSHANGLE_FAST(
    const double *xyz1,
    const double *xyz2,
    double *xyz1_aligned,
    long int NATOMS,
    const McMM *tors,
    long int ntors,
    bool compute_angles,
    double *RMSD,
    double *DMAX,
    double *AMAX,
    double LINEARTHLD,
    bool VER
) {
    double cent1[ 3 ] = {0.0, 0.0, 0.0};
    double cent2[ 3 ] = {0.0, 0.0, 0.0};

    for (long int i = 0; i < NATOMS; i++) {
        cent1[ 0 ] += *(xyz1 + 3 * i + 0);
        cent1[ 1 ] += *(xyz1 + 3 * i + 1);
        cent1[ 2 ] += *(xyz1 + 3 * i + 2);
        cent2[ 0 ] += *(xyz2 + 3 * i + 0);
        cent2[ 1 ] += *(xyz2 + 3 * i + 1);
        cent2[ 2 ] += *(xyz2 + 3 * i + 2);
    }
    double inv_n = 1.0 / (double)NATOMS;
    cent1[ 0 ] *= inv_n; cent1[ 1 ] *= inv_n; cent1[ 2 ] *= inv_n;
    cent2[ 0 ] *= inv_n; cent2[ 1 ] *= inv_n; cent2[ 2 ] *= inv_n;

    double H[ 3 ][ 3 ] = {{0.0}};
    for (long int i = 0; i < NATOMS; i++) {
        double p0 = *(xyz1 + 3 * i + 0) - cent1[ 0 ];
        double p1 = *(xyz1 + 3 * i + 1) - cent1[ 1 ];
        double p2 = *(xyz1 + 3 * i + 2) - cent1[ 2 ];

        double q0 = *(xyz2 + 3 * i + 0) - cent2[ 0 ];
        double q1 = *(xyz2 + 3 * i + 1) - cent2[ 1 ];
        double q2 = *(xyz2 + 3 * i + 2) - cent2[ 2 ];

        H[ 0 ][ 0 ] += p0 * q0; H[ 0 ][ 1 ] += p0 * q1; H[ 0 ][ 2 ] += p0 * q2;
        H[ 1 ][ 0 ] += p1 * q0; H[ 1 ][ 1 ] += p1 * q1; H[ 1 ][ 2 ] += p1 * q2;
        H[ 2 ][ 0 ] += p2 * q0; H[ 2 ][ 1 ] += p2 * q1; H[ 2 ][ 2 ] += p2 * q2;
    }

    double Sxx = H[ 0 ][ 0 ], Sxy = H[ 0 ][ 1 ], Sxz = H[ 0 ][ 2 ];
    double Syx = H[ 1 ][ 0 ], Syy = H[ 1 ][ 1 ], Syz = H[ 1 ][ 2 ];
    double Szx = H[ 2 ][ 0 ], Szy = H[ 2 ][ 1 ], Szz = H[ 2 ][ 2 ];

    double Nmat[ 4 ][ 4 ] = {
        { Sxx + Syy + Szz,  Syz - Szy,        Szx - Sxz,        Sxy - Syx       },
        { Syz - Szy,        Sxx - Syy - Szz,  Sxy + Syx,        Szx + Sxz       },
        { Szx - Sxz,        Sxy + Syx,       -Sxx + Syy - Szz,  Syz + Szy       },
        { Sxy - Syx,        Szx + Sxz,        Syz + Szy,       -Sxx - Syy + Szz }
    };

    double eval[ 4 ];
    double evec[ 4 ][ 4 ];
    jacobi_sym4(Nmat, eval, evec);

    int max_idx = 0;
    for (int k = 1; k < 4; k++) {
        if (eval[ k ] > eval[ max_idx ]) max_idx = k;
    }

    double q0 = evec[ 0 ][ max_idx ];
    double q1 = evec[ 1 ][ max_idx ];
    double q2 = evec[ 2 ][ max_idx ];
    double q3 = evec[ 3 ][ max_idx ];
    double qnorm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    if (qnorm > 1e-15) {
        q0 /= qnorm; q1 /= qnorm; q2 /= qnorm; q3 /= qnorm;
    } else {
        q0 = 1.0; q1 = q2 = q3 = 0.0;
    }

    double R[ 3 ][ 3 ];
    R[ 0 ][ 0 ] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R[ 0 ][ 1 ] = 2.0 * (q1*q2 + q0*q3);
    R[ 0 ][ 2 ] = 2.0 * (q1*q3 - q0*q2);

    R[ 1 ][ 0 ] = 2.0 * (q1*q2 - q0*q3);
    R[ 1 ][ 1 ] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    R[ 1 ][ 2 ] = 2.0 * (q2*q3 + q0*q1);

    R[ 2 ][ 0 ] = 2.0 * (q1*q3 + q0*q2);
    R[ 2 ][ 1 ] = 2.0 * (q2*q3 - q0*q1);
    R[ 2 ][ 2 ] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    double rmsd_sum = 0.0;
    double dmax = 0.0;

    for (long int i = 0; i < NATOMS; i++) {
        double p0 = *(xyz1 + 3 * i + 0) - cent1[ 0 ];
        double p1 = *(xyz1 + 3 * i + 1) - cent1[ 1 ];
        double p2 = *(xyz1 + 3 * i + 2) - cent1[ 2 ];

        double rx = p0 * R[ 0 ][ 0 ] + p1 * R[ 1 ][ 0 ] + p2 * R[ 2 ][ 0 ];
        double ry = p0 * R[ 0 ][ 1 ] + p1 * R[ 1 ][ 1 ] + p2 * R[ 2 ][ 1 ];
        double rz = p0 * R[ 0 ][ 2 ] + p1 * R[ 1 ][ 2 ] + p2 * R[ 2 ][ 2 ];

        if (xyz1_aligned) {
            *(xyz1_aligned + 3 * i + 0) = rx;
            *(xyz1_aligned + 3 * i + 1) = ry;
            *(xyz1_aligned + 3 * i + 2) = rz;
        }

        double dx = rx - (*(xyz2 + 3 * i + 0) - cent2[ 0 ]);
        double dy = ry - (*(xyz2 + 3 * i + 1) - cent2[ 1 ]);
        double dz = rz - (*(xyz2 + 3 * i + 2) - cent2[ 2 ]);
        double d2 = dx * dx + dy * dy + dz * dz;
        rmsd_sum += d2;
        if (compute_angles) {
            double d = sqrt(d2);
            if (d > dmax) dmax = d;
        }
    }

    *RMSD = sqrt(inv_n * rmsd_sum);
    if (DMAX) *DMAX = dmax;

    if (!compute_angles) {
        if (AMAX) *AMAX = 0.0;
        return;
    }

    double amax = 0.0;
    long int count = 0;
    for (long int i = 0; i < ntors; i++) {
        long int a0 = tors[ i ].A[ 0 ] - 1;
        long int a1 = tors[ i ].A[ 1 ] - 1;
        long int a2 = tors[ i ].A[ 2 ] - 1;
        long int a3 = tors[ i ].A[ 3 ] - 1;

        double p1_0[ 3 ] = { *(xyz1 + 3*a0 + 0), *(xyz1 + 3*a0 + 1), *(xyz1 + 3*a0 + 2) };
        double p1_1[ 3 ] = { *(xyz1 + 3*a1 + 0), *(xyz1 + 3*a1 + 1), *(xyz1 + 3*a1 + 2) };
        double p1_2[ 3 ] = { *(xyz1 + 3*a2 + 0), *(xyz1 + 3*a2 + 1), *(xyz1 + 3*a2 + 2) };
        double p1_3[ 3 ] = { *(xyz1 + 3*a3 + 0), *(xyz1 + 3*a3 + 1), *(xyz1 + 3*a3 + 2) };

        double p2_0[ 3 ] = { *(xyz2 + 3*a0 + 0), *(xyz2 + 3*a0 + 1), *(xyz2 + 3*a0 + 2) };
        double p2_1[ 3 ] = { *(xyz2 + 3*a1 + 0), *(xyz2 + 3*a1 + 1), *(xyz2 + 3*a1 + 2) };
        double p2_2[ 3 ] = { *(xyz2 + 3*a2 + 0), *(xyz2 + 3*a2 + 1), *(xyz2 + 3*a2 + 2) };
        double p2_3[ 3 ] = { *(xyz2 + 3*a3 + 0), *(xyz2 + 3*a3 + 1), *(xyz2 + 3*a3 + 2) };

        double A11 = MC_ANGLEB(p1_0, p1_1, p1_2);
        double A12 = MC_ANGLEB(p1_1, p1_2, p1_3);
        double A21 = MC_ANGLEB(p2_0, p2_1, p2_2);
        double A22 = MC_ANGLEB(p2_1, p2_2, p2_3);

        if ((fabs(180.0 - A11) < LINEARTHLD) || (fabs(180.0 - A12) < LINEARTHLD) ||
            (fabs(180.0 - A21) < LINEARTHLD) || (fabs(180.0 - A22) < LINEARTHLD)) {
            continue;
        }

        double A1 = MC_DIHEDRAL(p1_0, p1_1, p1_2, p1_3);
        double A2 = MC_DIHEDRAL(p2_0, p2_1, p2_2, p2_3);
        if (isnan(A1) || isnan(A2)) continue;

        double DIFF = A2 - A1;
        if (DIFF > 180.0) DIFF -= 360.0;
        else if (DIFF <= -180.0) DIFF += 360.0;
        DIFF = fabs(DIFF);

        if (DIFF > amax) {
            amax = DIFF;
            if (VER) {
                printf("i:%li   j:%li   k:%li   l:%li   A1:%f  A2:%f DIFF: %f count=%li\n",
                       tors[ i ].A[ 0 ], tors[ i ].A[ 1 ], tors[ i ].A[ 2 ], tors[ i ].A[ 3 ], A1, A2, DIFF, count);
            }
        }
        count++;
    }
    if (AMAX) *AMAX = amax;
}
