#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "mc_types.h"
#include "mc_atominfo.h"
#include "mc_inertia.h"

// Jacobi eigensolver for real symmetric 3x3 matrix (no GSL, no malloc)
static void jacobi_sym3(double A[ 3 ][ 3 ], double eval[ 3 ], double evec[ 3 ][ 3 ]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            evec[ i ][ j ] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < 50; sweep++) {
        double off = fabs(A[ 0 ][ 1 ]) + fabs(A[ 0 ][ 2 ]) + fabs(A[ 1 ][ 2 ]);
        if (off < 1e-15) break;

        for (int p = 0; p < 2; p++) {
            for (int q = p + 1; q < 3; q++) {
                double apq = A[ p ][ q ];
                if (fabs(apq) < 1e-16) continue;

                double app = A[ p ][ p ];
                double aqq = A[ q ][ q ];
                double tau = (aqq - app) / (2.0 * apq);
                double t = (tau >= 0.0)
                    ? (1.0 / (tau + sqrt(1.0 + tau * tau)))
                    : (-1.0 / (-tau + sqrt(1.0 + tau * tau)));
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;

                A[ p ][ p ] = app - t * apq;
                A[ q ][ q ] = aqq + t * apq;
                A[ p ][ q ] = A[ q ][ p ] = 0.0;

                int r = 3 - p - q;
                double arp = A[ r ][ p ];
                double arq = A[ r ][ q ];
                A[ r ][ p ] = A[ p ][ r ] = c * arp - s * arq;
                A[ r ][ q ] = A[ q ][ r ] = s * arp + c * arq;

                for (int k = 0; k < 3; k++) {
                    double vkp = evec[ k ][ p ];
                    double vkq = evec[ k ][ q ];
                    evec[ k ][ p ] = c * vkp - s * vkq;
                    evec[ k ][ q ] = s * vkp + c * vkq;
                }
            }
        }
    }

    eval[ 0 ] = A[ 0 ][ 0 ];
    eval[ 1 ] = A[ 1 ][ 1 ];
    eval[ 2 ] = A[ 2 ][ 2 ];
}

static double det3x3(double M[ 3 ][ 3 ]) {
    return M[ 0 ][ 0 ] * (M[ 1 ][ 1 ] * M[ 2 ][ 2 ] - M[ 1 ][ 2 ] * M[ 2 ][ 1 ])
         - M[ 0 ][ 1 ] * (M[ 1 ][ 0 ] * M[ 2 ][ 2 ] - M[ 1 ][ 2 ] * M[ 2 ][ 0 ])
         + M[ 0 ][ 2 ] * (M[ 1 ][ 0 ] * M[ 2 ][ 1 ] - M[ 1 ][ 1 ] * M[ 2 ][ 0 ]);
}

// This function calculates the principal moments of inertia of the molecule and sorts them
void MC_INERTIA(McAtom **coord, double II[ 3 ], long int NATOMS, double EF[ 3 ][ 3 ], bool VER) {
    // calculate the center of mass: xi = sum(mi*xi)/sum(mi)
    double com[ 3 ] = {0, 0, 0};

    // calculate the molar mass, M
    double M = 0;
    long int i = 0, j = 0, m1 = 0;
    for (i = 0; i < NATOMS; i++) {
    M = M + MC_ATOMINFO((*coord)[ i ].El).M;
    }
    if (M <= 0.0) M = 1.0;

    // calculate the center of mass
    for (i = 0; i < 3; i++) {
        for (j = 0; j < NATOMS; j++) {
            com[ i ] = com[ i ] + MC_ATOMINFO((*coord)[ j ].El).M * (*coord)[ j ].xyz[ i ];
        }
        com[ i ] = com[ i ] / M;
    }

    // adjust the molecule to the center of mass
    for (i = 0; i < NATOMS; i++) {
        (*coord)[ i ].xyz[ 0 ] = (*coord)[ i ].xyz[ 0 ] - com[ 0 ];
        (*coord)[ i ].xyz[ 1 ] = (*coord)[ i ].xyz[ 1 ] - com[ 1 ];
        (*coord)[ i ].xyz[ 2 ] = (*coord)[ i ].xyz[ 2 ] - com[ 2 ];
    }

    // MATRIX Iij {{Ixx, Ixy, Ixz},{},{}}
    double I[ 3 ][ 3 ] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            long int d = (i == j) ? 1 : 0;
            for (m1 = 0; m1 < NATOMS; m1++) {
                double Mi = MC_ATOMINFO((*coord)[ m1 ].El).M;
                double x = (*coord)[ m1 ].xyz[ 0 ];
                double y = (*coord)[ m1 ].xyz[ 1 ];
                double z = (*coord)[ m1 ].xyz[ 2 ];
                double li = (*coord)[ m1 ].xyz[ i ];
                double lj = (*coord)[ m1 ].xyz[ j ];
                I[ i ][ j ] = I[ i ][ j ] + Mi * ((x * x + y * y + z * z) * d - li * lj);
            }
        }
    }

    // Principal moments of inertia and eigenvectors
    double eval[ 3 ];
    double evec[ 3 ][ 3 ];
    jacobi_sym3(I, eval, evec);

    // Sort eigenvalues and eigenvectors in descending order of |eval| (GSL_EIGEN_SORT_ABS_DESC)
    int order[ 3 ] = {0, 1, 2};
    for (int a = 0; a < 2; a++) {
        for (int b = a + 1; b < 3; b++) {
            if (fabs(eval[ order[ b ] ]) > fabs(eval[ order[ a ] ])) {
                int tmp = order[ a ];
                order[ a ] = order[ b ];
                order[ b ] = tmp;
            }
        }
    }

    double PI[ 3 ];
    for (i = 0; i < 3; i++) {
        int src = order[ i ];
        EF[ 0 ][ i ] = evec[ 0 ][ src ];
        EF[ 1 ][ i ] = evec[ 1 ][ src ];
        EF[ 2 ][ i ] = evec[ 2 ][ src ];
        PI[ i ] = eval[ src ];
    }

    // Ensure proper rotation (no mirror images: det(EF) = +1)
    double det = det3x3(EF);
    if (VER) printf("Determinant rotation matrix: %f \n", det);
    if (round(det) == -1) {
        if (VER) printf("Improper rotation detected. Multiplication by -1 (last column).\n");
        for (i = 0; i < 3; i++) {
            EF[ i ][ 2 ] = EF[ i ][ 2 ] * (-1.0);
        }
    }

    for (i = 0; i < 3; i++) {
        if (VER) printf("eigenvalue = %g\n", PI[ i ]);
        II[ i ] = PI[ i ];
    }
}

