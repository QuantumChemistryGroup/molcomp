#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "maintypes.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

void KABSHDIST (Atom *coord1, Atom *coord2, Atom **coord1mod, long int NATOMS, double *RMSD, double *DAV, double *DMAX) {
long int i=0;
long int j=0;

// keep unperturbed structures
for (i=0; i<NATOMS; i++) {
  strcpy((*coord1mod)[i].El, coord1[i].El);
  for (j=0;j<3;j++) {
    (*coord1mod)[i].xyz[j] = coord1[i].xyz[j];
  }
  (*coord1mod)[i].N = coord1[i].N;
}

// create two coordinate matrixes P (coord1) and Q (coord2)
// P (Nx3) and Q(N3)
gsl_matrix *P = gsl_matrix_alloc (NATOMS, 3);
gsl_matrix *Q = gsl_matrix_alloc (NATOMS, 3);
gsl_matrix *H = gsl_matrix_alloc (3, 3);

// calculate centroid

double cent1[3] = {0,0,0};
double cent2[3] = {0,0,0};

for (i=0; i<NATOMS; i++) {
  for (j=0;j<3;j++) {
    cent1[j] = coord1[i].xyz[j] + cent1[j];
    cent2[j] = coord2[i].xyz[j] + cent2[j];   
  }
}

for (i=0;i<3;i++) {
  cent1[i] = cent1[i]/NATOMS;
  cent2[i] = cent2[i]/NATOMS;
}

//printf ("centroid 1 %f %f %f \n", cent1[0], cent1[1], cent1[2]);
//printf ("centroid 2 %f %f %f \n", cent2[0], cent2[1], cent2[2]);

// populate matrix

for (i=0; i<NATOMS; i++) {
  for (j=0;j<3;j++) {
    gsl_matrix_set (P, i, j, coord1[i].xyz[j]-cent1[j]);
    gsl_matrix_set (Q, i, j, coord2[i].xyz[j]-cent2[j]);
  }
}


// transpose first matrix
// P(3xN)
//gsl_matrix_transpose (P);

// calculate H = (P^T)Q
gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, Q, 0.0, H);

// Calculate SVD: 
// H = US(V^T)

gsl_vector *S = gsl_vector_alloc(3);
gsl_vector *w = gsl_vector_alloc(3);
gsl_matrix *V = gsl_matrix_alloc(3,3);
gsl_matrix *X = gsl_matrix_alloc(3,3);


//printf ("Original matrix H:\n");
//gsl_matrix_fprintf(stdout, H, "%g");

int status = gsl_linalg_SV_decomp (H, V, S, w);
//int status = gsl_linalg_SV_decomp_jacobi (H, V, S);
if (status != GSL_SUCCESS) {
  fprintf(stderr, "GSL SVD decomposition failed: %s\n", gsl_strerror(status));
}
/*
printf ("U (stored in H\n");
gsl_matrix_fprintf(stdout, H, "%g");
printf ("Singular Values (S)\n");
gsl_vector_fprintf(stdout, S, "%g");
printf("V:\n");
gsl_matrix_fprintf(stdout, V, "%g");
//printf("Note: To get the original matrix back, compute U * S_matrix * V^T\n");
*/
double d=0;
// finding determinants of U and V
// U

gsl_matrix *Utmp = gsl_matrix_alloc (3, 3);
gsl_matrix *Vtmp = gsl_matrix_alloc (3, 3);

for (i=0;i<3;i++) {
  for (j=0;j<3;j++) {
    gsl_matrix_set (Utmp, i, j, gsl_matrix_get (H, i, j));
    gsl_matrix_set (Vtmp, i, j, gsl_matrix_get (V, i, j));
  }
}

gsl_permutation *p = gsl_permutation_alloc(3);
int s;
gsl_linalg_LU_decomp(Utmp, p, &s);
double detU = gsl_linalg_LU_det(Utmp, s);
// V
gsl_permutation *p1 = gsl_permutation_alloc(3);
int s1;
gsl_linalg_LU_decomp(Vtmp, p1, &s1);
double detV = gsl_linalg_LU_det(Vtmp, s1);
// d = det(U)*det(V)
d = detU*detV;

gsl_matrix *M = gsl_matrix_alloc (3, 3);
for (i=0; i<3; i++) {
  for (j=0;j<3;j++) {
    if (i!=j) {
      gsl_matrix_set (M, i, j, 0);
    }
    else if (i<2 && j< 2) {
      gsl_matrix_set (M, i, j, 1);
    }
    else {
      gsl_matrix_set (M, i, j, d);
    }
  } 
}

//printf("M:\n");
//gsl_matrix_fprintf(stdout, M, "%g");
//printf("%f %f %f\n", detU, detV, d);


gsl_matrix *RTMP = gsl_matrix_alloc (3, 3);
gsl_matrix *R = gsl_matrix_alloc (3, 3);

gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, M, V, 0.0, RTMP);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H, RTMP, 0.0, R);
//gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H, V, 0.0, R);
//printf("ROT:\n");
//gsl_matrix_fprintf(stdout, R, "%g");
gsl_matrix *PROT = gsl_matrix_alloc (NATOMS, 3);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P, R, 0.0, PROT);
//printf("PROT:\n");
//gsl_matrix_fprintf(stdout, PROT, "%g");

(*RMSD)=0;

double Dav = 0; // average distance between atoms of two molecules
double Dmax = 0; // maximum distance between atoms of two molecules

for (i=0;i<NATOMS;i++) {
  double d = 0; // distance between atom pair
  for (j=0;j<3;j++) {
    (*coord1mod)[i].xyz[j] = gsl_matrix_get (PROT, i, j);
    double x0 = gsl_matrix_get (PROT, i, j);
    double x1 = gsl_matrix_get (Q, i, j);
    double dx2 = (x0-x1)*(x0-x1);
    (*RMSD) = dx2 + (*RMSD);
    d = dx2 + d;
  }
  d = sqrt (d);
  Dav = d + Dav;
  if (d>Dmax) {
    Dmax = d;
  }
  //printf("%s %10.6f %10.6f %10.6f\n", coord1[i].El, coord1[i].xyz[0], coord1[i].xyz[1], coord1[i].xyz[2]);
}
Dav = Dav/NATOMS;
(*RMSD) = sqrt ((1.0/NATOMS)*(*RMSD));
if (VER) printf ("RMSD = %f\n", (*RMSD));
(*DAV) = Dav;
(*DMAX) = Dmax;
/*
if (Dmax/Dav > 8.0) {
  printf (" Warning: Dav= %f; Dmax= %f \n", Dav, Dmax);
}
*/
// free
gsl_matrix_free (P);
gsl_matrix_free (Q);
gsl_matrix_free (Utmp);
gsl_matrix_free (Vtmp);
gsl_matrix_free (H);
gsl_matrix_free (V);
gsl_matrix_free (X);
gsl_matrix_free (RTMP);
gsl_matrix_free (R);
gsl_matrix_free (M);
gsl_matrix_free (PROT);
gsl_vector_free (S);
gsl_vector_free (w);
gsl_permutation_free(p);
gsl_permutation_free(p1);
}
