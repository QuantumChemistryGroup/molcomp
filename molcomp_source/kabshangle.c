#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "maintypes.h"
#include "stdfunc.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

double DIHEDRAL(double p1[3], double p2[3], double p3[3], double p4[3]);
double ANGLEB (double xyz1[3], double xyz2[3], double xyz3[3]);
void GETVEC(double atom2[3], double atom1[3], double v[3]);
void GETUNITVEC(double v[3], double u[3]);
double scalar_product(double v1[3], double v2[3]);
void vector_product(double u[3], double w[3], double v[3]);


void KABSHANGLE (Atom *coord1, Atom *coord2, Atom **coord1mod, long int NATOMS, double *RMSD, double *AAV, double *AMAX) {
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

// get all dihedrals

MM *tors = malloc(1*sizeof(*tors));
Atom *coord;

if (MOL2 == true) {
  coord = coord2;
}
else {
  coord = coord1;
}

long int l =0;
long int q = 0;
long int m=0;
for (i=0;i<NATOMS;i++) {
  //printf ("%li ", coord[i].nn);
  for (j=0;j<coord[i].nn;j++) {
    long int z = coord[i].n[j];
      for (m=0;m<coord[i].nn;m++) {
        for (q=0;q<coord[z-1].nn;q++) {
          if (coord[i].n[m] ==z || coord[z-1].n[q] ==i+1 ) {
          }
          else {
            //long int A0 = coord[i].n[m];
            //long int A1 = i+1;
            //long int A2 = z;
            //long int A3 = coord[z-1].n[q];
            //bool Copy = false;
            // check if this combination was previously found
            tors[l].A[0] = coord[i].n[m];
            tors[l].A[1] = i+1;
            tors[l].A[2] = z;
            tors[l].A[3] = coord[z-1].n[q];
            //tors[l].V = DIHEDRAL(coord[tors[l].A[0]-1].xyz, coord[tors[l].A[1]-1].xyz, coord[tors[l].A[2]-1].xyz, coord[tors[l].A[3]-1].xyz);
            // reallocate
            l++;
            tors = realloc(tors, (l+1)*sizeof(*tors));
            /*if (VER) {
              printf ("i:%li   j:%li   k:%li   l:%li   V:%f  %li\n", tors[l-1].A[0], tors[l-1].A[1], tors[l-1].A[2], tors[l-1].A[3], tors[l-1].V,  l-1);
            }*/
          }
        }
      }
    }	
}

//printf ("I am here\n");

(*RMSD)=0;

double Aav = 0; // average angle 
double Amax = 0; // maximum
double Dmax = 0; // maximum distance between alligned atoms

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
  if (d>Dmax) {
    Dmax = d;
  }
  //printf("%s %10.6f %10.6f %10.6f\n", coord1[i].El, coord1[i].xyz[0], coord1[i].xyz[1], coord1[i].xyz[2]);
}
long int count=0;
// compare angles
for (i=0;i<l;i++) {
  double A1=0;
  double A2=0;
  A1 = DIHEDRAL(coord1[tors[i].A[0]-1].xyz, coord1[tors[i].A[1]-1].xyz, coord1[tors[i].A[2]-1].xyz, coord1[tors[i].A[3]-1].xyz);
  A2 = DIHEDRAL(coord2[tors[i].A[0]-1].xyz, coord2[tors[i].A[1]-1].xyz, coord2[tors[i].A[2]-1].xyz, coord2[tors[i].A[3]-1].xyz);
  // compare that 1,2,3 or 2,3,4 of A1 or A2 are not close to ~180.0
  
  double A11 = ANGLEB(coord1[tors[i].A[0]-1].xyz, coord1[tors[i].A[1]-1].xyz, coord1[tors[i].A[2]-1].xyz);
  double A12 = ANGLEB(coord1[tors[i].A[1]-1].xyz, coord1[tors[i].A[2]-1].xyz, coord1[tors[i].A[3]-1].xyz);
  double A21 = ANGLEB(coord2[tors[i].A[0]-1].xyz, coord2[tors[i].A[1]-1].xyz, coord2[tors[i].A[2]-1].xyz);
  double A22 = ANGLEB(coord2[tors[i].A[1]-1].xyz, coord2[tors[i].A[2]-1].xyz, coord2[tors[i].A[3]-1].xyz);
  // threshold for linearity
  double thld = LINEARTHLD;
  if ( (fabs(180-A11) < thld) || (fabs(180-A12) < thld) || (fabs(180-A21) < thld) || (fabs(180-A22) < thld)) {
    continue;
  }
  double DIFF=(A2-A1);
  if (DIFF>180.0) {
    DIFF = DIFF - 360;
  }
  else if (DIFF <= -180.0 ) {
    DIFF = DIFF + 360;
  }
  DIFF=fabs(DIFF);
  if ((isnan(A1) == true) || (isnan(A2) == true)) {
    continue;
  }
  Aav = Aav + DIFF;
  if (DIFF>Amax) {
    Amax = DIFF;
  if (VER) {
    printf ("i:%li   j:%li   k:%li   l:%li   A1:%f  A2:%f DIFF: %f count=%li\n", tors[i].A[0], tors[i].A[1], tors[i].A[2], tors[i].A[3], A1,  A2, DIFF, count);
  }
  }
  count = count + 1;
}
if (count !=0) {
  Aav = Aav/count;
}
else {
  Aav=0;
}
(*RMSD) = sqrt ((1.0/NATOMS)*(*RMSD));
if (VER) printf ("RMSD = %f\n", (*RMSD));
//(*AAV) = Aav;
(*AAV) = Dmax; // now the value in the middle is the largest distance
(*AMAX) = Amax;

//printf (" Aav= %f; Amax= %f \n", Aav, Amax);

// free
coord=NULL;
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

/*

Given the coordinates of the four points, obtain the vectors b1, b2, and b3

by vector subtraction.

Let me use the nonstandard notation ⟨v⟩
to denote v/∥v∥, the unit vector in the direction of the vector v. 
Compute n1=⟨b1×b2⟩ and n2=⟨b2×b3⟩, the normal vectors to the planes
containing b1 and b2, and b2 and b3 respectively. The angle we seek
is the same as the angle between n1 and n2

.

The three vectors n1
, ⟨b2⟩, and m1:=n1×⟨b2⟩ form an orthonormal frame. Compute the
coordinates of n2 in this frame: x=n1⋅n2 and y=m1⋅n2.
(You don't need to compute ⟨b2⟩⋅n2

as it should always be zero.)

The dihedral angle, with the correct sign, is atan2(y,x)
.

*/

double DIHEDRAL(double p1[3], double p2[3], double p3[3], double p4[3]){
double angle;

double b1[3], b2[3], b3[3];
GETVEC(p2,p1,b1);
GETVEC(p3,p2,b2);
GETVEC(p4,p3,b3);

double b1xb2[3];
double b2xb3[3];

vector_product(b1,b2,b1xb2);
vector_product(b2,b3,b2xb3);

double n1[3], n2[3], b2x[3], m1[3];

GETUNITVEC(b1xb2, n1);
GETUNITVEC(b2xb3, n2);
GETUNITVEC(b2, b2x);

vector_product(n1,b2x,m1);

double x,y;
x=scalar_product(n1,n2);
y=scalar_product(m1,n2);

angle=atan2(y,x);

angle=angle*180/M_PI;

return angle;
}

void GETVEC(double atom2[3], double atom1[3], double v[3]) {

int i;
for (i=0;i<3;i++) {
        v[i] = *(atom2+i) - *(atom1+i);
}

}

void GETUNITVEC(double v[3], double u[3]) {

int i;
// get the vector length
double lenv= sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
// unit vector u
for (i=0;i<3;i++) {
        u[i] = v[i]/lenv;
}
// printf ("Unit vector %f %f %f\n", u[0], u[1], u[2]);

}

double scalar_product(double v1[3], double v2[3]){
double a1;
a1 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
return a1;
}

void vector_product(double u[3], double w[3], double v[3]) {

v[0] = u[1]*w[2] - u[2]*w[1];                             
v[1] = u[2]*w[0] - u[0]*w[2];                             
v[2] = u[0]*w[1] - u[1]*w[0];                             

}

double ANGLEB (double xyz1[3], double xyz2[3], double xyz3[3]) {
double angle= 0;

// get the vectors
double v1[3] = {xyz1[0]-xyz2[0], xyz1[1]-xyz2[1], xyz1[2]-xyz2[2]};
double v2[3] = {xyz3[0]-xyz2[0], xyz3[1]-xyz2[1], xyz3[2]-xyz2[2]};

double v1mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
double v1norm[3] = {v1[0]/v1mag, v1[1]/v1mag, v1[2]/v1mag};

double v2mag = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
double v2norm[3] = {v2[0]/v2mag, v2[1]/v2mag, v2[2]/v2mag};

// calculate the dot product

double res = v1norm[0]*v2norm[0] + v1norm[1]*v2norm[1] + v1norm[2]*v2norm[2];

angle = acos(res)*(180/M_PI);

return angle;
}

