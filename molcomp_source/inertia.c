#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "maintypes.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include "atominfo.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

// This function calculates the principal moment of inertia of the molecule and sorts them

void INERTIA (Atom **coord, double II[3], long int NATOMS, double EF[3][3]) {

// calculate the center of mass
/*
	xi = summ(mi*xi)/summ(mi)
*/
double com[3] = {0,0,0};


// calculate the molar mass, M

double M = 0;
long int i=0, j=0, m1=0;
for (i=0;i<NATOMS;i++) {
	M = M + ATOMINFO((*coord)[i].El).M;
}

// calculate the center of mass

for (i=0;i<3;i++) {
	for (j=0;j<NATOMS;j++) {
		com[i] = com[i] + ATOMINFO((*coord)[j].El).M*(*coord)[j].xyz[i];
	}
	com[i] = com[i]/M;
}

// adjust the molecule to the center of mass
//printf("COM %f %f %f \n", com[0], com[1], com[2]);
for (i=0; i<NATOMS; i++) {
  (*coord)[i].xyz[0] = (*coord)[i].xyz[0] - com[0];
  (*coord)[i].xyz[1] = (*coord)[i].xyz[1] - com[1];
  (*coord)[i].xyz[2] = (*coord)[i].xyz[2] - com[2];
}


// MATRIX Iij {{Ixx, Ixy, Ixz},{},{}}
double I[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

for (i=0;i<3;i++) {
	for (j=0;j<3;j++) {
		// calculate the element of Iij
		long int d = 0;
		if (i==j) {
			d = 1;
		}
		else {
			d = 0;
		}
		for (m1=0;m1<NATOMS;m1++) {
			double Mi = ATOMINFO((*coord)[m1].El).M;
			double x = (*coord)[m1].xyz[0];
			double y = (*coord)[m1].xyz[1];
			double z = (*coord)[m1].xyz[2];
			double li = (*coord)[m1].xyz[i];
			double lj = (*coord)[m1].xyz[j];
			I[i][j] = I[i][j] + Mi*(  ( pow(x, 2) + (pow(y, 2) + pow(z,2))  )*d - (li)*(lj) ); 
		} 
	}
}



/*
for (i=0;i<3;i++) {
	for (j=0;j<3;j++) {
		printf ("%f ", I[i][j]);
	}
	printf ("\n");
}
*/
//printf ("com: %f %f %f\n", com[0], com[1], com[2]);


// we need to feed a vector to gsl library 
// I[3][3] -> data[9]
/*
double data[9];
m1 = 0;
for (i=0;i<3;i++) {
	for (j=0;j<3;j++) {
		data[m1] = I[i][j];
		m1++;
	}
}
*/

gsl_matrix *m = gsl_matrix_alloc(3,3);
for (i=0;i<3;i++) {
  for (j=0;j<3;j++) {
    gsl_matrix_set (m, i, j, I[i][j]);
  }
}

// Try to find the principal moments of inertia 
// V,D = linalg.eig (matr_N)
// use function from gsl

// Principal moments of inertia

double PI[3];
//gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3);
gsl_vector *eval = gsl_vector_alloc (3);
//printf("I am here ... \n");
gsl_matrix *evec = gsl_matrix_alloc (3, 3);
//printf("I am here ... \n");
gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  
gsl_eigen_symmv (m, eval, evec, w);


gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

// we need to ensure proper rotation (no mirror images)
gsl_matrix *evectmp = gsl_matrix_alloc (3, 3);

for (i=0;i<3;i++) {
  for (j=0;j<3;j++) {
    gsl_matrix_set (evectmp, i, j, gsl_matrix_get (evec, i, j));
  }
}
gsl_permutation *p = gsl_permutation_alloc(3);
int s;
gsl_linalg_LU_decomp(evectmp, p, &s);
double det = gsl_linalg_LU_det(evectmp, s);
if (VER) printf("Determinant rotation matrix: %f \n", det);
if (round(det) == -1) {
  if (VER) printf ("Improper rotation detected. Multiplication by -1 (last column).\n");
  for (i=0;i<3;i++) {
    gsl_matrix_set (evec, i, 2, gsl_matrix_get (evec, i, 2)*(-1.0));
  }
}

for (i=0;i<3;i++) {
  for (j=0;j<3;j++) {
    gsl_matrix_set (evectmp, i, j, gsl_matrix_get (evec, i, j));
  }
}
gsl_linalg_LU_decomp(evectmp, p, &s);
det = gsl_linalg_LU_det(evectmp, s);
//printf("%f \n", det);
 
{
	long int i;
	for (i = 0; i < 3; i++)
	{
		double eval_i = gsl_vector_get (eval, i);
		//gsl_vector_view evec_i = gsl_matrix_column (evec, i);
		EF[0][i] = gsl_matrix_get (evec, 0, i);
		EF[1][i] = gsl_matrix_get (evec, 1, i);
		EF[2][i] = gsl_matrix_get (evec, 2, i);
		//printf ("%f %f %f --\n", EF[0][i], EF[1][i], EF[2][i]);
		PI[i] = eval_i;
		if (VER) printf ("eigenvalue = %g\n", eval_i);
		//printf ("eigenvector = \n");
		//gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
	}
}
/*
for (i=0;i<3;i++){
  printf ("%f %f %f\n", EF[i][0], EF[i][1], EF[i][2]);
}
*/

for (i=0;i<3;i++) {
	II[i] = PI[i];
}

gsl_vector_free (eval);
gsl_matrix_free (evec);
gsl_matrix_free (evectmp);
gsl_eigen_symmv_free (w);
gsl_permutation_free (p);
}

