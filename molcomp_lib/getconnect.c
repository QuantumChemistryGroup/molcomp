#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "maintypes.h"
#include <math.h>
#include <stdlib.h>
#include "atominfo.h"


bool isvalueinarray(long int val, long int *arr, long int size, long int start){
        bool Belong = false;
        long int i;
        for (i=start; i < size; i++) {
                if (arr[i] == val){
                        Belong = true;
                }       
        }

        return Belong;
}

double ANGLE (double xyz1[3], double xyz2[3], double xyz3[3]) {
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

void GETCONNECT (Atom **coordinates, long int NATOMS) {

long int i=0, j=0, i1=0;
long int maxbond = 0;
long int *bonds;
// array with the distances corresponding to atom numbers in bonds
// we need it to rule out the distances w
double *dist;
// distance
double d;

for (i=0;i<NATOMS;i++) {
	maxbond = 0;
	bonds = malloc(1*sizeof(long int));
	dist = malloc(1*sizeof(double));
//	printf("%s %f %s...\n", (*coordinates)[i].El, ATOMINFO((*coordinates)[i].El).R, ATOMINFO((*coordinates)[i].El).ENAME);
	for (j=0;j<NATOMS;j++) {
		// skip comparison with itself
		if (i != j) {
			// calculate the distance
			d = sqrt (pow(((*coordinates)[i].xyz[0] - (*coordinates)[j].xyz[0]),2)
		   +  pow(((*coordinates)[i].xyz[1] - (*coordinates)[j].xyz[1]),2)
		    + pow(((*coordinates)[i].xyz[2] - (*coordinates)[j].xyz[2]),2) );
		        // check whether the distance is less than sum of covalent radii + 0.45A
 			if (d < ATOMINFO((*coordinates)[i].El).R + ATOMINFO((*coordinates)[j].El).R + 0.45  ) {
				// j+1 since the numeration starts from 0
 				bonds[maxbond] = j+1;
 				dist[maxbond] = d;
				maxbond++;
 				bonds = (long int *) realloc(bonds, (maxbond+1)*sizeof(long int));
 				dist = (double *) realloc(dist, (maxbond+1)*sizeof(double));
// 				printf ("%i %i \n", maxbond, bonds[maxbond-1]);
	
 			} 
		    
		}
	
	}
	// check with maximum valency
	if (maxbond <= ATOMINFO((*coordinates)[i].El).V ) {
	// everything is perfect
	// organize copying
		(*coordinates)[i].nn=maxbond;
		for (j=0;j<maxbond;j++) {
			(*coordinates)[i].n[j]=bonds[j];
		}
	}
	else {
	// we have to delete the largest distances
	// buble sorting
	// then the largest distances are removed.
	// Another algorithm would be to run a loop over connected atoms
	// and detelete there which already have the saturated valcies
	// and keep the atoms which do hae the unsaturated valencies
//	printf ("Something is wrong\n");
		for (j=0;j<maxbond;j++) {
			for (i1=0;i1<(maxbond-j-1);i1++) {
				if (dist[i1] > dist[i1+1]) {
					long int N;
					double dist1;
					dist1 = dist[i1];
					N = bonds[i1];
					dist[i1] = dist[i1+1];
					bonds[i1] = bonds[i1+1];
					dist[i1+1] = dist1;
					bonds[i1+1] = N;	
				}
			}
		}

	// organize copying
		(*coordinates)[i].nn=ATOMINFO((*coordinates)[i].El).V;
		for (j=0;j<(*coordinates)[i].nn;j++) {
			(*coordinates)[i].n[j]=bonds[j];
		}
	
	}
	
	free (bonds);
	free (dist);	
		
}
//printf ("PART 1 is done...\n");

// get rid of some bonds which give the angle less than 50
for (i=0;i<NATOMS;i++) {
	long int pairs[1000][2];
	double angle;
	long int count=0;
	for (i1=0;i1<(*coordinates)[i].nn;i1++) {
		for (j=i1+1;j<(*coordinates)[i].nn;j++) {
			angle = ANGLE((*coordinates)[ (*coordinates)[i].n[i1]-1].xyz, (*coordinates)[i].xyz, (*coordinates)[ (*coordinates)[i].n[j]-1].xyz);
			if (angle < 50.0) {
//					printf ("count %i ... %i %i %i\n", count, (*coordinates)[i].n[i1], i+1, (*coordinates)[i].n[j]);
					// if angle is less then 60 - rememeber the atom numbers.
					pairs[count][0] = (*coordinates)[i].n[i1];
					pairs[count][1] = (*coordinates)[i].n[j];		
					// calculate bond 1 and criteria 1
					long int z = pairs[count][0];
					long int z1 = pairs[count][1];
					d = sqrt (pow(((*coordinates)[i].xyz[0] - (*coordinates)[z-1].xyz[0]),2)
					   +  pow(((*coordinates)[i].xyz[1] - (*coordinates)[z-1].xyz[1]),2)
					    + pow(((*coordinates)[i].xyz[2] - (*coordinates)[z-1].xyz[2]),2) );
					// creteria to check the distance	
					double cr = ATOMINFO((*coordinates)[i].El).R + ATOMINFO((*coordinates)[z-1].El).R + 0.45;
					// calculate bond 2 and criteria 2
					double d1 = sqrt (pow(((*coordinates)[i].xyz[0] - (*coordinates)[z1-1].xyz[0]),2)
					   +  pow(((*coordinates)[i].xyz[1] - (*coordinates)[z1-1].xyz[1]),2)
					    + pow(((*coordinates)[i].xyz[2] - (*coordinates)[z1-1].xyz[2]),2) );
					// creteria to check the distance	
					double cr1 = ATOMINFO((*coordinates)[i].El).R + ATOMINFO((*coordinates)[z1-1].El).R + 0.45;
					if ((cr-d) < (cr1-d1)) {
						pairs[count][0] = z;
					}
					else {
						pairs[count][0] = z1;
					}
//					printf ("count %i ... %i %i %i\n", count, (*coordinates)[i].n[i1], i+1, (*coordinates)[i].n[j]);
					count++;
			}
		}
	}
//	printf ("PART 2 is done... atom %i\n", i+1);
	// calculated all the angles for atom
	// run over the all strange angles and find the atom to be disconnected
	// unless both atoms forming the strange angle are in freeze group
/*	for (i1=0;i1<count;i1++) {
		// calculate bond 1 and criteria 1
		long int z = pairs[i1][0];
		long int z1 = pairs[i1][1];
		d = sqrt (pow(((*coordinates)[i].xyz[0] - (*coordinates)[z-1].xyz[0]),2)
		   +  pow(((*coordinates)[i].xyz[1] - (*coordinates)[z-1].xyz[1]),2)
		    + pow(((*coordinates)[i].xyz[2] - (*coordinates)[z-1].xyz[2]),2) );
		// creteria to check the distance	
		double cr = ATOMINFO((*coordinates)[i].El).R + ATOMINFO((*coordinates)[z-1].El).R + 0.45;
		// calculate bond 2 and criteria 2
		double d1 = sqrt (pow(((*coordinates)[i].xyz[0] - (*coordinates)[z1-1].xyz[0]),2)
		   +  pow(((*coordinates)[i].xyz[1] - (*coordinates)[z1-1].xyz[1]),2)
		    + pow(((*coordinates)[i].xyz[2] - (*coordinates)[z1-1].xyz[2]),2) );
		// creteria to check the distance	
		double cr1 = ATOMINFO((*coordinates)[i].El).R + ATOMINFO((*coordinates)[z1-1].El).R + 0.45;
		if ((cr-d) < (cr1-d1)) {
			pairs[i1][0] = z;
		}
		else {
			pairs[i1][0] = z1;
		}
		// now pairs[i1][0] is the element to be disconnected

	}*/
//	printf ("PART 3 is done... atom %i\n", i+1);

	// detelete the atoms
	
	// allocate the new array
//	long int atoms[(*coordinates)[i].nn-count];
	long int atoms[(*coordinates)[i].nn];
	long int count1=0;
	for (i1=0;i1<(*coordinates)[i].nn;i1++) {
		bool IN = false;
		for (j=0;j<count;j++) {
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
//	printf ("PART 4 is done... atom %i\n", i+1);

	// replace the current atom connections
	(*coordinates)[i].nn = count1;
	for (i1=0;i1<(*coordinates)[i].nn;i1++) {
		(*coordinates)[i].n[i1] = atoms[i1];
	}	
//	printf ("PART 5 is done... atom %i count %i\n", i+1, count);
	
	// de-connect the atoms which were connected to these element
	for (i1=0;i1<count;i1++) {
//		printf("This is %i", pairs[i1][0]);
		// check the atom
		long int z;
		z = pairs[i1][0];
//		printf ("This is i+1 %i this is z %i\n", i+1, z);
		if (isvalueinarray(i+1, (*coordinates)[z-1].n, (*coordinates)[z-1].nn, 0)==true) {
//			printf("YES\n");
			// there is an element which has to be removed
			long int atoms1[(*coordinates)[z-1].nn-1];
			long int count1 = 0;
			for (j=0;j<(*coordinates)[z-1].nn;j++) {
				if ( (*coordinates)[z-1].n[j] != i+1 ) {
					atoms1[count1] = (*coordinates)[z-1].n[j];
					count1++;
				}	
			}
			// change the existing bonding info of the atom
			(*coordinates)[z-1].nn = (*coordinates)[z-1].nn - 1;
			for (j=0;j<(*coordinates)[z-1].nn;j++) {
				(*coordinates)[z-1].n[j] = atoms1[j];
			}
		
		}
	}
}

// check symmetry, i.e. if 1 connected to 2, then 2 should be connected to 1

for (i=0;i<NATOMS;i++) {
	long int N=0;
	long int DEL[23];
	for (j=0;j<(*coordinates)[i].nn;j++) {
		long int ind = (*coordinates)[i].n[j]-1;
		//printf("%li ;", ind);
		bool F=false;
		for (i1=0;i1<(*coordinates)[ind].nn;i1++) {
			if ((i+1) == (*coordinates)[ind].n[i1]) {
				F=true;
				break;
			}
		}
		//printf("F=%i %li ", F,i1);
		// if there is no element i in ind
		if (F==false) {
			DEL[N]=ind+1;
			N=N+1;
		}
	}
	if (N>0) {
		// modification is needed
		// printf("Cleaning up bonding.\n");
		long int nn = (*coordinates)[i].nn-N;
		long int tmp[nn];
		long int x=0;
		for (j=0;j<(*coordinates)[i].nn;j++) {
			bool Found=isvalueinarray((*coordinates)[i].n[j], DEL, N, 0);
			if (Found==false) {
				tmp[x] = (*coordinates)[i].n[j];
				x=x+1;
			}
		}		
		(*coordinates)[i].nn=(*coordinates)[i].nn-N;
		for (j=0;j<x;j++) {
			(*coordinates)[i].n[j]=tmp[j];
		}
		
	}
	//printf ("\n");
}


}


void SYMMETRYCHECK (Atom **coordinates, long int NATOMS) {

// check symmetry, i.e. if 1 connected to 2, then 2 should be connected to 1
long int i=0, j=0, i1=0;
for (i=0;i<NATOMS;i++) {
	long int N=0;
	long int DEL[23];
	for (j=0;j<(*coordinates)[i].nn;j++) {
		long int ind = (*coordinates)[i].n[j]-1;
		//printf("%li ;", ind);
		bool F=false;
		for (i1=0;i1<(*coordinates)[ind].nn;i1++) {
			if ((i+1) == (*coordinates)[ind].n[i1]) {
				F=true;
				break;
			}
		}
		//printf("F=%i %li ", F,i1);
		// if there is no element i in ind
		if (F==false) {
			DEL[N]=ind+1;
			N=N+1;
		}
	}
	if (N>0) {
		// modification is needed
		//printf("Cleaning up bonding.\n");
		long int nn = (*coordinates)[i].nn-N;
		long int tmp[nn];
		long int x=0;
		for (j=0;j<(*coordinates)[i].nn;j++) {
			bool Found=isvalueinarray((*coordinates)[i].n[j], DEL, N, 0);
			if (Found==false) {
				tmp[x] = (*coordinates)[i].n[j];
				x=x+1;
			}
		}		
		(*coordinates)[i].nn=(*coordinates)[i].nn-N;
		for (j=0;j<x;j++) {
			(*coordinates)[i].n[j]=tmp[j];
		}
		
	}
	//printf ("\n");
}

}

void SYMMETRYCHECKDEF (Atom *coorddef, Atom **coord1, long int NATOMS) {
  // start comparison
  long int i=0;
  long int j=0;
  long int k=0;
  for (i=0;i<NATOMS;i++) {
    // check if there are any new elements bound compared to previous
    // if yes, then to ensure that new elements are bound to this one
    for (j=0;j<(*coord1)[i].nn;j++) {
      // F - found 
      bool F = true;
      F = isvalueinarray((*coord1)[i].n[j], coorddef[i].n, coorddef[i].nn, 0);
      //printf ("%i %li\n", F, i);
      if (F == false) {
        long int m = (*coord1)[i].n[j]-1;
        //printf ("%li \n", i);
        // go to this element & check & fix if needed
        bool F1 = true;
        F1 = isvalueinarray(i+1, (*coord1)[m].n, (*coord1)[m].nn, 0);
        if (F1 == false) {
          (*coord1)[m].nn = (*coord1)[m].nn + 1;
          long int n = (*coord1)[m].nn;
          (*coord1)[m].n[n-1] = i+1; 
        }
      }
    }
    // check if some elements were deleted compared to previous
    // if yes, ensure that these bonds were removed
    for (j=0;j<(coorddef)[i].nn;j++) {
      // F - found 
      bool F = true;
      F = isvalueinarray((coorddef)[i].n[j], (*coord1)[i].n, (*coord1)[i].nn, 0);
      if (F == false) {
        long int m = (coorddef)[i].n[j]-1;
        // go to this element & check & fix if needed
        bool F1 = true;
        F1 = isvalueinarray(i+1, (*coord1)[m].n, (*coord1)[m].nn, 0);
        if (F1 == true) {
          /*
          for (k=0;k<(*coord1)[m].nn;k++) {
            printf ("%li ", (*coord1)[m].n[k]);
          }*/
          //printf("--\n");
          (*coord1)[m].nn = (*coord1)[m].nn - 1;
          long int n = (*coord1)[m].nn;
          long int new[n];
          long int count = 0;
          // go through previous array of bonded
          // and do not copy elements that has to be deleted
          for (k=0;k<(*coord1)[m].nn+1;k++) {
            if ((*coord1)[m].n[k] != i+1) {
              new[count] = (*coord1)[m].n[k];
              count = count+1;
            }
          }
          // copy new to (*coord1)[m].n
          for (k=0;k<(*coord1)[m].nn;k++) {
              (*coord1)[m].n[k] = new[k];
            }
          /* 
          for (k=0;k<(*coord1)[m].nn;k++) {
            printf ("%li ", (*coord1)[m].n[k]);
          }*/
          //printf("-\n");
        }
      }
    }
  }
}
