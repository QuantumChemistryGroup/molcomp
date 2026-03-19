#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "maintypes.h"
#include <math.h>
#include <string.h>
#include "kabsh.h"
#include "shuffle.h"

void RENUMB (Atom **coord1tmp, Atom *coord2, long int NATOMS) {

  long int i=0;
  long int j=0;
  // temporary array for renumeration
  Atom *coord = malloc(NATOMS*sizeof(Atom));
  // another tmp array to store best numeration
  Atom *coord1best = malloc(NATOMS*sizeof(Atom));
  for (i=0;i<NATOMS;i++) {
    (*coord1tmp)[i].N = -1; // will store its correspondence to sequence in
                            // the second structure
    // second molecule
    coord2[i].N = -1; // will be changed to positive (itself) if already matched 
  }
  long int k=0;
  // summ of maximum interatom distances
  double DMAX = 1e20;
//  double RMSDM = 1e20;
  for (k=0;k<1;k++) {
    long int IT [NATOMS];
    // summ of the distances
    double dd = 0;
    for (i=0;i<NATOMS;i++) {
      IT[i] = i;
    }
    if (k != 0) {
      shuffle (IT, NATOMS);
      }
    long int v=0; 
    for (v=0;v<NATOMS;v++) {
      // take i from shuffled array
      i = IT[v];
      double Dmin=1e10; // store minimum dist from atom i in the first str 
                        // to atom j in the second str 
      long int N = -1;  // corresponding Dmin index from the second str (j)
      for (j=0;j<NATOMS;j++) {
        // if the same element & not previously assigned
        if ( (strcmp((*coord1tmp)[i].El, coord2[j].El) == 0) && (coord2[j].N == -1)  ) {
          double x = (*coord1tmp)[i].xyz[0] - coord2[j].xyz[0];
          double y = (*coord1tmp)[i].xyz[1] - coord2[j].xyz[1];
          double z = (*coord1tmp)[i].xyz[2] - coord2[j].xyz[2];
          double d = sqrt(x*x + y*y + z*z);
          if (d<Dmin) {
            Dmin = d;
            N = j;
          }
        }
      }
      // subscribe new index to 1 mol, designate 2 mol  
      if (N>=0) {
        (*coord1tmp)[i].N = N;
        coord2[N].N = N;
        dd = Dmin + dd; // store the minimum (final) distance for atom pair
      }
      else {
        printf ("N is negative. Stop.\n");
        exit(0);
      }
    }
    // re-populate arrays coord2[j].N
    for (j=0;j<NATOMS;j++) {
      // re-populate arrays coord2[j].N
      coord2[j].N = -1;
    }
//    double RMSD = 0;
    // populate coord (temporary array) according to renumbering
//    if RMSD is used.
//    for (j=0;j<NATOMS;j++) {
//      long int it = (*coord1tmp)[j].N;
//      strcpy(coord[it].El,(*coord1tmp)[j].El);
//      for (long int mm=0;mm<3;mm++) {
//        coord[it].xyz[mm] = (*coord1tmp)[j].xyz[mm];
//      }
      // coord[it].N - will be former number 
//      coord[it].N = j;
//    }
//    Atom *coord555 = malloc(NATOMS*sizeof(Atom));
//    KABSH (coord, coord2, &coord555, NATOMS, &RMSD);
//    free (coord555);
//    printf ("RMSD (in renumb) = %f\n", RMSD);
//    if (RMSD<RMSDM) {
    if (dd<DMAX) {
//      printf ("RMSD (in renumb) = %f %f\n", RMSD, RMSDM);
      for (j=0;j<NATOMS;j++) {
        strcpy(coord1best[j].El, (*coord1tmp)[j].El);
        coord1best[j].xyz[0] = (*coord1tmp)[j].xyz[0];
        coord1best[j].xyz[1] = (*coord1tmp)[j].xyz[1];
        coord1best[j].xyz[2] = (*coord1tmp)[j].xyz[2];
        coord1best[j].N = (*coord1tmp)[j].N;
      }
      DMAX = dd;
//      RMSDM = RMSD;
    }
  //printf ("k= %li DMAX = %f\n", k, DMAX);
  }
  
  // populate coord (temporary array) according to renumbering
  for (i=0;i<NATOMS;i++) {
    long int it = coord1best[i].N;
    strcpy(coord[it].El,coord1best[i].El);
    for (j=0;j<3;j++) {
      coord[it].xyz[j] = coord1best[i].xyz[j];
    }
    // coord[it].N - will be former number 
    coord[it].N = i;
  }
  // overwrite coord1tmp
  for (i=0;i<NATOMS;i++) {
    strcpy((*coord1tmp)[i].El, coord[i].El);
    for (j=0;j<3;j++) {
      (*coord1tmp)[i].xyz[j] = coord[i].xyz[j];
    }
    // coord[it].N - will be former number 
    (*coord1tmp)[i].N = coord[i].N; 
  }

//  double RMSD = 0;
//  KABSH (*coord1tmp, coord2, &coord, NATOMS, RMSD);

/*  
  printf ("First best fit\n");
  for (i=0;i<NATOMS;i++){
    printf("%s %10.6f %10.6f %10.6f %li\n", (*coord1tmp)[i].El, (*coord1tmp)[i].xyz[0], (*coord1tmp)[i].xyz[1], (*coord1tmp)[i].xyz[2], (*coord1tmp)[i].N );
  }
*/  
  free (coord);
  free (coord1best);
}
