#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include "maintypes.h"
#include <math.h>
#include <string.h>
#include "kabsh.h"
#include "hungarian.h"
#include "getconnect.h"
#include "stdfunc.h"


void RENUMBHUNGCONNECT (Atom **coord1tmp, Atom *coord2, long int NATOMS, char (*UnEl)[3], long int UE, bool IDN[][NATOMS]) {

  long int i=0;
  long int j=0;
  // temporary array for renumeration
  Atom *coord = malloc(NATOMS*sizeof(Atom));
  // another tmp array to store best numeration
  Atom *coord1best = malloc(NATOMS*sizeof(Atom));
  // array with the number of each element
  long int UnElN[UE];
  for (i=0;i<UE;i++) {
    UnElN[i] = 0;
  }
  for (i=0;i<NATOMS;i++) {
    (*coord1tmp)[i].N = -1; // will store its correspondence to sequence in
    for (j=0;j<UE;j++) {
      if (strcmp(coord2[i].El,UnEl[j]) == 0) {
        UnElN[j] = UnElN[j]+1;
      }
    }                        // the second structure
    // second molecule
    coord2[i].N = -1; // will be changed to positive (itself) if already matched 
  }
  // check if any of unique element is polyvalent
  //printf ("Start renumb\n");
  bool POLYVAL[UE];
  for (i=0;i<UE;i++) {
    POLYVAL[i] = false;
    for (j=0;j<NATOMS;j++) {
      if (strcmp(UnEl[i], (*coord1tmp)[j].El) ==0 ) {
        if ((*coord1tmp)[j].nn != 1) {
          POLYVAL[i] = true;
          break;  
        }
      }      
    }
  }
  for (i=0;i<UE;i++) {
    //printf ("%s %i\n", UnEl[i], POLYVAL[i]);
  }
  /*
  for (i=0;i<NATOMS;i++){
    printf("%3li %2s %10.6f %10.6f %10.6f %8.3f %6s ", i+1, (*coord1tmp)[i].El, (*coord1tmp)[i].xyz[0], (*coord1tmp)[i].xyz[1], (*coord1tmp)[i].xyz[2], (*coord1tmp)[i].Z, (*coord1tmp)[i].tu);
    for (j=0;j<(*coord1tmp)[i].nn;j++) {
      printf("%3li ", (*coord1tmp)[i].n[j]);              
    }
  printf("\n");
  }
  */
//  for (i=0;i<UE;i++) {
//    printf ("%li ", UnElN[i]);
//  }
//  printf("\n");
  // for each element we construc the cost matrix
  for (i=0;i<UE;i++) {
    // for polyvalent atoms first
    if ( POLYVAL[i] != false ) {
      long long int IND = UnElN[i];
      //long int COST[IND][IND];
      long long int (*COST)[IND] = malloc(IND*sizeof(*COST));
      long int p=0;
      long int q=0;
      long int m=0;
      long int l=0;
      //printf ("Start %li %li\n", i, IND);
      // populate matrix with interatomic distances
      for (p=0;p<NATOMS;p++) {
        if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
          // start index for i-th El
          // iterate via another molecule
          for (q=0;q<NATOMS;q++) {
            if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
              // start index for i-th El in the second mol
              // compute the distance
              //printf ("m l %li %li\n", m, l);
              double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
              double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
              double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
              // in cost matrix we cannot operate with doubles
              //printf ("%li - %li %i \n", p+1, q+1, BONDEQUAL(p, q, (*coord1tmp), NATOMS));
              long long int d = sqrt(x*x + y*y + z*z)*1e6;
              if (IDN[p][q]== true ) {
                d = d/1000;
              }
              COST[m][l] = d;
              l=l+1;
            }
          }
        l=0; 
        m=m+1;
        }
      }
      /*
      printf ("Cost matrix (%lli %lli):\n", IND, IND);
      for (p=0;p<IND;p++) {
        for (q=0;q<IND;q++) {
          printf ("%lli ", COST[p][q]);
        }
      printf ("\n");
      }
      */
      // array storing matched elements in mol2
      //long int MATCH[IND];
      long int *MATCH = malloc(IND*sizeof(long int));  
      //printf ("HUNGMAIN START\n");
      HUNGMAIN(IND, IND, COST, &MATCH);
      //printf ("HUNGMAIN FINISH\n");
      /*
      for (p=0;p<IND;p++) {
        printf("%li ", MATCH[p]);
      }
      printf("\n");
      */
      m=0;
      l=0;
      for (p=0;p<NATOMS;p++) {
        if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
          // start index for i-th El
          // iterate via another molecule
        
          for (q=0;q<NATOMS;q++) {
            if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
              if (MATCH[m] == l) {
                // matched
                (*coord1tmp)[p].N = q;
                if (p != q) {
                  // reodering takes place
                  //printf ("%li %li -- pq\n", p, q);
                }
              }
              l=l+1;
            }
          }
        // assigment
        l=0; 
        m=m+1;
        }
      }
    free (MATCH);
    free (COST);
    }
  }
  // deal with univalent atoms
  for (i=0;i<UE;i++) {
    // for polyvalent atoms first
    if ( POLYVAL[i] == false ) {
      long long int IND = UnElN[i];
      //long int COST[IND][IND];
      long long int (*COST)[IND] = malloc(IND*sizeof(*COST));
      long int p=0;
      long int q=0;
      long int m=0;
      long int l=0;
      //printf ("Start %li %lli\n", i, IND);
      // populate matrix with interatomic distances
      for (p=0;p<NATOMS;p++) {
        if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
          // start index for i-th El
          // iterate via another molecule
          for (q=0;q<NATOMS;q++) {
            if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
              // start index for i-th El in the second mol
              // compute the distance
              //printf ("m l %li %li\n", m, l);
              double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
              double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
              double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
              // in cost matrix we cannot operate with doubles
              long long int d = sqrt(x*x + y*y + z*z)*1e6;
              // check if not connected to the same atom
              long int ino = (*coord1tmp)[p].n[0]-1; // old atom this is connected to
              long int inw = (*coord1tmp)[ino].N;  // old atom new index
//              if ((ino == inw) && (coord2[q].n[0]-1 == inw) ) {

              if (MOL2 == false) {
                if ( ((*coord1tmp)[q].n[0]-1 == inw) && (IDN[ino][inw]== true) ) {
                  d = d/1e3;
                  //printf ("%li %li %li %li %li -- \n", p, q, ino, inw, d);
                }
              }
              else {
                if ( (coord2[q].n[0]-1 == inw) && (IDN[ino][inw]== true ) ) {
                  d = d/1e3;
                  //printf ("%li %li %li %li %li -- \n", p, q, ino, inw, d);
                }
              }
              COST[m][l] = d;
              l=l+1;
            }
          }
        l=0; 
        m=m+1;
        }
      }
      /*
      printf ("Cost matrix (%lli %lli):\n", IND, IND);
      for (p=0;p<IND;p++) {
        for (q=0;q<IND;q++) {
          printf ("%lli ", COST[p][q]);
        }
      printf ("\n");
      }
      */
      // array storing matched elements in mol2
      //long int MATCH[IND];
      long int *MATCH = malloc(IND*sizeof(long int));  
      //printf ("HUNGMAIN START\n");
      HUNGMAIN(IND, IND, COST, &MATCH);
      //printf ("HUNGMAIN FINISH\n");
      
      for (p=0;p<IND;p++) {
        //printf("%li ", MATCH[p]);
      }
      //printf("\n");
      
      m=0;
      l=0;
      for (p=0;p<NATOMS;p++) {
        if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
          // start index for i-th El
          // iterate via another molecule
        
          for (q=0;q<NATOMS;q++) {
            if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
              if (MATCH[m] == l) {
                // matched
                (*coord1tmp)[p].N = q;
              }
              l=l+1;
            }
          }
        // assigment
        l=0; 
        m=m+1;
        }
      }
    free (MATCH);
    free (COST);
    }
  }
  
  
  // populate coord (temporary array) according to renumbering
  for (i=0;i<NATOMS;i++) {
    long int it = (*coord1tmp)[i].N;
    strcpy(coord[it].El,(*coord1tmp)[i].El);
    for (j=0;j<3;j++) {
      coord[it].xyz[j] = (*coord1tmp)[i].xyz[j];
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

void RENUMBHUNG (Atom **coord1tmp, Atom *coord2, long int NATOMS, char (*UnEl)[3], long int UE) {

  long int i=0;
  long int j=0;
  // temporary array for renumeration
  Atom *coord = malloc(NATOMS*sizeof(Atom));
  // another tmp array to store best numeration
  Atom *coord1best = malloc(NATOMS*sizeof(Atom));
  // array with the number of each element
  long int UnElN[UE];
  for (i=0;i<UE;i++) {
    UnElN[i] = 0;
  }
  for (i=0;i<NATOMS;i++) {
    (*coord1tmp)[i].N = -1; // will store its correspondence to sequence in
    for (j=0;j<UE;j++) {
      if (strcmp(coord2[i].El,UnEl[j]) == 0) {
        UnElN[j] = UnElN[j]+1;
      }
    }                        // the second structure
    // second molecule
    coord2[i].N = -1; // will be changed to positive (itself) if already matched 
  }
//  for (i=0;i<UE;i++) {
//    printf ("%li ", UnElN[i]);
//  }
//  printf("\n");
  // for each element we construc the cost matrix
  for (i=0;i<UE;i++) {
    long long int IND = UnElN[i];
    //long int COST[IND][IND];
    long long int (*COST)[IND] = malloc(IND*sizeof(*COST));
    long int p=0;
    long int q=0;
    long int m=0;
    long int l=0;
    //printf ("Start %li %li\n", i, IND);
    // populate matrix with interatomic distances
    for (p=0;p<NATOMS;p++) {
      if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
        // start index for i-th El
        // iterate via another molecule
        for (q=0;q<NATOMS;q++) {
          if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
            // start index for i-th El in the second mol
            // compute the distance
            //printf ("m l %li %li\n", m, l);
            double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
            double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
            double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
            // in cost matrix we cannot operate with doubles
            long long int d = sqrt(x*x + y*y + z*z)*1e6;
            COST[m][l] = d;
            l=l+1;
          }
        }
      l=0; 
      m=m+1;
      }
    }
    /*
    printf ("Cost matrix (%lli %lli):\n", IND, IND);
    for (p=0;p<IND;p++) {
      for (q=0;q<IND;q++) {
        printf ("%lli ", COST[p][q]);
      }
    printf ("\n");
    }
    */
    // array storing matched elements in mol2
    //long int MATCH[IND];
    long int *MATCH = malloc(IND*sizeof(long int));  
    //printf ("HUNGMAIN START\n");
    HUNGMAIN(IND, IND, COST, &MATCH);
    //printf ("HUNGMAIN FINISH\n");
    /*
    for (p=0;p<IND;p++) {
      printf("%li ", MATCH[p]);
    }
    printf("\n");
    */
    m=0;
    l=0;
    for (p=0;p<NATOMS;p++) {
      if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
        // start index for i-th El
        // iterate via another molecule
        
        for (q=0;q<NATOMS;q++) {
          if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
            if (MATCH[m] == l) {
              // matched
              (*coord1tmp)[p].N = q;
            }
            l=l+1;
          }
        }
      // assigment
      l=0; 
      m=m+1;
      }
    }
  free (MATCH);
  free (COST);
  }

  /*
  for (i=0;i<NATOMS;i++) {
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
    }
    else {
      printf ("N is negative. Stop.\n");
      exit(0);
    }
  }
  */ 
  // populate coord (temporary array) according to renumbering
  for (i=0;i<NATOMS;i++) {
    long int it = (*coord1tmp)[i].N;
    strcpy(coord[it].El,(*coord1tmp)[i].El);
    for (j=0;j<3;j++) {
      coord[it].xyz[j] = (*coord1tmp)[i].xyz[j];
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

void RENUMBHUNGGRAPH (Atom **coord1tmp, Atom *coord2, long int NATOMS, char (*UnEl)[3], long int UE, bool IDN[][NATOMS]) {

  long int i=0;
  long int j=0;
  // temporary array for renumeration
  Atom *coord = malloc(NATOMS*sizeof(Atom));
  // another tmp array to store best numeration
  Atom *coord1best = malloc(NATOMS*sizeof(Atom));
  // array with the number of each element
  long int UnElN[UE];
  for (i=0;i<UE;i++) {
    UnElN[i] = 0;
  }
  for (i=0;i<NATOMS;i++) {
    (*coord1tmp)[i].N = -1; // will store its correspondence to sequence in
    for (j=0;j<UE;j++) {
      if (strcmp(coord2[i].El,UnEl[j]) == 0) {
        UnElN[j] = UnElN[j]+1;
      }
    }                        // the second structure
    // second molecule
    coord2[i].N = -1; // will be changed to positive (itself) if already matched 
  }
//  for (i=0;i<UE;i++) {
//    printf ("%li ", UnElN[i]);
//  }
//  printf("\n");
  // for each element we construc the cost matrix
  for (i=0;i<UE;i++) {
    long long int IND = UnElN[i];
    //long int COST[IND][IND];
    long long int (*COST)[IND] = malloc(IND*sizeof(*COST));
    long int p=0;
    long int q=0;
    long int m=0;
    long int l=0;
    //printf ("Start %li %li\n", i, IND);
    // populate matrix with interatomic distances
    for (p=0;p<NATOMS;p++) {
      if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
        // start index for i-th El
        // iterate via another molecule
        for (q=0;q<NATOMS;q++) {
          if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
            // start index for i-th El in the second mol
            // compute the distance
            //printf ("m l %li %li\n", m, l);
            double x = (*coord1tmp)[p].xyz[0] - coord2[q].xyz[0];
            double y = (*coord1tmp)[p].xyz[1] - coord2[q].xyz[1];
            double z = (*coord1tmp)[p].xyz[2] - coord2[q].xyz[2];
            // in cost matrix we cannot operate with doubles
            long long int d = sqrt(x*x + y*y + z*z)*1e6;
            // check if bonding pattern is similar
            if (IDN[p][q] == true) {
              d = d/1000 ;
              //printf ("%li %li - similar\n", p+1, q+1);
            }
            COST[m][l] = d;
            l=l+1;
          }
        }
      l=0; 
      m=m+1;
      }
    }
    /*
    printf ("Cost matrix (%lli %lli):\n", IND, IND);
    for (p=0;p<IND;p++) {
      for (q=0;q<IND;q++) {
        printf ("%lli ", COST[p][q]);
      }
    printf ("\n");
    }
    */
    // array storing matched elements in mol2
    //long int MATCH[IND];
    long int *MATCH = malloc(IND*sizeof(long int));  
    //printf ("HUNGMAIN START\n");
    HUNGMAIN(IND, IND, COST, &MATCH);
    //printf ("HUNGMAIN FINISH\n");
    /*
    for (p=0;p<IND;p++) {
      printf("%li ", MATCH[p]);
    }
    printf("\n");
    */
    m=0;
    l=0;
    for (p=0;p<NATOMS;p++) {
      if (   strcmp((*coord1tmp)[p].El, UnEl[i]) == 0 ) {
        // start index for i-th El
        // iterate via another molecule
        
        for (q=0;q<NATOMS;q++) {
          if (  strcmp( coord2[q].El, UnEl[i]) == 0  ) {
            if (MATCH[m] == l) {
              // matched
              (*coord1tmp)[p].N = q;
            }
            l=l+1;
          }
        }
      // assigment
      l=0; 
      m=m+1;
      }
    }
  free (MATCH);
  free (COST);
  }

  /*
  for (i=0;i<NATOMS;i++) {
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
    }
    else {
      printf ("N is negative. Stop.\n");
      exit(0);
    }
  }
  */ 
  // populate coord (temporary array) according to renumbering
  for (i=0;i<NATOMS;i++) {
    long int it = (*coord1tmp)[i].N;
    strcpy(coord[it].El,(*coord1tmp)[i].El);
    for (j=0;j<3;j++) {
      coord[it].xyz[j] = (*coord1tmp)[i].xyz[j];
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

