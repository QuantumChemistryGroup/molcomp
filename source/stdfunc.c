#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "maintypes.h"
#include "atominfo.h"


int comp(const void* a,const void* b) {
  	return *(long int*)a - *(long int*)b;
}


bool BONDEQUAL (long int p, long int q, Atom *coord1, Atom *coord2, long int NATOMS, bool MOL2) {

bool EQ = true;

if (MOL2 == false) {
  if (coord1[p].nn != coord1[q].nn) {
    EQ = false;
    return EQ;
  }
}
else {
  if (coord1[p].nn != coord2[q].nn) {
    EQ = false;
    return EQ;
  }
}

long int N=coord1[p].nn;
long int array1[N];
long int array2[N];
// populate it with element symbols (numbers)
long int i=0;
for (i=0;i<N;i++) {
  long int seq1 = coord1[p].n[i]-1;
  array1[i] = ATOMINFO(coord1[seq1].El).Z;
  //printf ("%li ", array1[i]); 
  long int seq2=-1;
  if (MOL2 == false) {
    seq2 = coord1[q].n[i]-1;
    array2[i] = ATOMINFO(coord1[seq2].El).Z;
  }
  else {
    seq2 = coord2[q].n[i]-1;
    array2[i] = ATOMINFO(coord2[seq2].El).Z;
  }
  //printf ("%li ", array2[i]); 
}

qsort(array1, N, sizeof(long int), comp);
qsort(array2, N, sizeof(long int), comp);

for (i=0;i<N;i++) {
  if (array1[i] != array2[i]) {
    EQ = false;
    break;
  }
}

// sort both array1 and array2

 

return EQ;

}

void POPADJ (long int NATOMS, long int ADJ[][NATOMS], Atom *coord) {

long int i=0;
long int j=0;

for (i=0;i<NATOMS;i++) {
  for (j=0;j<NATOMS;j++) {
    long int l=0;
    for (l=0;l<coord[j].nn;l++) {
      if (coord[j].n[l]-1 == i) {
        ADJ[i][j] = 1;
      }
    }
  }
}

}

bool RADIAL (long int p, long int q, Atom *coord1, Atom *coord2, long int NATOMS, bool MOL2, long int ADJ1[][NATOMS], long int ADJ2[][NATOMS]) {

bool EQ = true;

/* FIRST QUICK CHECK - immediate surrounding coord number*/

if (MOL2 == false) {
  if (coord1[p].nn != coord1[q].nn) {
    EQ = false;
    return EQ;
  }
}
else {
  if (coord1[p].nn != coord2[q].nn) {
    EQ = false;
    return EQ;
  }
}

/* SECOND QUICK CHECK - immediate surrounding */

long int N=coord1[p].nn;
long int array1[N];
long int array2[N];
// populate it with element symbols (numbers)
long int i=0;
for (i=0;i<N;i++) {
  long int seq1 = coord1[p].n[i]-1;
  array1[i] = ATOMINFO(coord1[seq1].El).Z;
  //printf ("%li ", array1[i]); 
  long int seq2=-1;
  if (MOL2 == false) {
    seq2 = coord1[q].n[i]-1;
    array2[i] = ATOMINFO(coord1[seq2].El).Z;
  }
  else {
    seq2 = coord2[q].n[i]-1;
    array2[i] = ATOMINFO(coord2[seq2].El).Z;
  }
  //printf ("%li ", array2[i]); 
}

qsort(array1, N, sizeof(long int), comp);
qsort(array2, N, sizeof(long int), comp);

for (i=0;i<N;i++) {
  if (array1[i] != array2[i]) {
    EQ = false;
    break;
  }
}
long int j=0;
/* BFS search */
if (EQ != false) {
// 1 is for p, 2 is for q
long int n;

long int visited1[NATOMS];
long int visited2[NATOMS];

long int ind = 0;
for (ind=0;ind<NATOMS;ind++) { 
  visited1[ind] = 0;
  visited2[ind] = 0; 
}

long int queue1[NATOMS];
long int front1 = 0;
long int rear1 = 0;
long int start1, i;

long int queue2[NATOMS];
long int front2 = 0;
long int rear2 = 0;
long int start2;
  
n=NATOMS;

start1=p;
visited1[start1] = 1;
queue1[rear1] = start1;
rear1=rear1+1;

start2=q;
visited2[start2] = 1;
queue2[rear2] = start2;
rear2=rear2+1;

long int level1[NATOMS];
long int maxl1 = 0;
long int level2[NATOMS];
long int maxl2 = 0;

for (i=0;i<NATOMS;i++) {
  level1[i] = 0;
  level2[i] = 0;
}
/* for p */
long int pl1 = 0;
long int pl2 = 0;
long int pl1arr[NATOMS];
long int plind = 0;
long int pl2arr[NATOMS];
for (i=0;i<NATOMS;i++) {
  pl1arr[i] = 0;
  pl2arr[i] = 0;
}

while ( (front1 < rear1) || (front2 < rear2) ) {
  if ( ( (front1<rear1) != true)  || ((front2<rear2) != true)  ) {
    EQ = false;
    return EQ;
  }
  long int u1 = queue1[front1];
  front1 = front1 + 1;
  long int u2 = queue2[front2];
  front2 = front2 + 1;
  // check inconsistency
  if ( (pl1!=pl2) || (level1[u1] != level2[u2] ) ) {
    EQ = false;
    return EQ;
  }
  if (pl1 != level1[u1]) {
    // start comparison
    // modify seq numbers to Z
    for (j=0;j<plind;j++) {
      // current seq numbers
      long int k = pl1arr[j];
      long int m = pl2arr[j];
      // replace to Z
      pl1arr[j] = ATOMINFO(coord1[k].El).Z;
      if (MOL2 == true) {
        pl2arr[j] = ATOMINFO(coord2[m].El).Z;
      }
      else {
        pl2arr[j] = ATOMINFO(coord1[m].El).Z;
      }
    }
    // sorting by Z
    qsort(pl1arr, plind, sizeof(long int), comp);
    qsort(pl2arr, plind, sizeof(long int), comp);
    // final check 
    for (j=0;j<plind;j++) {
      if (pl1arr[j] != pl2arr[j]) {
        EQ = false;
        return EQ;
      }
    }
    // clean
    for (j=0;j<plind;j++) {
      pl1arr[j] = 0;
      pl2arr[j] = 0;
    }
    plind=0;
    // start populating - new level
    pl1arr[plind] = u1;
    pl2arr[plind] = u2;
    plind = plind+1;
  }
  else {
    pl1arr[plind] = u1;
    pl2arr[plind] = u2;
    plind = plind+1;
  }
  // p
  for (i=0;i<n;i++) {
    if (ADJ1[u1][i] == 1 && visited1[i] == 0) {
      visited1[i] = 1;
      level1[i] = level1[u1]+1;
      if (level1[i] > maxl1) {
        maxl1 = level1[i];
      }
      queue1[rear1] = i;
      rear1 = rear1+1;
    }
  }
  pl1 = level1[u1];
  // q
  if (MOL2 == false) {
    for (i=0;i<n;i++) {
      if (ADJ1[u2][i] == 1 && visited2[i] == 0) {
        visited2[i] = 1;
        level2[i] = level2[u2]+1;
        if (level2[i] > maxl2) {
          maxl2 = level2[i];
        }
        queue2[rear2] = i;
        rear2 = rear2+1;
      }
    }
  }
  else {
    for (i=0;i<n;i++) {
      if (ADJ2[u2][i] == 1 && visited2[i] == 0) {
        visited2[i] = 1;
        level2[i] = level2[u2]+1;
        if (level2[i] > maxl2) {
          maxl2 = level2[i];
        }
        queue2[rear2] = i;
        rear2 = rear2+1;
      }
    }
  }
  pl2 = level2[u2];
}

//printf ("Survivor: %li %li Level: %li previous: %li\n", p+1, q+1, level1[p], pl1);

return EQ;
}

return EQ;

}
