#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "maintypes.h"
#include <string.h>
#include <ctype.h>
#include <omp.h>
#include <time.h>
#include "fread.h"
#include "inertia.h"
#include "kabsh.h"
#include "kabshdist.h"
#include "kabshangle.h"
#include "renumb.h"
#include "renumb-hung.h"
#include "getconnect.h"
#include "stdfunc.h"
#include "txyzread.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>


char *strlwr(char *str);
// program to compute RMSD between two sets of coordinates
bool VER = false;
bool PRINT = false; // print the first modifed structure to file (best fit) 
char ALG[100];
char DIST[100]; // print average dist & max dist or angle?
bool DISTRATIO = false; // consider also the ratio between the largest and average distance
double DISTLR = 1.0; // largest distance above which conformers are different 
double DISTSM = 0.75; // largest distance for which ratio is cheked
double DISTRT = 4.0; // ratio between largest distance and average dist
bool STEREO = false;
bool CONNECT = true; // use connectivity information
bool BINDH2HEAVY = true; // keep hydrogens bonded to heavy atoms upon renumbering
bool MOL2 = false; // use also the connectivity of the second molecule
bool PLPOPIDN = false; // parallel populate IDN 
bool WRITETOPO = false; // write 
bool READTOPO = false; // read reference topology
double LINEARTHLD = 5.0; // threshold for linearity. Torsions with linear atoms are excluded.
int main(int argc, char *argv[])
{

// check for arguments

// coordinates of the first atom
if (argv[1] != NULL && argv[2] != NULL) {
}
else {
  printf ("Files with coordinates are not specified.\n");
  exit(0);
}

// allocate the initial molecule.
long int NAT1 = 0;
long int NAT2 = 0;
long int NATDEF = 0;

Atom *coord1 = malloc(1*sizeof(Atom));
Atom *coord2 = malloc(1*sizeof(Atom));
// temporary coord to get the default connect
// and compare it with the modified
// to ensure the right connectivity is symmetrically kept
Atom *coorddef = malloc(1*sizeof(Atom));


// * Files to read and write  
FILE *fp;

// iterators
long int i=0;
long int j=0;

// get the number of arguments apart from files with xyz
for (i=3;i<argc;i++) {
  // populate VER
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--verbose=", 10) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      VER = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      VER = false;
    }
    else {
      printf ("Your verbosity --verbose=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

// open the file to write
char fileout[300];
snprintf(fileout, sizeof(fileout), "%.*s-mod.xyz", (int)strlen(argv[1])-4, argv[1] );

for (i=3;i<argc;i++) {
  // populate VER
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--print=", 8) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      PRINT = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      PRINT = false;
    }
    else {
      printf ("Your printing --print=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

for (i=3;i<argc;i++) {
  // populate STEREO
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--keepenantio=", 14) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      STEREO = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      STEREO = false;
    }
    else {
      printf ("Your --keepenantio=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

for (i=3;i<argc;i++) {
  // populate CONNECT
  // defauls is true
  if (  strncmp(strlwr(argv[i]), "--connect=", 10) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      CONNECT = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      CONNECT = false;
    }
    else {
      printf ("Your --connect=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

for (i=3;i<argc;i++) {
  // populate MOL2
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--writetopo=", 12) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      WRITETOPO = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      WRITETOPO = false;
    }
    else {
      printf ("Your --writetopo=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

for (i=3;i<argc;i++) {
  // populate MOL2
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--readtopo=", 11) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      READTOPO = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      READTOPO = false;
    }
    else {
      printf ("Your --readtopo=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

// linearity threshold

for (i=3;i<argc;i++) {
  // populate LINEARTHLD
  // default is 5 degrees
  if (  strncmp(strlwr(argv[i]), "--linearthld=", 13) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    LINEARTHLD = atof (ptr);
  break;
  }
}


for (i=3;i<argc;i++) {
  // populate MOL2
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--mol2=", 7) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      MOL2 = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      MOL2 = false;
    }
    else {
      printf ("Your --mol2=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}


for (i=3;i<argc;i++) {
  // populate MOL2
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--bindh2heavy=", 14) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      BINDH2HEAVY = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      BINDH2HEAVY = false;
    }
    else {
      printf ("Your --bindh2heavy=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

for (i=3;i<argc;i++) {
  // populate MOL2
  // defauls is false
  if (  strncmp(strlwr(argv[i]), "--parallelidn=", 14) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "true") == 0 ) {
      PLPOPIDN = true;
    }
    else if (strcmp (ptr, "false") == 0) {
      PLPOPIDN = false;
    }
    else {
      printf ("Your --parallelidn=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}


// get ALG
strcpy(ALG, "std"); // default
for (i=3;i<argc;i++) {
  // populate ALG
  // default is --std
  if (  strncmp(strlwr(argv[i]), "--alg=", 6) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "std") == 0 ) {
      strcpy(ALG, "std");
    }
    else if (strcmp (ptr, "adv") == 0) {
      strcpy(ALG, "adv");
    }
    else if (strcmp (ptr, "kabsch") == 0) {
      strcpy(ALG, "kabsch");
    }
    else {
      printf ("Your --alg=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

// get DIST
strcpy(DIST, "false"); // default
for (i=3;i<argc;i++) {
  // populate DIST
  // default is false
  if (  strncmp(strlwr(argv[i]), "--compare=", 10) == 0) {
    char *ptr;
    ptr = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr = strtok(NULL, "=");
    if ( strcmp (ptr, "false") == 0 ) {
      strcpy(DIST, "false");
    }
    else if (strcmp (ptr, "dist") == 0) {
      strcpy(DIST, "dist");
    }
    else if (strcmp (ptr, "angle") == 0) {
      strcpy(DIST, "angle");
      CONNECT=true;
    }
    else {
      printf ("Your --compare=%s is not clear\n", ptr);
      exit(0);
    }
  break;
  }
}

if (strcmp(DIST, "dist") == 0) {
  // then check if DISTRATIO is true
  for (i=3;i<argc;i++) {
    // populate DISTRATIO
    // defauls is false
    if (  strncmp(strlwr(argv[i]), "--distratio=", 12) == 0) {
      char *ptr;
      ptr = strtok(strlwr(argv[i]), "=");
      // take second argument
      ptr = strtok(NULL, "=");
      if ( strcmp (ptr, "true") == 0 ) {
        DISTRATIO = true;
      }
      else if (strcmp (ptr, "false") == 0) {
        DISTRATIO = false;
      }
      else {
        printf ("Your --distratio=%s is not clear\n", ptr);
        exit(0);
      }
    break;
    }
  }
}

if (DISTRATIO == true) {
  // then we need to populate 3 releated criteria
  // DISTLM 
  for (i=3;i<argc;i++) {
  if (  strncmp(strlwr(argv[i]), "--distlr=", 9) == 0) {
    char *ptr1;
    ptr1 = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr1 = strtok(NULL, "=");
    DISTLR = atof (ptr1);
  }
  // DISTSM
  else if (strncmp(strlwr(argv[i]), "--distsm=", 9) == 0) {
    char *ptr1;
    ptr1 = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr1 = strtok(NULL, "=");
    DISTSM = atof (ptr1);
  }
  // DISTRT
  else if (strncmp(strlwr(argv[i]), "--distrt=", 9) == 0) {
    char *ptr1;
    ptr1 = strtok(strlwr(argv[i]), "=");
    // take second argument
    ptr1 = strtok(NULL, "=");
    DISTRT = atof (ptr1);
  }
  }
}

// coordinates of the first atom
if (argv[1] != NULL) {
  fp = fopen (argv[1], "r");
  if (fp != NULL) {
    if (VER) printf ("Start reading first input file: %s \n", argv[1]);
    // just to compare initial & modified connectivity
    if ((READTOPO==true) && (MOL2 == false) ) {
      FREAD(fp, &coorddef, &NATDEF);
    } 
    FREAD(fp, &coord1, &NAT1);
    fclose(fp);
  }
  else {
    printf ("File 1 with coordinates does not exist.\n");
  }
}
else {
  printf ("File 1 with coordinates is not specified.\n");
  exit(0);
}

// coordinates of the first atom
if (argv[2] != NULL) {
  fp = fopen (argv[2], "r");
  if (fp != NULL) {
    if (VER) printf ("Start reading second input file: %s \n", argv[2]);
    // just to compare initial & modified connectivity
    if ((READTOPO==true) && (MOL2 == true)) {
      FREAD(fp, &coorddef, &NATDEF);
    } 
    FREAD(fp, &coord2, &NAT2);
    fclose(fp);
  }
  else {
    printf ("File 2 with coordinates does not exist.\n");
  }
}
else {
  printf ("File 2 with coordinates is not specified.\n");
  exit(0);
}



// check if NATOMS are equal
long int NATOMS=0;
if (NAT1 != NAT2) {
  printf ("NAT1 is not equal to NAT2. Stop. %li %li \n", NAT1, NAT2);
  exit(0);
}
else {
  NATOMS = NAT1;
}

long int (*ADJ1)[NATOMS] = malloc(NATOMS*sizeof(*ADJ1));
long int (*ADJ2)[NATOMS] = malloc(NATOMS*sizeof(*ADJ2));

if (CONNECT == true) {
  // get connectivity for the first molecule
  if (READTOPO==false) {
    GETCONNECT (&(coord1), NATOMS);
  }
  else {
    // for comparison get the defalt connect
    GETCONNECT (&(coorddef), NATOMS);
    fp = fopen ("topo.txt", "r");
    if (fp == NULL) { 
      printf ("FILE WITH TOPOLOGY %s DOES NOT EXIST. STOP.\n", "topo.txt");
      exit(0);
    }
    TXYZREAD (fp, &coord1, NATOMS);
    fclose(fp);
/*
It is dangerous to read topo.txt by multiple files: their reference topology is different. Hence, if we add some new bond to the first Li atom in topo.txt, and do not add connection to respective atoms, and if default topo has the same Li bond pattern, it will further discover in topo.txt missing connections to 1, and by checking it with default topo, will go to 1 and remove the bonds. 
*/    
    //SYMMETRYCHECKDEF (coorddef, &(coord1), NATOMS);
  }
  SYMMETRYCHECK (&(coord1), NATOMS);
  // compute the connectivity matrix
  POPADJ (NATOMS, ADJ1, coord1);
  if (WRITETOPO == true) {
    // writing topology to file
    fp = fopen ("topo.txt", "w");
    fprintf (fp, "%3li MOLECULAR CONNECTIVITY \n", NATOMS);
    for (i=0;i<NATOMS;i++) {
      fprintf(fp, "%5li %5s %10.6f %10.6f %10.6f %8s", i+1, coord1[i].El, coord1[i].xyz[0], coord1[i].xyz[1], coord1[i].xyz[2], "UDF");
      for (j=0;j<coord1[i].nn;j++) {
        fprintf(fp, "%5li ", coord1[i].n[j]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  /*for (i=0;i<NATOMS;i++) {
    for (j=0;j<NATOMS;j++) {
      if (ADJ1[i][j]!=0) {
        printf ("  %li %li %li", i+1, j+1, ADJ1[i][j]);
      }
    }
  printf ("\n");
  }*/
  if (MOL2 == true) {
    if (READTOPO == false) {
      GETCONNECT (&(coord2), NATOMS);
    }
    else {
      fp = fopen ("topo.txt", "r");
      if (fp == NULL) { 
        printf ("FILE WITH TOPOLOGY %s DOES NOT EXIST. STOP.\n", "topo.txt");
        exit(0);
      }
      TXYZREAD (fp, &coord2, NATOMS);
      fclose(fp);
/*
It is dangerous to read topo.txt by multiple files: their reference topology is different. Hence, if we add some new bond to the first Li atom in topo.txt, and do not add connection to respective atoms, and if default topo has the same Li bond pattern, it will further discover in topo.txt missing connections to 1, and by checking it with default topo, will go to 1 and remove the bonds. 
*/    
      //SYMMETRYCHECKDEF (coorddef, &(coord2), NATOMS);
    }
    SYMMETRYCHECK (&(coord2), NATOMS);
    POPADJ (NATOMS, ADJ2, coord2);
    if (WRITETOPO == true) {
      // writing topology to file
      printf ("As --mol2=true connectivity of second molecule is written.\n");
      fp = fopen ("topo.txt", "w");
      fprintf (fp, "%3li MOLECULAR CONNECTIVITY \n", NATOMS);
      for (i=0;i<NATOMS;i++) {
        fprintf(fp, "%5li %5s %10.6f %10.6f %10.6f %8s", i+1, coord2[i].El, coord2[i].xyz[0], coord2[i].xyz[1], coord2[i].xyz[2], "UDF");
        for (j=0;j<coord2[i].nn;j++) {
          fprintf(fp, "%5li ", coord2[i].n[j]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }
}
// matrix of identity of atoms p and q
bool (*IDN)[NATOMS] = malloc (NATOMS*sizeof(*IDN));
if (CONNECT == true && PLPOPIDN == true) {
  #pragma omp parallel private (j)
  {
  // fill the matrix
   //printf ("Start IDN.\n");
  #pragma omp for
  for (i=0;i<NATOMS;i++) {
    //printf ("Layer: %li\n", i);
    for (j=0;j<NATOMS;j++) {
      bool EqEl = false;
      if (MOL2 == true) {
        if (strcmp(coord1[i].El, coord2[j].El) == 0) EqEl = true;
      }
      else {
        if (strcmp(coord1[i].El, coord1[j].El) == 0) EqEl = true;
      }
      if (EqEl == false) {
        IDN[i][j] = false;
      }
      else {
        IDN[i][j] = RADIAL(i,j,(coord1),coord2,NATOMS,MOL2,ADJ1,ADJ2);
      }    
    }
  }
  } // close pragma
   //printf ("Finish IDN.\n");
}
else if (CONNECT == true && PLPOPIDN == false) {
  // fill the matrix
   //printf ("Start IDN.\n");
  for (i=0;i<NATOMS;i++) {
    for (j=0;j<NATOMS;j++) {
      bool EqEl = false;
      if (MOL2 == true) {
        if (strcmp(coord1[i].El, coord2[j].El) == 0) EqEl = true;
      }
      else {
        if (strcmp(coord1[i].El, coord1[j].El) == 0) EqEl = true;
      }
      if (EqEl == false) {
        IDN[i][j] = false;
      }
      else {
        IDN[i][j] = RADIAL(i,j,(coord1),coord2,NATOMS,MOL2,ADJ1,ADJ2);
      }    
    }
  }
}


// first molecule initial as KABSH changes it
// to store renumbered structure
Atom *coord1mod = malloc(NAT1*sizeof(Atom));
// to store simple Kabsch (no renumbering)
Atom *coord1modk = malloc(NAT1*sizeof(Atom));
// best coord1 after renumbering
Atom *coord1best = malloc(NAT1*sizeof(Atom));

// check that everything is stored in coord1 and coord2

// principal moments of inertia of the first mol
double IP1[3] = {0,0,0};
double IP2[3] = {0,0,0};

// its eigenfunctions
double EF1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double EF2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

// place molecule in the center of mass
// get principle moments of inertia and eigenvetors



INERTIA (&coord1, IP1, NAT1, EF1);
INERTIA (&coord2, IP2, NAT2, EF2);

// array of unique elements
char (*UnEl)[3] = malloc(1*sizeof(*UnEl));
char (*UnEl0)[3] = malloc(1*sizeof(*UnEl));
// iterator of unique elements
long int UE = 0;
long int UE0 = 0;
// number of elements of each kind
long int *ElNum = malloc(1*sizeof(long int));
long int *ElNum0 = malloc(1*sizeof(long int));

for (i=0;i<NAT1;i++) {
  // first molecule
  if (i == 0) {
    strcpy(UnEl[0],coord2[i].El);
    strcpy(UnEl0[0],coord1[i].El);
    UE = UE + 1;
    UE0 = UE0 + 1;
    UnEl = realloc (UnEl, (UE+1)*sizeof(*UnEl));
    UnEl0 = realloc (UnEl0, (UE0+1)*sizeof(*UnEl0));
    ElNum = realloc (ElNum, (UE+1)*sizeof(long int));
    ElNum0 = realloc (ElNum0, (UE0+1)*sizeof(long int));
    ElNum[UE-1] = 1;
    ElNum0[UE0-1] = 1;
  }
  else {
    // check if we have this element
    bool MATCH = false;
    for (j=0;j<UE;j++) {
      if (strcmp(coord2[i].El, UnEl[j]) == 0) {
        MATCH = true;
        ElNum[j] = ElNum[j] + 1;
        //printf ("%s %li\n", coord2[i].El, ElNum[j]);
        break;
      }
    }
    if (MATCH == false) {
      strcpy(UnEl[UE],coord2[i].El);
      UE=UE+1;
      UnEl = realloc (UnEl, (UE+1)*sizeof(*UnEl));
      ElNum = realloc (ElNum, (UE+1)*sizeof(long int));
      ElNum[UE-1] = 1;
    }
    // another molecule
    MATCH = false;
    for (j=0;j<UE0;j++) {
      if (strcmp(coord1[i].El, UnEl0[j]) == 0) {
        MATCH = true;
        ElNum0[j] = ElNum0[j] + 1;
        break;
      }
    }
    if (MATCH == false) {
      strcpy(UnEl0[UE0],coord1[i].El);
      UE0=UE0+1;
      UnEl0 = realloc (UnEl0, (UE0+1)*sizeof(*UnEl0));
      ElNum0 = realloc (ElNum0, (UE0+1)*sizeof(long int));
      ElNum0[UE0-1] = 1;
    }  
  }
}

// check if the number of unique elements is the same 
// in both structures

if (UE != UE0) {
  printf ("Number of unique elements is not the same in both structures. Exit.\n");
  exit(0);
} 
/*
for (i=0;i<UE;i++) {
  printf ("%s\n", UnEl[i]);
  printf ("%li\n", ElNum[i]);
}
for (i=0;i<UE0;i++) {
  printf ("%s\n", UnEl0[i]);
  printf ("%li\n", ElNum0[i]);
}
*/
for (i=0;i<UE;i++) {
  for (j=0;j<UE0;j++) {
    if ( strcmp(UnEl[i],UnEl0[j]) == 0) {
      // check the number of unique elements
      if (ElNum[i] != ElNum0[j]) {
        // mismatch
        printf("Mismatch in %s mol1: %li mol2: %li. Close.\n", UnEl[i], ElNum0[j], ElNum[i]);
        exit(0);
      }
    }
  }
}






// ORIENT ALONG PRINCIPAL AXES

for (i=0;i<NAT1;i++){
  double xyz[3] = {coord1[i].xyz[0], coord1[i].xyz[1], coord1[i].xyz[2]};
  coord1[i].xyz[0] = xyz[0]*EF1[0][0] + xyz[1]*EF1[1][0] + xyz[2]*EF1[2][0];
  coord1[i].xyz[1] = xyz[0]*EF1[0][1] + xyz[1]*EF1[1][1] + xyz[2]*EF1[2][1];
  coord1[i].xyz[2] = xyz[0]*EF1[0][2] + xyz[1]*EF1[1][2] + xyz[2]*EF1[2][2];
}


for (i=0;i<NAT2;i++){
  double xyz[3] = {coord2[i].xyz[0], coord2[i].xyz[1], coord2[i].xyz[2]};
  coord2[i].xyz[0] = xyz[0]*EF2[0][0] + xyz[1]*EF2[1][0] + xyz[2]*EF2[2][0];
  coord2[i].xyz[1] = xyz[0]*EF2[0][1] + xyz[1]*EF2[1][1] + xyz[2]*EF2[2][1];
  coord2[i].xyz[2] = xyz[0]*EF2[0][2] + xyz[1]*EF2[1][2] + xyz[2]*EF2[2][2];
}


/*
printf ("First mol\n");
for (i=0;i<NAT1;i++){
  printf("%s %10.6f %10.6f %10.6f\n", coord1[i].El, coord1[i].xyz[0], coord1[i].xyz[1], coord1[i].xyz[2]);
}


printf ("Second mol\n");
for (i=0;i<NAT2;i++){
  printf("%s %10.6f %10.6f %10.6f\n", coord2[i].El, coord2[i].xyz[0], coord2[i].xyz[1], coord2[i].xyz[2]);
}
*/

/* FINALLY OUR TWO MOLECULES ARE ORIENTED ALONG THEIR PRINCIPAL ROTATION AXES*/

// try to get all mirror images
/*
  #  0 +x +y +z    
  #  1 +x -y +z    
  #  2 -x -y +z    
  #  3 -x +y +z    
  #  4 -x +y -z    
  #  5 +x +y -z    
  #  6 +x -y -z    
  #  7 -x -y -z    
*/

double Imag[8][3] = {
{1,1,1},
{1,-1,1},
{-1,-1,1},
{-1,1,1},
{-1,1,-1},
{1,1,-1},
{1,-1,-1},
{-1,-1,-1},
};

long int m=0;
// check only for the first molecule
double RMSDmin = 1e10;

if (strcmp(ALG, "std")==0) {
  for (i=0;i<8;i++) {
    // introduce modified first molecule
    //printf ("Start %li \n", i);
    if (STEREO == true) {
      // check if the overall multiplication is not -1
      long int REFL = Imag[i][0]*Imag[i][1]*Imag[i][2];
      if (REFL == -1) {
        continue;
      }
    }
    //printf ("Continue %li \n", i);
    Atom *coord1tmp = malloc(NAT1*sizeof(Atom));
    for (j=0;j<NAT1;j++) {
      strcpy (coord1tmp[j].El, coord1[j].El);
      for (m=0;m<3;m++) {
        coord1tmp[j].xyz[m] = coord1[j].xyz[m]*Imag[i][m];
      }
      if (CONNECT==true) {
        coord1tmp[j].nn = coord1[j].nn;
        long int ind = coord1tmp[j].nn;
        for (m=0;m<ind;m++) {
          coord1tmp[j].n[m] = coord1[j].n[m];
        }
      }
    }
    // modified first molecule
    // try to get the best match in numbering of the first molecule
    double RMSD = 0;
    //printf ("Before renumb\n");
    //KABSH (coord1tmp, coord1, &coord1mod, NAT1, &RMSD);
  //  RENUMB (&coord1tmp, coord2, NAT1);
    if (CONNECT != true) {
      RENUMBHUNG (&coord1tmp, coord2, NAT1, UnEl, UE);
    }
    else {
      // OLD FUNCTION, USES ONLY IMMEDIATE SURROUNDING INSTEAD OF GRAPH
      if (BINDH2HEAVY == true) { 
        RENUMBHUNGCONNECT (&coord1tmp, coord2, NAT1, UnEl, UE, IDN); // using graph + keep hyndrogens bonded
      }
      else {
        RENUMBHUNGGRAPH (&coord1tmp, coord2, NAT1, UnEl, UE, IDN);
      }
    }
    // run standard KABSH
    //printf ("After renumb\n");
    KABSH (coord1tmp, coord2, &coord1mod, NAT1, &RMSD);
    //printf ("RMSD = %f\n", RMSD);
    if (RMSD < RMSDmin) {
      RMSDmin = RMSD;
      for (j=0;j<NAT1;j++) {
        strcpy(coord1best[j].El, coord1tmp[j].El);
        for (m=0;m<3;m++) {
          coord1best[j].xyz[m] = coord1tmp[j].xyz[m];
        }
        if (CONNECT==true) {
          coord1best[j].nn = coord1tmp[j].nn;
          long int ind = coord1best[j].nn;
          for (m=0;m<ind;m++) {
            coord1best[j].n[m] = coord1tmp[j].n[m];
          }
        }
        coord1best[j].N = coord1tmp[j].N;
      }
    } 
    free(coord1tmp);
  }
  if (VER) printf ("Standard Kabsch, permutation.\n");
  double RMSD = 0;
  double RMSDST = 0;
  if (strcmp(DIST, "false") == 0) {
    KABSH (coord1best, coord2, &coord1mod, NAT1, &RMSD);
    if (VER) {
      printf ("Kabsch (renumbering) %f\n", RMSD);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    KABSH (coord1, coord2, &coord1modk, NAT1, &RMSDST);
    if (VER) {
      printf ("Kabsch (no renumbering) %f\n", RMSDST);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1modk[i].El, coord1modk[i].xyz[0], coord1modk[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    if (RMSD<RMSDST) {
      printf ("%f \n", RMSD);
      if (PRINT) {
        FPRINT (fileout, coord1mod, NATOMS, RMSD);
      }
    }
    else {
      printf ("%f \n", RMSDST);
      if (PRINT) {
        FPRINT (fileout, coord1modk, NATOMS, RMSDST);
      }
    }
  }   
  else if (strcmp(DIST, "dist") == 0) {
    double DAV=0;
    double DMAX=0;
    double DAVST=0;
    double DMAXST=0;
    KABSHDIST (coord1best, coord2, &coord1mod, NAT1, &RMSD, &DAV, &DMAX);
    if (VER) {
      printf ("Kabsch (renumbering) %f %f %f\n", RMSD, DAV, DMAX);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    KABSHDIST (coord1, coord2, &coord1modk, NAT1, &RMSDST, &DAVST, &DMAXST);
    if (VER) {
      printf ("Kabsch (no renumbering) %f %f %f\n", RMSDST, DAVST, DMAXST);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1modk[i].El, coord1modk[i].xyz[0], coord1modk[i].xyz[1], coord1modk[i].xyz[2]);
      }
    }
    if (RMSD<RMSDST) {
      if (PRINT) {
        FPRINT (fileout, coord1mod, NATOMS, RMSD);
      }
      if (DISTRATIO==false) {
        printf ("%f %f %f\n", RMSD, DAV, DMAX);
      }
      else {
        long int UN = 0; // uniqiue
        if (DMAX>DISTLR) {
          UN = 1; 
        }
        else {
          if (DMAX>DISTSM) {
            // check for ratio
            double rat = DMAX/DAV;
            if (rat > DISTRT) {
              UN = 1;
            }
          }
        }
      printf ("%li %f %f\n", UN, DAV, DMAX);
      }
    }
    else {
      if (PRINT) {
        FPRINT (fileout, coord1modk, NATOMS, RMSDST);
      }
      if (DISTRATIO==false) {
        printf ("%f %f %f \n", RMSDST, DAVST, DMAXST);
      }
      else {
        long int UN = 0; // uniqiue
        if (DMAXST>DISTLR) {
          UN = 1; 
        }
        else {
          if (DMAXST>DISTSM) {
            // check for ratio
            double rat = DMAXST/DAV;
            if (rat > DISTRT) {
              UN = 1;
            }
          }
        }
      printf ("%li %f %f\n", UN, DAVST, DMAXST);
      }
    }
  }
  else if (strcmp(DIST, "angle") == 0) {
    double AAV=0;
    double AMAX=0;
    double AAVST=0;
    double AMAXST=0;
    KABSHANGLE (coord1best, coord2, &coord1mod, NAT1, &RMSD, &AAV, &AMAX);
    if (VER) {
      printf ("Kabsch (renumbering) %f %f %f\n", RMSD, AAV, AMAX);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    KABSHANGLE (coord1, coord2, &coord1modk, NAT1, &RMSDST, &AAVST, &AMAXST);
    if (VER) {
      printf ("Kabsch (no renumbering) %f %f %f\n", RMSDST, AAVST, AMAXST);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1modk[i].El, coord1modk[i].xyz[0], coord1modk[i].xyz[1], coord1modk[i].xyz[2]);
      }
    }
    if (RMSD<RMSDST) {
        if (PRINT) {
          FPRINT (fileout, coord1mod, NATOMS, RMSD);
        }
        printf ("%f %f %f\n", RMSD, AAV, AMAX);
      }
    else {
        if (PRINT) {
          FPRINT (fileout, coord1modk, NATOMS, RMSDST);
        }
        printf ("%f %f %f \n", RMSDST, AAVST, AMAXST);
      }
    }
  else {
    printf ("Your --compare=%s is not clear.\n", DIST);
    exit(0);    
  }
}
else if (strcmp(ALG, "adv")==0) {
  // advanced, permuation of all axes in the first molecule
  long int PERM[6][3] = {
    {0, 1, 2},
    {0, 2, 1},
    {1, 0, 2},
    {1, 2, 0},
    {2, 0, 1},
    {2, 1, 0}
  }; 
  long int p=0;

  for (p=0;p<6;p++) {
    for (i=0;i<8;i++) {
      // introduce modified first molecule
      if (STEREO == true) {
        // check if the overall multiplication is not -1
        long int REFL = Imag[i][0]*Imag[i][1]*Imag[i][2];
        if (REFL == -1) {
          continue;
        }
      }
      Atom *coord1tmp = malloc(NAT1*sizeof(Atom));
      for (j=0;j<NAT1;j++) {
        strcpy (coord1tmp[j].El, coord1[j].El);
        for (m=0;m<3;m++) {
          long int q = PERM[p][m];
          coord1tmp[j].xyz[m] = coord1[j].xyz[q]*Imag[i][m];
        }
       if (CONNECT == true) {
         coord1tmp[j].nn = coord1[j].nn;
         long int ind = coord1tmp[j].nn;
         for (m=0;m<ind;m++) {
           coord1tmp[j].n[m] = coord1[j].n[m];
         }
       }
      }
      // modified first molecule
      // try to get the best match in numbering of the first molecule
      double RMSD = 0;
      //printf ("Before renumb\n");
      //KABSH (coord1tmp, coord2, &coord1init, NAT1, RMSD);
    //  RENUMB (&coord1tmp, coord2, NAT1);
    if (CONNECT != true) {
      RENUMBHUNG (&coord1tmp, coord2, NAT1, UnEl, UE);
    }
    else {
      if (BINDH2HEAVY == true) { 
        RENUMBHUNGCONNECT (&coord1tmp, coord2, NAT1, UnEl, UE, IDN); // using graph + keep hyndrogens bonded
      }
      else {
        RENUMBHUNGGRAPH (&coord1tmp, coord2, NAT1, UnEl, UE, IDN);
      }
    }
      // run standard KABSH
      //printf ("After renumb\n");
      KABSH (coord1tmp, coord2, &coord1mod, NAT1, &RMSD);
    //  printf ("%li %li %f \n", p, i, RMSD);
    //  printf ("RMSD = %f\n", RMSD);
      if (RMSD < RMSDmin) {
        RMSDmin = RMSD;
        for (j=0;j<NAT1;j++) {
          strcpy(coord1best[j].El, coord1tmp[j].El);
          for (m=0;m<3;m++) {
            coord1best[j].xyz[m] = coord1tmp[j].xyz[m];
          }
          if (CONNECT==true) {
            coord1best[j].nn = coord1tmp[j].nn;
            long int ind = coord1best[j].nn;
            for (m=0;m<ind;m++) {
              coord1best[j].n[m] = coord1tmp[j].n[m];
            }
          }
          coord1best[j].N = coord1tmp[j].N;
        }
      } 
      free(coord1tmp);
    }
  }
  if (VER) printf ("Standard Kabsch, permutation.\n");
  double RMSD = 0;
  double RMSDST = 0;
  if (strcmp(DIST, "false") == 0) {
    KABSH (coord1best, coord2, &coord1mod, NAT1, &RMSD);
    if (VER) {
      printf ("Kabsch (renumbering) %f\n", RMSD);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    KABSH (coord1, coord2, &coord1modk, NAT1, &RMSDST);
    if (VER) {
      printf ("Kabsch (no renumbering) %f\n", RMSDST);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1modk[i].El, coord1modk[i].xyz[0], coord1modk[i].xyz[1], coord1modk[i].xyz[2]);
      }
    }
    if (RMSD<RMSDST) {
      if (PRINT) {
        FPRINT (fileout, coord1mod, NATOMS, RMSD);
      }
      printf ("%f \n", RMSD);
    }
    else {
      if (PRINT) {
        FPRINT (fileout, coord1modk, NATOMS, RMSDST);
      }
      printf ("%f \n", RMSDST);
    }
  }   
  else if (strcmp(DIST, "dist") == 0) {
    double DAV=0;
    double DMAX=0;
    double DAVST=0;
    double DMAXST=0;
    KABSHDIST (coord1best, coord2, &coord1mod, NAT1, &RMSD, &DAV, &DMAX);
    if (VER) {
      printf ("Kabsch (renumbering) %f %f %f\n", RMSD, DAV, DMAX);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    KABSHDIST (coord1, coord2, &coord1modk, NAT1, &RMSDST, &DAVST, &DMAXST);
    if (VER) {
      printf ("Kabsch (no renumbering) %f %f %f\n", RMSDST, DAVST, DMAXST);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1modk[i].El, coord1modk[i].xyz[0], coord1modk[i].xyz[1], coord1modk[i].xyz[2]);
      }
    }
    if (RMSD<RMSDST) {
      if (PRINT) {
        FPRINT (fileout, coord1mod, NATOMS, RMSD);
      }
      if (DISTRATIO==false) {
        printf ("%f %f %f\n", RMSD, DAV, DMAX);
      }
      else {
        long int UN = 0; // uniqiue
        if (DMAX>DISTLR) {
          UN = 1; 
        }
        else {
          if (DMAX>DISTSM) {
            // check for ratio
            double rat = DMAX/DAV;
            if (rat > DISTRT) {
              UN = 1;
            }
          }
        }
      printf ("%li %f %f\n", UN, DAV, DMAX);
      }
    }
    else {
      if (PRINT) {
        FPRINT (fileout, coord1modk, NATOMS, RMSDST);
      }
      if (DISTRATIO==false) {
        printf ("%f %f %f \n", RMSDST, DAVST, DMAXST);
      }
      else {
        long int UN = 0; // uniqiue
        if (DMAXST>DISTLR) {
          UN = 1; 
        }
        else {
          if (DMAXST>DISTSM) {
            // check for ratio
            double rat = DMAXST/DAV;
            if (rat > DISTRT) {
              UN = 1;
            }
          }
        }
      printf ("%li %f %f\n", UN, DAVST, DMAXST);
      }
    }
  }
  else if (strcmp(DIST, "angle") == 0) {
    double AAV=0;
    double AMAX=0;
    double AAVST=0;
    double AMAXST=0;
    KABSHANGLE (coord1best, coord2, &coord1mod, NAT1, &RMSD, &AAV, &AMAX);
    if (VER) {
      printf ("Kabsch (renumbering) %f %f %f\n", RMSD, AAV, AMAX);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
      }
    }
    KABSHANGLE (coord1, coord2, &coord1modk, NAT1, &RMSDST, &AAVST, &AMAXST);
    if (VER) {
      printf ("Kabsch (no renumbering) %f %f %f\n", RMSDST, AAVST, AMAXST);
      printf ("First mol\n");
      for (i=0;i<NAT1;i++){
        printf("%s %10.6f %10.6f %10.6f\n", coord1modk[i].El, coord1modk[i].xyz[0], coord1modk[i].xyz[1], coord1modk[i].xyz[2]);
      }
    }
    if (RMSD<RMSDST) {
        if (PRINT) {
          FPRINT (fileout, coord1mod, NATOMS, RMSD);
        }
        printf ("%f %f %f\n", RMSD, AAV, AMAX);
      }
    else {
        if (PRINT) {
          FPRINT (fileout, coord1modk, NATOMS, RMSDST);
        }
        printf ("%f %f %f \n", RMSDST, AAVST, AMAXST);
      }
    }
    else {
    printf ("Your --compare=%s is not clear.\n", DIST);
    exit(0);    
  }
  //printf ("RMSD = %f\n", RMSD);
}
else if (strcmp(ALG, "kabsch")==0) {
  // simple KABSH, no permutation
  if (VER) printf ("Standard Kabsch, no permutation.\n");
  double RMSD = 0;
  if (strcmp(DIST, "false") == 0) {
    KABSH (coord1, coord2, &coord1mod, NAT1, &RMSD);
    printf ("%f\n", RMSD);
  }
  else if (strcmp(DIST, "dist") == 0) {
    double DAV=0;
    double DMAX=0;
    KABSHDIST (coord1, coord2, &coord1mod, NAT1, &RMSD, &DAV, &DMAX);
      if (DISTRATIO==false) {
        printf ("%f %f %f\n", RMSD, DAV, DMAX);
      }
      else {
        long int UN = 0; // uniqiue
        if (DMAX>DISTLR) {
          UN = 1; 
        }
        else {
          if (DMAX>DISTSM) {
            // check for ratio
            double rat = DMAX/DAV;
            if (rat > DISTRT) {
              UN = 1;
            }
          }
        }
      printf ("%li %f %f\n", UN, DAV, DMAX);
      }
  }
  else if (strcmp(DIST, "angle") == 0) {
    double AAV=0;
    double AMAX=0;
    KABSHANGLE (coord1, coord2, &coord1mod, NAT1, &RMSD, &AAV, &AMAX);
    printf ("%f %f %f\n", RMSD, AAV, AMAX);
  }
  else {
    printf ("Your --compare=%s is not clear.\n", DIST);
    exit(0);    
  }
  if (VER) {
    printf ("First mol\n");
    for (i=0;i<NAT1;i++){
      printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
    }
  }
  if (PRINT) {
    FPRINT (fileout, coord1mod, NATOMS, RMSD);
  }
}

/*
printf ("First mol\n");
for (i=0;i<NAT1;i++){
  printf("%s %10.6f %10.6f %10.6f\n", coord1mod[i].El, coord1mod[i].xyz[0], coord1mod[i].xyz[1], coord1mod[i].xyz[2]);
}
printf ("Second mol\n");
for (i=0;i<NAT2;i++){
  printf("%s %10.6f %10.6f %10.6f\n", coord2[i].El, coord2[i].xyz[0], coord2[i].xyz[1], coord2[i].xyz[2]);
}
*/

// free all allocated structures & arrays

free(coord1);
free(coord2);
free(coorddef);
free(UnEl);
free(ElNum);
free(ElNum0);
free(ADJ1);
free(ADJ2);
free(IDN);
free(coord1mod);
free(coord1modk);
free(coord1best);
}

// convert string to lower case

char *strlwr(char *str)
{
  unsigned char *p = (unsigned char *)str;

  while (*p) {
     *p = tolower(*p);
      p++;
  }

  return str;
}

