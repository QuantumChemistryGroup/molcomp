#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "maintypes.h"
#include "fread.h"
#include <ctype.h>

//int isalpha(int c);

int isBlank (char * line)
{
  char * ch;
  int is_blank = -1;

  // Iterate through each character.
  for (ch = line; *ch != '\0'; ++ch)
  {
    if (!isspace(*ch))
    {
      // Found a non-whitespace character.
      is_blank = 0;
      break;
    }
  }

  return is_blank;
}
void FREAD(FILE * fp, Atom **coordinates, long int *NATOMS) {
rewind(fp);
// the length of the array of char
long int N=0;
// reading the whole file and putting it into array of char
char **list = LINES(fp, &N);
//printf("I am here %li %li\n", N, *NATOMS);

long int i;
long int NA = 0; // natoms given
long int j = 0;
for (i=0;i<N;i++) {
  char s[300];
  strcpy(s, list[i]);
  //printf ("%li %s %li %li\n", i, list[i], *NATOMS, j);
  if (i == 0) {
    NA = atoi(s);
    if (NA == 0) {
      if (VER) printf ("Warning: NATOMS in your xyz at first line is either 0 or missing.\n");
    }
    else {
      //continue;
      i=1;
      //printf ("NATOMS = %li\n", NA);
    }	
  }
  //printf ("%li %s %li %li\n", i, (*coordinates)[(*NATOMS)].El, *NATOMS, j);
//	strlwr(s);
  if (i>=0) {
    if (j==0) {
    // check if this contains only alphabet
      if (isalpha(list[i][0])) {
        strcpy((*coordinates)[(*NATOMS)].El, list[i]);
        //printf ("%s %li %li\n", (*coordinates)[(*NATOMS)].El, *NATOMS, j);
        j++;
      }
      else {
        j=0;
      }
    }
    else if (j>0 && j<4) {
      // if there is a character -> we stop the code
      if (isalpha(list[i][0])) {
        // stop the code
        printf ("Something is wrong with your coordinates.\n");
        printf ("Stop the code\n");
        printf ("%s\n", list[i]);
        exit(0);
      }
      else {
        (*coordinates)[*NATOMS].xyz[j-1] = atof(list[i]);
        //printf ("%f... %li %li\n", (*coordinates)[*NATOMS].xyz[j-1], *NATOMS, j);
        j++;
      }
    }
    if (j==4) {
      (*NATOMS)++;
      // reallocate
      (*coordinates) = realloc((*coordinates), ((*NATOMS)+1)*sizeof(Atom));
      // set j to 0
      j=0;
    }
    
  }
}
if (NA != (*NATOMS)) {
  if (NA == 0) {
    if (VER) printf ("Warning: No number of atoms in the first line.\n");
  }
  else {
    printf ("NA = %li NATOMS = %li\n", NA, (*NATOMS));
    printf ("Stop the code\n");
    exit(0);
  }
}
for (i=0;i<N+1;i++) {
  free(list[i])	;
}

free(list);

}

// the function which reads the file and puts everything into char array

char **LINES (FILE *fp, long int *N) {
char **list = malloc(1*sizeof(char*));
list[0] = malloc(100*sizeof(char));
//strcpy((*lines), *start);
// the string length
char str[300];
// ptr to a string
char *ptr;
long int IT =0;
bool Skip2line = false;
while( fgets(str, sizeof(str), fp) != NULL ){
	// just skip the second line - molecule name
	  bool Skip= false;
	  if ((Skip2line == true) && (IT == 1))
	  	{
	  		Skip=true;
	  	}
	  if (Skip == false) 
	  {
	if (isBlank(str) == 0) {
	ptr = strtok(str, " \n\t\r");
	while (ptr !=NULL) {
		strcpy((list[*N]), ptr);
		// check if the first symbol is the number of atoms
		if (IT == 0 && (*N) == 0) {
		  if (atoi (ptr) !=0) {
		    Skip2line = true;
		  }
		}
		(*N)++;
		list = realloc(list, ((*N)+1)*sizeof(char*) );
		list[(*N)] = malloc(100*sizeof(char));
		ptr = strtok(NULL, " \n\t\r");
	}
	
	}
	else {
		strcpy((list[*N]), "\n");
		(*N)++;
		list = realloc(list, ((*N)+1)*sizeof(char*) );
		list[(*N)] = malloc(100*sizeof(char));
	}
	}
	IT=IT+1;
}
//printf ("Everything is copied from the input file\n");
return list;
}

// function to print modified first structure
void FPRINT(char *fileout, Atom *coordinates, long int NATOMS, double RMSD) {
FILE *fw;
fw = fopen(fileout, "w");
long int i=0;
fprintf(fw, " %li\n", NATOMS);
fprintf(fw, " RMSD = %f\n", RMSD);
for (i=0;i<NATOMS;i++) {
  fprintf(fw, "%3s %12.8f %12.8f %12.8f\n", coordinates[i].El, coordinates[i].xyz[0], coordinates[i].xyz[1], coordinates[i].xyz[2]);  
}
fclose(fw);
}
