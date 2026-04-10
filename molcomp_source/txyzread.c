#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "maintypes.h"
#include <string.h>
#include <ctype.h>
#include "maintypes.h"
#include "txyzread.h"

void TXYZREAD (FILE *fp, Atom **coordinates, long int NATOMS) {

// the string length
char str[300];
// ptr to a string
char *ptr;
long int i=0, index=0;

// for some strange reasons we need to read the first line twice
if (fgets(str, 300, fp)!=NULL) {
	//puts(str);
}
ptr = strtok(str, " ");
ptr = strtok(NULL, " ");
ptr = strtok(NULL, " ");
char *T;
T= ptr;
//T = strdup(ptr);
if (fgets(str, 300, fp)!=NULL) {
	//puts(str);
}
// get the index=0
// read the atoms and out them into array
for (i=0; i< NATOMS; i++) {
	// set to 0
	(*coordinates)[i].nn = 0;
	ptr = strtok(str, " ");         // split our findings around the " "
	//printf( "%2i %2s", index,"This is the number of lines\n");
	index = 0;
	while(ptr != NULL)  // while there's more to the string
	{
		// when splitting every string, the first element is the atom number, etc. We put it into structure.    
		if (index == 5 && strcmp(T, "MM2") == 0){
			(*coordinates)[i].t = atoi(ptr);
			}
		else {
			(*coordinates)[i].tu = ptr;
			// hard copy
			//(*coordinates)[i].tu = strdup((*coordinates)[i].tu);
		}
		// check whether we do not deal with the blank line	
		if (index > 5 && atoi(ptr) !=0){
			(*coordinates)[i].nn = index-5;
			(*coordinates)[i].n[index-6] = atoi(ptr);
			//printf(" %li ? \n", atoi(ptr));
			}
		ptr = strtok(NULL, " ");                
		// and keep splitting
		index++;
	}
	// read another line            
	if (fgets(str, 300, fp)!=NULL) {
		//puts(str);
	}
}

}

