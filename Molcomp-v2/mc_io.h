#ifndef MC_IO_H
#define MC_IO_H

#include <stdio.h>
#include <stdbool.h>
#include "mc_types.h"

extern void MC_FREAD(FILE *fp, McAtom **coordinates, long int *NATOMS, bool VER);
extern void MC_FPRINT(char *fileout, McAtom *coordinates, long int NATOMS, double RMSD);
extern void MC_TXYZREAD(FILE *fp, McAtom **coordinates, long int NATOMS);

#endif /* MC_IO_H */
