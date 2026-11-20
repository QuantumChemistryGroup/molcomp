#ifndef MC_STDFUNC_H
#define MC_STDFUNC_H

#include <stdbool.h>
#include "mc_types.h"

extern bool MC_BONDEQUAL(long int p, long int q, McAtom *coord1, McAtom *coord2, long int NATOMS, bool MOL2);
extern void MC_POPADJ(long int NATOMS, long int ADJ[][NATOMS], McAtom *coord);
extern bool MC_RADIAL(long int p, long int q, McAtom *coord1, McAtom *coord2, long int NATOMS, bool MOL2, long int ADJ1[][NATOMS], long int ADJ2[][NATOMS]);

#endif /* MC_STDFUNC_H */
