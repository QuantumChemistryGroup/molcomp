#ifndef MC_INERTIA_H
#define MC_INERTIA_H

#include <stdbool.h>
#include "mc_types.h"

extern void MC_INERTIA(McAtom **coord, double II[3], long int NATOMS, double EF[3][3], bool VER);

#endif /* MC_INERTIA_H */
