#ifndef MC_GETCONNECT_H
#define MC_GETCONNECT_H

#include "mc_types.h"

extern void MC_GETCONNECT(McAtom **coordinates, long int NATOMS);
extern void MC_SYMMETRYCHECK(McAtom **coordinates, long int NATOMS);
extern void MC_SYMMETRYCHECKDEF(McAtom *coorddef, McAtom **coord1, long int NATOMS);

#endif /* MC_GETCONNECT_H */
