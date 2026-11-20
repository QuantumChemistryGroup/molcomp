#ifndef MC_RENUMB_HUNG_H
#define MC_RENUMB_HUNG_H

#include <stdbool.h>
#include "mc_types.h"

extern void MC_RENUMBHUNGCONNECT(McAtom **coord1tmp, McAtom *coord2, long int NATOMS, char (*UnEl)[4], long int UE, bool IDN[][NATOMS], bool MOL2);
extern void MC_RENUMBHUNG(McAtom **coord1tmp, McAtom *coord2, long int NATOMS, char (*UnEl)[4], long int UE);
extern void MC_RENUMBHUNGGRAPH(McAtom **coord1tmp, McAtom *coord2, long int NATOMS, char (*UnEl)[4], long int UE, bool IDN[][NATOMS]);

#endif /* MC_RENUMB_HUNG_H */
