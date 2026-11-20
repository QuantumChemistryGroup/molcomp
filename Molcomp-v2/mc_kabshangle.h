#ifndef MC_KABSHANGLE_H
#define MC_KABSHANGLE_H

#include <stdbool.h>
#include "mc_types.h"

extern void MC_KABSHANGLE(
    McAtom *coord1,
    McAtom *coord2,
    McAtom **coord1mod,
    long int NATOMS,
    double *RMSD,
    double *AAV,
    double *AMAX,
    bool MOL2,
    double LINEARTHLD,
    bool VER
);

/* Fast version: uses precomputed torsions and flat arrays; skips torsion loop if compute_angles == false */
extern void MC_KABSHANGLE_FAST(
    const double *xyz1,
    const double *xyz2,
    double *xyz1_aligned,
    long int NATOMS,
    const McMM *tors,
    long int ntors,
    bool compute_angles,
    double *RMSD,
    double *DMAX,
    double *AMAX,
    double LINEARTHLD,
    bool VER
);

#endif /* MC_KABSHANGLE_H */


