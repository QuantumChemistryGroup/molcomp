#ifndef MC_CORE_H
#define MC_CORE_H

#include <stdbool.h>
#include "mc_types.h"

typedef struct McPrecomp McPrecomp;

extern McOptions MC_DEFAULT_OPTIONS(void);

extern McResult MC_COMPARE(
    const McAtom *coord1_in,
    const McAtom *coord2_in,
    long int NATOMS,
    const McOptions *opts,
    bool skip_connect1,
    bool skip_connect2,
    McAtom *coord1best_out
);

extern McResult MC_COMPARE_FLAT(
    long int NATOMS,
    const char (*elements)[ 4 ],
    const long int *nn,
    const long int (*nbrs)[ 23 ],
    const double *xyz1,
    const double *xyz2,
    const McOptions *opts
);

/* Precomputed molecular topology context (computed once for all conformers) */
extern McPrecomp *MC_PRECOMP_CREATE(
    long int NATOMS,
    const char (*elements)[ 4 ],
    const double *xyz_ref,
    const long int *nn_in,
    const long int (*nbrs_in)[ 23 ],
    const McOptions *opts
);

extern void MC_PRECOMP_FREE(McPrecomp *pc);

/* Ultra-fast comparison using precomputed IDN, masses, and torsions (thread-safe) */
extern McResult MC_COMPARE_PRECOMP(
    const McPrecomp *pc,
    const double *xyz1,
    const double *xyz2,
    const McOptions *opts
);

#endif /* MC_CORE_H */
