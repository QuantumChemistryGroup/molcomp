#ifndef MC_TYPES_H
#define MC_TYPES_H

#include <inttypes.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

typedef struct {
    // element
    char El[4];
    // x, y, z coordinates
    double xyz[3];
    // atom MM2 type
    long int t;
    // atom UFF type
    char *tu;
    // the array containing the number of atoms connected to this one
    long int n[23];
    long int nn; // the number of the neighbouring atoms
    double Z;
    double C6; // C6 coefficient for van der Waals
    double x;  // equilibrium distance for van der Waals
    long int N; // old sequential number
} McAtom;

/* Structure containing Mendeleev Table Data */
typedef struct {
    char El[4];
    long int Z;
    double A;
    double R;
    double RB;
    double RvdW;
    long int V;
    double M;
    double EN;
    double IP;
    double EA;
    char *ENAME;
} McMPT;

typedef struct {
    long int A[4]; // array containing the atom numbers
    double V;      // value of the distance or angle
} McMM;

typedef struct {
    bool VER;
    bool STEREO;       // --keepenantio
    bool CONNECT;      // --connect
    bool BINDH2HEAVY;  // --bindh2heavy
    bool MOL2;         // --mol2
    double LINEARTHLD; // --linearthld
    char ALG[16];      // "std", "adv", "kabsch"
} McOptions;

typedef struct {
    double rmsd;
    double dmax;
    double amax;
    int status;        // 0 = OK, < 0 = error
} McResult;

#endif /* MC_TYPES_H */
