#include <inttypes.h>

#ifndef MAINTYPES_H
#define MAINTYPES_H

extern bool VER;
extern bool MOL2;
extern double LINEARTHLD;
static const double M_PI = 3.14159265358979323846264338327950288;

typedef struct {
        // element
        char El[3];
        // x, y, z coordinates
        double xyz[3];
        // atom MM2 type
        long int t;
        // atom UFF type
        char *tu;    
        // the array containing the number of atoms connected to this one
        // the first element gives the overall number of atoms connected.
        long int n[23]; // this is the array. The max number of connected elements is specified.
        long int nn; // the number of the neighbouring atoms
        double Z;
        double C6; // C6 coefficient for van der Waals
        double x; // equilibrium distance for van der Waals
        long int N; // old sequential number
} Atom;

/* Structure containing Mendeleev Table Data*/

typedef struct {
        // Element name
        char El[3];
        // Z
        long int Z;
        // ARENeg
        double A;
        // Covalent radii
        double R;
        // RBO
        double RB;
        // Van der Waals Radii
        double RvdW;
        // Maximum valency (coordination number) or maximum nuber of bonds
        long int V;
        // Mass
        double M;
        // Electronegativity
        double EN;
        // Ionization Potential
        double IP;
        // Electron Affinity
        double EA;
        // Element full name
        char *ENAME;
} MPT;

typedef struct {
        long int A[4]; // array containing the atom numbers
        double V; // value of the distance or angle
} MM;


#endif
