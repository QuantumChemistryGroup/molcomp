bool BONDEQUAL (long int p, long int q, Atom *coord1, Atom *coord2, long int NATOMS, bool MOL2);
void POPADJ (long int NATOMS, long int ADJ[][NATOMS], Atom *coord);
bool RADIAL (long int p, long int q, Atom *coord1, Atom *coord2, long int NATOMS, bool MOL2, long int ADJ1[][NATOMS], long int ADJ2[][NATOMS]);
