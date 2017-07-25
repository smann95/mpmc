#pragma once

#include "Atom.h"
#include "Molecule.h"

// Each atom in a molecule has a pair list, which holds data relevant to each pair of atoms in the 
// system (so that it doesn't need to be recalculated for each simulation step--or perhaps a flag 
// indicating that the energetics for this pair *never* needs to be computed). If one of the pairs
// moves, the pair is flagged for recalculation, indicating that (some of) the stored data is no
// longer valid and needs to be recomputed.

// The first atom in the first molecule of the system will have a pair list with a pair entry for 
// every other atom in the system (across all other molecules). The 2nd atom in the system (i.e. 
// the 2nd atom in molecule 1, or [if "molecule" 1 only had 1 atom] atom 1 in molecule 2) will 
// have a pair list that contains a pair entry for every atom in the system, except atom 1. The 
// pairing of atom 1 with atom 2 is covered in the pair list of atom 1.

class Pair
{
public:
	Pair(){};
	~Pair(){};

	
	int      frozen,               //are they both MOF atoms, for instance
	         rd_excluded, 
	         es_excluded,
	         attractive_only,
	         recalculate_energy;
	double   lrc,                  //LJ long-range correction
	         last_volume,          //what was the volume when we last calculated LRC? needed for NPT
	         epsilon, sigma,       //LJ
	         r, rimg, dimg[3],     //separation and separation with nearest image
	         d_prev[3],            //last known position
	         rd_energy, 
	         es_real_energy,
	         es_self_intra_energy,
	         sigrep,
	         c6, c8, c10;
	Atom     * atom;               // the other atom in the pairing
	Molecule * molecule;           // the molecule to which the other atom belongs
	Pair     * next;
	
};

