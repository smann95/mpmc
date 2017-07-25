#pragma once
#include <random>

#include "constants.h"
#include "Atom.h"
class PeriodicBoundary;

class Molecule
{
public:

	Molecule();
	Molecule( const  Molecule &orig );

	~Molecule();

	void free_pairs();
	void rotate_rand_pbc(    double scale, const PeriodicBoundary &pbc, std::mt19937 *mt_rand );
	void translate_rand_pbc( double scale, const PeriodicBoundary &pbc, std::mt19937 *mt_rand );
	void displace_gwp( double scale, std::mt19937 *mt_rand );
	
	int        id;
	char       moleculetype[maxLine];
	double     mass;
	int        frozen, 
	           adiabatic,
	           spectre,
	           target;
	double     com[3],
	           wrapped_com[3]; //center of mass
	double     iCOM       [3]; // initial Center of Mass
	int        nuclear_spin;
	double     rot_partfunc_g, 
	           rot_partfunc_u,
	           rot_partfunc;
	Atom     * atoms;
	// vector<Atom> atoms; ?
	Molecule * next;


	#ifdef QM_ROTATION
		double     * quantum_rotational_energies;
		complex_t ** quantum_rotational_eigenvectors;
		int        * quantum_rotational_eigensymmetry;
		double       quantum_rotational_potential_grid[QUANTUM_ROTATION_GRID][QUANTUM_ROTATION_GRID];
	#endif // QM_ROTATION
	#ifdef XXX
		// XXX - vib work in progress
		double     * quantum_vibrational_energies;
		complex_t ** quantum_vibrational_eigenvectors;
		int        * quantum_vibrational_eigensymmetry;
	#endif // XXX
};

