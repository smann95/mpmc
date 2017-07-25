#include <random>
#include <cstring>

#include "Molecule.h"
#include "PeriodicBoundary.h"
#include "Pair.h"
#include "Quaternion.h"
#include "SafeOps.h"


Molecule::Molecule()
{
	id = 0;
	moleculetype[0] = (char) 0;
	mass            = 0.0;
	frozen          = 0; 
	adiabatic       = 0;
	spectre         = 0;
	target          = 0;
	nuclear_spin    = 0;
	rot_partfunc_g  = 0;
	rot_partfunc_u  = 0;
	rot_partfunc    = 0;
	for( int i=0; i<3; i++ ) {
		com        [i] = 0.0;
	    wrapped_com[i] = 0.0;  //center of mass
	    iCOM       [i] = 0.0;  // initial Center of Mass
	}

	#ifdef QM_ROTATION
		double     *quantum_rotational_energies      = nullptr;
		complex_t **quantum_rotational_eigenvectors  = nullptr;
		int        *quantum_rotational_eigensymmetry = nullptr;
		double      quantum_rotational_potential_grid[QUANTUM_ROTATION_GRID][QUANTUM_ROTATION_GRID] = {0};
	#endif // QM_ROTATION
	#ifdef XXX
		// XXX - vib work in progress
		double     *quantum_vibrational_energies      = nullptr;
		complex_t **quantum_vibrational_eigenvectors  = nullptr;
		int        *quantum_vibrational_eigensymmetry = nullptr;
	#endif // XXX
}


Molecule::Molecule( const Molecule &other ) {

	Atom *atom_dst_ptr       = nullptr,
	     *prev_atom_dst_ptr  = nullptr,
	     *atom_src_ptr       = nullptr;
	Pair *pair_dst_ptr       = nullptr, 
	     *prev_pair_dst_ptr  = nullptr,
	     *pair_src_ptr       = nullptr;

	// copy molecule attributes
	id             = other.id;
	mass           = other.mass;
	frozen         = other.frozen;
	adiabatic      = other.adiabatic;
	spectre        = other.spectre;
	target         = other.target;
	nuclear_spin   = other.nuclear_spin;
	rot_partfunc_g = other.rot_partfunc_g;
	rot_partfunc_u = other.rot_partfunc_u;
	rot_partfunc   = other.rot_partfunc;

	strcpy( moleculetype, other.moleculetype );
	
	for( int i=0; i<3; i++ ) {
		com        [i] = other.com        [i];
		wrapped_com[i] = other.wrapped_com[i];
	}

	#ifdef QM_ROTATION
		// untouched--needs to be converted - bt
		int i, j;
		if(system->quantum_rotation) {
			dst->quantum_rotational_energies = calloc(system->quantum_rotation_level_max, sizeof(double));
			memnullcheck(dst->quantum_rotational_energies, system->quantum_rotation_level_max*sizeof(double),__LINE__-1, __FILE__);
			dst->quantum_rotational_eigenvectors = calloc(system->quantum_rotation_level_max, sizeof(complex_t *));
			memnullcheck(dst->quantum_rotational_eigenvectors,system->quantum_rotation_level_max*sizeof(complex_t *),__LINE__-1, __FILE__);
			for(i = 0; i < system->quantum_rotation_level_max; i++) {
				dst->quantum_rotational_eigenvectors[i] = calloc((system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1), sizeof(complex_t));
				memnullcheck(dst->quantum_rotational_eigenvectors[i],(system->quantum_rotation_l_max+1)*(system->quantum_rotation_l_max+1)*sizeof(complex_t),__LINE__-1, __FILE__);
			}
			dst->quantum_rotational_eigensymmetry = calloc(system->quantum_rotation_level_max, sizeof(int));
			memnullcheck(dst->quantum_rotational_eigensymmetry,system->quantum_rotation_level_max*sizeof(int),__LINE__-1, __FILE__);

			memcpy(dst->quantum_rotational_energies, src->quantum_rotational_energies, system->quantum_rotation_level_max*sizeof(double));
			for(i = 0; i < system->quantum_rotation_level_max; i++) {
				for(j = 0; j < (system->quantum_rotation_l_max + 1)*(system->quantum_rotation_l_max + 1); j++) {
					dst->quantum_rotational_eigenvectors[i][j].real = src->quantum_rotational_eigenvectors[i][j].real;
					dst->quantum_rotational_eigenvectors[i][j].imaginary = src->quantum_rotational_eigenvectors[i][j].imaginary;
				}
			}
			memcpy(dst->quantum_rotational_eigensymmetry, src->quantum_rotational_eigensymmetry, system->quantum_rotation_level_max*sizeof(int));
		}
	#endif // QM_ROTATION
	next = nullptr;


	// new atoms list
	SafeOps::calloc( atoms, 1, sizeof(Atom), __LINE__, __FILE__ );
	prev_atom_dst_ptr = atoms;


	for( atom_dst_ptr = atoms, atom_src_ptr = other.atoms; atom_src_ptr; atom_dst_ptr = atom_dst_ptr->next, atom_src_ptr = atom_src_ptr->next ) {

		strcpy(atom_dst_ptr->atomtype, atom_src_ptr->atomtype);
		atom_dst_ptr->id             = atom_src_ptr->id;
		atom_dst_ptr->frozen         = atom_src_ptr->frozen;
		atom_dst_ptr->adiabatic      = atom_src_ptr->adiabatic;
		atom_dst_ptr->spectre        = atom_src_ptr->spectre;
		atom_dst_ptr->target         = atom_src_ptr->target;
		atom_dst_ptr->mass           = atom_src_ptr->mass;
		atom_dst_ptr->charge         = atom_src_ptr->charge;
		atom_dst_ptr->gwp_alpha      = atom_src_ptr->gwp_alpha;
		atom_dst_ptr->gwp_spin       = atom_src_ptr->gwp_spin;
		atom_dst_ptr->polarizability = atom_src_ptr->polarizability;
		atom_dst_ptr->omega          = atom_src_ptr->omega;
		atom_dst_ptr->epsilon        = atom_src_ptr->epsilon;
		atom_dst_ptr->sigma          = atom_src_ptr->sigma;
		atom_dst_ptr->gwp_spin       = atom_src_ptr->gwp_spin;
		atom_dst_ptr->c6             = atom_src_ptr->c6;
		atom_dst_ptr->c8             = atom_src_ptr->c8;
		atom_dst_ptr->c10            = atom_src_ptr->c10;
		atom_dst_ptr->c9             = atom_src_ptr->c9;

		memcpy( atom_dst_ptr->pos,         atom_src_ptr->pos,         3*sizeof(double) );
		memcpy( atom_dst_ptr->wrapped_pos, atom_src_ptr->wrapped_pos, 3*sizeof(double) );
		memcpy( atom_dst_ptr->ef_static,   atom_src_ptr->ef_static,   3*sizeof(double) );
		memcpy( atom_dst_ptr->ef_induced,  atom_src_ptr->ef_induced,  3*sizeof(double) );
		memcpy( atom_dst_ptr->mu,          atom_src_ptr->mu,          3*sizeof(double) );
		memcpy( atom_dst_ptr->old_mu,      atom_src_ptr->old_mu,      3*sizeof(double) );
		memcpy( atom_dst_ptr->new_mu,      atom_src_ptr->new_mu,      3*sizeof(double) );

		SafeOps::calloc( atom_dst_ptr->pairs, 1, sizeof(Pair), __LINE__, __FILE__ );
		pair_dst_ptr = atom_dst_ptr->pairs;
		prev_pair_dst_ptr = pair_dst_ptr;
		for(pair_src_ptr = atom_src_ptr->pairs; pair_src_ptr; pair_src_ptr = pair_src_ptr->next) {

			pair_dst_ptr->rd_energy            = pair_src_ptr->rd_energy;
			pair_dst_ptr->es_real_energy       = pair_src_ptr->es_real_energy;
			pair_dst_ptr->es_self_intra_energy = pair_src_ptr->es_self_intra_energy;


			pair_dst_ptr->frozen               = pair_src_ptr->frozen;
			pair_dst_ptr->rd_excluded          = pair_src_ptr->rd_excluded;
			pair_dst_ptr->es_excluded          = pair_src_ptr->es_excluded;
//			pair_dst_ptr->charge               = pair_src_ptr->charge;
//			pair_dst_ptr->gwp_alpha            = pair_src_ptr->gwp_alpha;
//			pair_dst_ptr->gwp_spin             = pair_src_ptr->gwp_spin;
			pair_dst_ptr->epsilon              = pair_src_ptr->epsilon;
			pair_dst_ptr->lrc                  = pair_src_ptr->lrc;
			pair_dst_ptr->sigma                = pair_src_ptr->sigma;
			pair_dst_ptr->r                    = pair_src_ptr->r;
			pair_dst_ptr->rimg                 = pair_src_ptr->rimg;

			SafeOps::calloc( pair_dst_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__ );
			prev_pair_dst_ptr = pair_dst_ptr;
			pair_dst_ptr = pair_dst_ptr->next;

		}
		prev_pair_dst_ptr->next = nullptr;
		free(pair_dst_ptr);
		// handle an empty list
		if( !atom_src_ptr->pairs )
			atom_dst_ptr->pairs = nullptr;

		prev_atom_dst_ptr = atom_dst_ptr;
		SafeOps::calloc( atom_dst_ptr->next, 1, sizeof(Atom), __LINE__, __FILE__ );
	}

	prev_atom_dst_ptr->next = nullptr;
	free(atom_dst_ptr);
}


// free a molecule and all of it's associated stuff
Molecule::~Molecule() {
	
	#ifdef QM_ROTATION
		// untouched--needs to be converted - bt
		int i;
		if(system->quantum_rotation && !molecule->frozen) {

			free(molecule->quantum_rotational_energies);
			free(molecule->quantum_rotational_eigensymmetry);
			for(i = 0; i < system->quantum_rotation_level_max; i++)
				free(molecule->quantum_rotational_eigenvectors[i]);
			free(molecule->quantum_rotational_eigenvectors);

		}
	#endif // QM_ROTATION

	// free this molecule's pairs	
	Pair  ** ptr_array = nullptr,
	       * pptr      = nullptr;
	Atom   * aptr      = nullptr,
	      ** aarray    = nullptr;
	int      i         = 0;


	//build an array of pairs
	for( aptr = atoms; aptr; aptr = aptr->next ) {
		for ( pptr = aptr->pairs; pptr; pptr = pptr->next ) {
			SafeOps::realloc( ptr_array, sizeof( Pair *) * (i+1), __LINE__, __FILE__ );
			ptr_array[i] = pptr;
			++i;
		}
	}

	// free the pairs
	for( --i; i>=0; i-- )
		free(ptr_array[i]);

	// zero out the heads
	for( aptr = atoms; aptr; aptr = aptr->next )
		aptr->pairs = nullptr;

	// free the temp array
	free(ptr_array);
	

	// free_my_atoms(molecule);
	aptr           = nullptr;
	i              = 0;
		
	// build an array of atoms
	for ( aptr = atoms; aptr; aptr = aptr->next ) {
		SafeOps::realloc( aarray, sizeof(Atom *)*(i+1), __LINE__, __FILE__ );
		aarray[i] = aptr;
		i++;
	}

	//free the atoms
	while ( i-- )
		free(aarray[i]);

	//free the temp array
	free(aarray);
}


// perform a general random rotation now with quaternions - AH 
void Molecule::rotate_rand_pbc( double scale, const PeriodicBoundary &pbc, std::mt19937 *mt_rand ) {

	//std::uniform_real_distribution<double> dist{0,1};
	std::normal_distribution<double> dist{0.5,1.0};

	Atom   * atom_ptr        = nullptr;
	double * new_coord_array = nullptr,
	         x               =  dist( *mt_rand )  -  0.5,           // create a random axis and random angle to rotate around
	         y               =  dist( *mt_rand )  -  0.5,           // constructor will ensure axis is a unit vector
	         z               =  dist( *mt_rand )  -  0.5,
	         angle           =  dist( *mt_rand )  *  360  *  scale; // random angle, can be affected by scale 
	         
	Quaternion rnd_rotation(x,y,z,angle, Quaternion::AXIS_ANGLE_DEGREE );
	Quaternion rnd_rotation_conjugate = rnd_rotation.conjugate(); // needed to transform coordinates

	// count the number of atoms in a molecule, and allocate new coords array 
	int n=0;
	for( atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next)
		++n;
	SafeOps::calloc( new_coord_array, n*3, sizeof(double), __LINE__, __FILE__ );
	
	// translate the molecule to the origin 
	for(atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] -= com[0];
		atom_ptr->pos[1] -= com[1];
		atom_ptr->pos[2] -= com[2];
	}

	// quaternion multiply
	atom_ptr = atoms;
	for( int i=0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		int ii = i*3;
    
		//position_vector = position
		Quaternion position_vector( atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2], 0.0, Quaternion::XYZW );
		//quaternion_construct_xyzw(&position_vector,atom_ptr->pos[0],atom_ptr->pos[1],atom_ptr->pos[2],0.);
		Quaternion answer = position_vector * rnd_rotation_conjugate;
		answer = rnd_rotation * answer;
		//answer = rnd_rotation*(position*rnd_rotation_conjugate)
		//quaternion_multiplication(&position_vector,&rnd_rotation_conjugate,&answer);
		//quaternion_multiplication(&rnd_rotation,&answer,&answer);

		//set the new coords
		new_coord_array[ii+0] = answer.X();
		new_coord_array[ii+1] = answer.Y();
		new_coord_array[ii+2] = answer.Z();
	}

	// set the new coordinates and then translate back from the origin
	atom_ptr = atoms;
	for( int i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		int ii = i*3;
		atom_ptr->pos[0] = new_coord_array[ii+0];
		atom_ptr->pos[1] = new_coord_array[ii+1];
		atom_ptr->pos[2] = new_coord_array[ii+2];

		atom_ptr->pos[0] += com[0];
		atom_ptr->pos[1] += com[1];
		atom_ptr->pos[2] += com[2];

	}

	// free our temporary array
	free(new_coord_array);

}



// perform a general random translation  (formerly "translate(...)" )
void Molecule::translate_rand_pbc( double scale,  const PeriodicBoundary &pbc, std::mt19937 *mt_rand ) {

	std::uniform_real_distribution<double> dist{0,1};

	Atom   * atom_ptr;
	double   trans_x = 0,
	         trans_y = 0,
	         trans_z = 0;

	trans_x = scale * dist( * mt_rand ) * pbc.cutoff;
	trans_y = scale * dist( * mt_rand ) * pbc.cutoff;
	trans_z = scale * dist( * mt_rand ) * pbc.cutoff;
	if( dist( * mt_rand ) < 0.5) trans_x *= -1.0;
	if( dist( * mt_rand ) < 0.5) trans_y *= -1.0;
	if( dist( * mt_rand ) < 0.5) trans_z *= -1.0;

	com[0] += trans_x;
	com[1] += trans_y;
	com[2] += trans_z;

	for(atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] += trans_x;
		atom_ptr->pos[1] += trans_y;
		atom_ptr->pos[2] += trans_z;
	}
}



// perform a perturbation to a gaussian width
void Molecule::displace_gwp( double scale, std::mt19937 *mt_rand ) {

	std::uniform_real_distribution<double> dist{0,1};

	Atom   * atom_ptr;
	double   perturb  = 0;

	for(atom_ptr = atoms; atom_ptr; atom_ptr = atom_ptr->next) {

		if(atom_ptr->gwp_spin) {
			perturb = scale * ( dist(* mt_rand) - 0.5);
			atom_ptr->gwp_alpha += perturb;
			// make sure the width remains positive 
			atom_ptr->gwp_alpha = fabs(atom_ptr->gwp_alpha);
		}
	}
}