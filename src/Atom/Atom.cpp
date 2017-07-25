#include "Atom.h"



Atom::Atom()
{
	
	id                   = 0;
	bond_id              = 0;
	atomtype[0]          = (char) 0;
	frozen               = 0;
	adiabatic            = 0;
	spectre              = 0;
	target               = 0;
	mass                 = 0.0;
	charge               = 0.0;
	polarizability       = 0.0;
	epsilon              = 0.0;
	sigma                = 0.0;
	omega                = 0.0;
	c6                   = 0;
	c8                   = 0;
	c10                  = 0;
	c9                   = 0;
	dipole_rrms          = 0.0;
	rank_metric          = 0.0;
	gwp_spin             = 0;
	gwp_alpha            = 0.0;
	site_neighbor_id     = 0; // dr fluctuations will be applied along the vector from this atom to the atom identified by this variable
	lrc_self             = 0.0;
	last_volume          = 0.0; // currently only used in disp_expansion.c
	es_self_point_energy = 0.0;

	for( int i=0; i<3; i++ ) {
		pos              [i] = 0.0;
		wrapped_pos      [i] = 0.0; //absolute and wrapped (into main unit cell) position
		ef_static        [i] = 0.0;
		ef_static_self   [i] = 0.0;
		ef_induced       [i] = 0.0;
		ef_induced_change[i] = 0.0;
		mu               [i] = 0.0;
		old_mu           [i] = 0.0;
		new_mu           [i] = 0.0;
	}
	

	// pair_t *pairs;

}


Atom::~Atom()
{
}
