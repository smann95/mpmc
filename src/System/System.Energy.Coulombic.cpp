#include "Atom.h"
#include "Molecule.h"
#include "Output.h"
#include "Pair.h"
#include "System.h"
#include "UsefulMath.h"

// No ewald summation - regular accumulation of Coulombic terms without out consideration of PBC
// Only used by surface module 
double System::coulombic_nopbc( Molecule * molecules ) {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double pe       = 0;
	double total_pe = 0;

	total_pe = 0;
	for(  molecule_ptr = molecules;   molecule_ptr;   molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms;   atom_ptr;   atom_ptr = atom_ptr->next ) {
			for( pair_ptr = atom_ptr->pairs;   pair_ptr;   pair_ptr = pair_ptr->next ) {
				if( ! pair_ptr->es_excluded ) {
					pe = atom_ptr->charge * pair_ptr->atom->charge / pair_ptr->r;
					total_pe += pe;
				}
			}
		}
	}

	return total_pe;
}



double System::coulombic_nopbc_gwp() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double pe       = 0, 
	       total_pe = 0;
	double qi       = 0,
	       qj       = 0,
	       ai       = 0,
	       aj       = 0,
	       r        = 0;

	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {

				r  = pair_ptr->rimg;
				qi = atom_ptr->charge;
				qj = pair_ptr->atom->charge;
				ai = atom_ptr->gwp_alpha;
				aj = pair_ptr->atom->gwp_alpha;

				if(atom_ptr->gwp_spin || pair_ptr->atom->gwp_spin) {
					pe = qi * qj *  erf(  sqrt(3.0/2.0*(ai*ai+aj*aj)) * r  )  /  r;
				} else {
					pe = qi*qj/r;
				}

				total_pe += pe;
			}
		}
	}

	return total_pe;
}



double System::coulombic_kinetic_gwp() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	double     ai      = 0, 
	           mass    = 0,
	           energy  = 0;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			if( atom_ptr->gwp_spin ) {
				ai      = atom_ptr->gwp_alpha;
				mass    = atom_ptr->mass;
				energy += 9.0 * hBar*hBar  /  (8.0* (ai/METER2ANGSTROM) * (ai/METER2ANGSTROM) * (AMU2KG*mass))  /  kB;
			}

		} // atom
	} // molecule

	return(energy);
}


// total ES energy term 
double System::coulombic() {

	double real       = 0,
	       reciprocal = 0,
	       self       = 0,
	       potential  = 0;

	// construct the relevant ewald terms 
	if( wolf )
		potential  = coulombic_wolf();
	else {
		real       = coulombic_real();
		reciprocal = coulombic_reciprocal();
		self       = coulombic_self();

		// return the total electrostatic energy
		potential = real + reciprocal + self;
	}

	return potential;
}



double System::coulombic_wolf() {

	Molecule * mptr;
	Atom     * aptr;
	Pair     * pptr;

	double pot        = 0,
	       alpha      = ewald_alpha,
	       R          = pbc.cutoff,
	       iR         = 1.0/R,
	       erfaRoverR = erf(alpha*R)/R,
	       r          = 0,
	       ir         = 0;

	for( mptr = molecules; mptr; mptr = mptr->next ) {
		for( aptr = mptr->atoms; aptr; aptr = aptr->next ) {
			for( pptr = aptr->pairs; pptr; pptr = pptr->next ) {

				if ( pptr->recalculate_energy ) {
					pptr->es_real_energy = 0;

					r = pptr->rimg;
					ir = 1.0/r;
					if( (!pptr->frozen) && (!pptr->es_excluded) && (r < R) ) {
						pptr->es_real_energy = 
							aptr->charge * pptr->atom->charge * (ir-erfaRoverR-iR*iR*(R-r));
									
						// get feynman-hibbs contribution
						if( feynman_hibbs ) {
							Output::err( "COULOMBIC: FH + es_wolf is not implemented\n" );
							throw incompatible_settings;
						} 
					}  // r<cutoff
				} //recalculate

				pot += pptr->es_real_energy;

			} //pair
		} //atom
	} //molecule

	return(pot);
}


// real space sum 
double System::coulombic_real() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	Pair     * pair_ptr     = nullptr;

	double alpha               = ewald_alpha,
	       r                   = 0,
	       erfc_term           = 0,
	       gaussian_term       = 0,
	       potential           = 0,
	       potential_classical = 0;


	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {

				if( pair_ptr->recalculate_energy ) {
					pair_ptr->es_real_energy = 0;

					if( ! pair_ptr->frozen ) {

						r = pair_ptr->rimg;
						if(  !(   (r > pbc.cutoff) || pair_ptr->es_excluded   )  ) { // unit cell part

							//calculate potential contribution
							erfc_term = erfc(alpha*r);
							gaussian_term = exp(-alpha*alpha*r*r);
							potential_classical = atom_ptr->charge * pair_ptr->atom->charge * erfc_term/r;
							//store for pair pointer, so we don't always have to recalculate
							pair_ptr->es_real_energy += potential_classical;

							if( feynman_hibbs )
								pair_ptr->es_real_energy += coulombic_real_FH( molecule_ptr, pair_ptr, gaussian_term, erfc_term );		

						} else if(pair_ptr->es_excluded) // calculate the charge-to-screen interaction
							pair_ptr->es_self_intra_energy  =  atom_ptr->charge * pair_ptr->atom->charge * erf(alpha*pair_ptr->r)  /  pair_ptr->r;

					} // frozen 
				} // recalculate 

				// sum all of the pairwise terms
				potential += pair_ptr->es_real_energy - pair_ptr->es_self_intra_energy;

			} // pair
		} // atom
	} // molecule

	return potential;
}


// feynman-hibbs for real space
double System::coulombic_real_FH( Molecule * molecule_ptr, Pair *pair_ptr, double gaussian_term, double erfc_term ) {

	double du           = 0,
	       d2u          = 0,
	       d3u          = 0,
	       d4u          = 0; //derivatives of the pair term
	double fh_2nd_order = 0;
	double fh_4th_order = 0;
	double r            = pair_ptr->rimg;
	double rr           = r*r;
	double ir           = 1.0/r;
	double ir2          = ir*ir;
	double ir3          = ir*ir2;
	double ir4          = ir2*ir2;
	double order        = feynman_hibbs_order;
	double alpha        = ewald_alpha;
	double a2           = alpha*alpha;
	double a3           = a2*alpha;
	double a4           = a3*alpha;
	double reduced_mass = AMU2KG*molecule_ptr->mass*pair_ptr->molecule->mass/(molecule_ptr->mass+pair_ptr->molecule->mass);

	du  = -2.0 * alpha * gaussian_term / (r*sqrt(pi)) - erfc_term * ir2;
	d2u = (4.0/sqrt(pi))*gaussian_term*(a3 + 1.0*ir2) + 2.0*erfc_term*ir3;

	fh_2nd_order = (M2A2)  *  (hBar2/(24.0 * kB * temperature * reduced_mass))  *  (d2u + 2.0*du/r);

	if( order >= 4 ) {

		d3u = (gaussian_term/sqrt(pi))*(-8.0*(a3*a2)*r - 8.0*(a3)/r - 12.0*alpha*ir3) - 6.0*erfc(alpha*r)*ir4;
		d4u = (gaussian_term/sqrt(pi))*( 8.0*a3*a2 + 16.0*a3*a4*rr + 32.0*a3*ir2 + 48.0*ir4 ) + 24.0*erfc_term*(ir4*ir);
		fh_4th_order = M2A4 * (hBar4  /  (1152.0*(kB*kB*temperature*temperature*reduced_mass*reduced_mass)))  *  (15.0*du*ir3 + 4.0*d3u/r + d4u);
	}
	else
		fh_4th_order = 0.0;

	return fh_2nd_order + fh_4th_order;
}


// fourier space sum 
double System::coulombic_reciprocal() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;

	int    kmax             = ewald_kmax,
	       l[3]             = {0};
	double alpha            = ewald_alpha,
	       k[3]             = {0},
	       k_squared        =  0,
	       position_product =  0,
	       SF_re            =  0,
	       SF_im            =  0, // structure factor 
	       potential        =  0;

	// perform the fourier sum over a hemisphere (skipping certain points to avoid overcounting the face) 
	for(l[0] = 0; l[0] <= kmax; l[0]++) {
		for(l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for(l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				//if norm is out of the sphere, skip
				if( UsefulMath::iidotprod(l,l) > kmax*kmax )
					continue;

				// get the reciprocal lattice vectors 
				for( int p = 0; p < 3; p++ ) {
					k[p] = 0;
					for( int q = 0;   q < 3;   q++ )
						k[p] += 2.0 * pi * pbc.reciprocal_basis[p][q] * l[q];
				}
				k_squared = UsefulMath::dddotprod(k,k);

				// structure factor 
				SF_re = 0;
				SF_im = 0;
				for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
					for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

						if(atom_ptr->frozen) 
							continue; //skip frozen
						if(atom_ptr->charge == 0.0)
							continue; //skip if no charge

						// the inner product of the position vector and the k vector 
						position_product = UsefulMath::dddotprod(k, atom_ptr->pos);

						SF_re += atom_ptr->charge * cos(position_product);
						SF_im += atom_ptr->charge * sin(position_product);

					} // atom
				} // molecule

				potential += exp(-k_squared/(4.0*alpha*alpha))/k_squared*(SF_re*SF_re + SF_im*SF_im);

			} // end for n
		} // end for m
	} // end for l

	potential  *=  4.0 * pi / pbc.volume;

	return potential;
}



double System::coulombic_self() {

	Molecule * molecule_ptr   = nullptr;
	Atom     * atom_ptr       = nullptr;
	double     alpha          = ewald_alpha,
	           self_potential = 0.0;
	
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			if ( atom_ptr->frozen )
				continue;
			atom_ptr->es_self_point_energy = alpha*atom_ptr->charge*atom_ptr->charge/sqrt(pi);
			self_potential -= atom_ptr->es_self_point_energy;
		}
	}

	return self_potential;
}

