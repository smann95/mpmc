//Copyright 2013-2015 Adam Hogan

#include "Atom.h"
#include "Molecule.h"
#include "Pair.h"
#include "System.h"
#include "UsefulMath.h"



double System::disp_expansion()
{
	double potential = 0.0;

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					// pair LRC
					if( rd_lrc )
						pair_ptr->lrc = disp_expansion_lrc( pair_ptr, pbc.cutoff );

					// make sure we're not excluded or beyond the cutoff
					if(!(pair_ptr->rd_excluded || pair_ptr->frozen)) {
						const double r = pair_ptr->rimg;
						const double r2 = r*r;
						const double r4 = r2*r2;
						const double r6 = r4*r2;
						const double r8 = r6*r2;
						const double r10 = r8*r2;

						double c6 = pair_ptr->c6;
						const double c8 = pair_ptr->c8;
						const double c10 = pair_ptr->c10;

						if( disp_expansion_mbvdw==1)
							c6 = 0.0;

						double repulsion = 0.0;

						if( pair_ptr->epsilon  !=  0.0    &&    pair_ptr->sigma  !=  0.0 )
							repulsion = 315.7750382111558307123944638 * exp(-pair_ptr->epsilon*(r-pair_ptr->sigma)); // K = 10^-3 H ~= 316 K

						if( damp_dispersion )
							pair_ptr->rd_energy = -tt_damping(6,pair_ptr->epsilon*r)*c6/r6-tt_damping(8,pair_ptr->epsilon*r)*c8/r8-tt_damping(10,pair_ptr->epsilon*r)*c10/r10+repulsion;
						else
							pair_ptr->rd_energy = -c6/r6-c8/r8-c10/r10+repulsion;

						if( cavity_autoreject )
						{
							if( r < cavity_autoreject_scale * pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
							if( cavity_autoreject_repulsion != 0.0   &&   repulsion > cavity_autoreject_repulsion )
								pair_ptr->rd_energy = MAXVALUE;
						}
					}

				}
				potential += pair_ptr->rd_energy + pair_ptr->lrc;
			}
		}
	}

	if( disp_expansion_mbvdw == 1 )
	{
		thole_amatrix();
		potential += vdw();
	}

	// calculate self LRC interaction 
	if( rd_lrc )
	{
		for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		{
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			{
				atom_ptr->lrc_self = disp_expansion_lrc_self( atom_ptr,pbc.cutoff );
				potential += atom_ptr->lrc_self;
			}
		}
	}

	return potential;
}



double System::disp_expansion_lrc( Pair * pair_ptr, const double cutoff ) // ignoring the exponential repulsion bit because it decays exponentially
{
	if( !( pair_ptr->frozen ) &&  // disqualify frozen pairs
		((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != pbc.volume) ) { // LRC only changes if the volume change

		pair_ptr->last_volume = pbc.volume;

		return -4.0 * pi * (pair_ptr->c6/(3.0*cutoff*cutoff*cutoff) + pair_ptr->c8/(5.0*cutoff*cutoff*cutoff*cutoff*cutoff)+pair_ptr->c10/(7.0*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff))/pbc.volume;
	}

	else return pair_ptr->lrc; // use stored value
}



double System::tt_damping(int n, double br)
{
	double sum = 0.0;
	
	for( int i=0; i<=n; i++ )
	{
		sum += pow(br,i)/UsefulMath::factorial(i);
	}

	const double result = 1.0 - exp(-br)*sum;

	if( result > 0.000000001 )
		return result;
	else
		return 0.0; // This is so close to zero lets just call it zero to avoid rounding error and the simulation blowing up
}



double System::disp_expansion_lrc_self( Atom * atom_ptr, const double cutoff )
{
	if( !(atom_ptr->frozen) && // disqualify frozen atoms 
 		((atom_ptr->lrc_self == 0.0) || atom_ptr->last_volume != pbc.volume) ) { // LRC only changes if the volume change) 

		atom_ptr->last_volume = pbc.volume;

		if ( extrapolate_disp_coeffs )
		{
			double c10;
			if (atom_ptr->c6!=0.0&&atom_ptr->c8!=0.0)
				c10 = 49.0/40.0*atom_ptr->c8*atom_ptr->c8/atom_ptr->c6;
			else
				c10 = 0.0;

			return -4.0 * pi * (atom_ptr->c6/(3.0*cutoff*cutoff*cutoff)+atom_ptr->c8/(5.0*cutoff*cutoff*cutoff*cutoff*cutoff)+c10/(7.0*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff))/pbc.volume;
		}
		else
			return -4.0 * pi * (atom_ptr->c6/(3.0*cutoff*cutoff*cutoff)+atom_ptr->c8/(5.0*cutoff*cutoff*cutoff*cutoff*cutoff)+atom_ptr->c10/(7.0*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff*cutoff))/pbc.volume;
	}

	return atom_ptr->lrc_self; // use stored value
}
