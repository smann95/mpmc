// Space Research Group
// Department of Chemistry
// University of South Florida

#include "Atom.h"
#include "Molecule.h"
#include "Pair.h"
#include "System.h"


// Lennard-Jones repulsion/dispersion 
double System::lj()
{

	Molecule * molecule_ptr        = nullptr;
	Atom     * atom_ptr            = nullptr;
	Pair     * pair_ptr            = nullptr;
	double     sigma_over_r        =  0,
	           term12              =  0,
	           term6               =  0,
	           sigma_over_r6       =  0,
	           sigma_over_r12      =  0,
	           r                   =  0, 
	           potential           =  0,
	           potential_classical =  0,
	           cutoff              =  0;
	int        i[3]                = {0};
	double     a[3]                = {0};

	//set the cutoff
	if( rd_crystal )
		cutoff = 2.0 * pbc.cutoff * ((double)rd_crystal_order - 0.5);
	else
		cutoff = pbc.cutoff;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;

					// pair LRC 
					if( rd_lrc )
						pair_ptr->lrc = lj_lrc_corr(atom_ptr,pair_ptr,cutoff);

					// to include a contribution, we require
					if(  ( pair_ptr->rimg - SMALL_dR < cutoff     ) &&   // inside cutoff?
					     ( ! pair_ptr->rd_excluded  || rd_crystal ) &&   // either not excluded OR rd_crystal is ON
					     ( ! pair_ptr->frozen )
					) { //not frozen

						//loop over unit cells
						if( rd_crystal ) {
							sigma_over_r6 = 0;
							sigma_over_r12 = 0;
							for ( i[0] = -(rd_crystal_order-1);  i[0] <= rd_crystal_order-1;  i[0]++ )
								for ( i[1] = -(rd_crystal_order-1);  i[1] <= rd_crystal_order-1;  i[1]++ )
									for ( i[2] = -(rd_crystal_order-1);  i[2] <= rd_crystal_order-1;  i[2]++ ) {
										if ( !i[0] && !i[1] && !i[2] && pair_ptr->rd_excluded ) 
											continue; //no i=j=k=0 for excluded pairs (intra-molecular)
										//calculate pair separation (atom with it's image)
										for( int p=0; p<3; p++ ) {
											a[p] = 0;
											for( int q=0; q<3; q++ )
												a[p] += pbc.basis[q][p] * i[q];
											a[p] += atom_ptr->pos[p] - pair_ptr->atom->pos[p];
										}
										r = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

										if ( r > cutoff )
											continue;
										sigma_over_r    = fabs(pair_ptr->sigma)/r;
										sigma_over_r6  += pow(sigma_over_r,6);
										sigma_over_r12 += pow(sigma_over_r,12);
									}
						}
						else { //otherwise, calculate as normal
							sigma_over_r = fabs(pair_ptr->sigma)/pair_ptr->rimg;
							sigma_over_r6 = sigma_over_r*sigma_over_r*sigma_over_r;
							sigma_over_r6 *= sigma_over_r6;
							sigma_over_r12 = sigma_over_r6 * sigma_over_r6;
						}

						// the LJ potential 
						if( spectre ) {
							term6 = 0;
							term12 = sigma_over_r12;
							potential_classical = term12;
						} 
						else {

							if( polarvdw ) 
								term6=0; //vdw calc'd by vdw.c	
							else term6 = sigma_over_r6;


							if(pair_ptr->attractive_only)
								term12 = 0;	
							else
								term12 = sigma_over_r12;


							if( cdvdw_sig_repulsion )
								potential_classical = pair_ptr->sigrep*term12; //C6*sig^6/r^12
							else
								potential_classical = 4.0*pair_ptr->epsilon*(term12 - term6);
						}

						pair_ptr->rd_energy += potential_classical;

						if( feynman_hibbs ) 
							pair_ptr->rd_energy += lj_fh_corr( molecule_ptr,pair_ptr, feynman_hibbs_order, term12, term6 );

						// if cavity_autoreject is on (cavity_autoreject_absolute is performed in energy.c)
						if( cavity_autoreject )
							if(pair_ptr->rimg < cavity_autoreject_scale*fabs(pair_ptr->sigma))
								pair_ptr->rd_energy = MAXVALUE;

					} //count contributions

				} // if recalculate

				// sum all of the pairwise terms 
				potential += pair_ptr->rd_energy + pair_ptr->lrc;

			} // pair
		} // atom
	} // molecule

	// molecule self-energy for rd_crystal -> energy of molecule interacting with its periodic neighbors

	if( rd_crystal )
		for( molecule_ptr = molecules;   molecule_ptr;   molecule_ptr = molecule_ptr->next )
			for( atom_ptr = molecule_ptr->atoms;   atom_ptr;   atom_ptr = atom_ptr->next)
				potential += rd_crystal_self( atom_ptr, cutoff );

	// calculate self LRC interaction
	if( rd_lrc ) 
		for( molecule_ptr = molecules;   molecule_ptr;   molecule_ptr = molecule_ptr->next )
			for( atom_ptr = molecule_ptr->atoms;   atom_ptr;   atom_ptr = atom_ptr->next)
				potential += lj_lrc_self( atom_ptr, cutoff );

	return potential;

}



double System::lj_lrc_corr( Atom * atom_ptr,  Pair * pair_ptr, double cutoff ) 
{

	double sig_cut  = 0,
	       sig3     = 0,
	       sig_cut3 = 0,
	       sig_cut9 = 0;

	// include the long-range correction.    I'm  not sure that I'm handling spectre pairs correctly
	// we can't use rd_excluded flag, since that disqualifies inter-molecular, but that DOES contribute to LRC
	// ALL OF THESE MUST BE TRUE TO PERFORM LRC CALCULATION
	if(  ( pair_ptr->epsilon != 0 && pair_ptr->sigma != 0 ) &&  //if these are zero, then we won't waste our time
	    !( atom_ptr->spectre && pair_ptr->atom->spectre   ) &&  //i think we want to disqualify s-s pairs 
	    !( pair_ptr->frozen                               ) &&  //disqualify frozen pairs
	    ((pair_ptr->lrc == 0.0) || pair_ptr->last_volume != pbc.volume)
	) { //LRC only changes if the volume change

		pair_ptr->last_volume = pbc.volume;
		sig_cut = fabs(pair_ptr->sigma)/cutoff;
		sig3 = fabs(pair_ptr->sigma);
		sig3 *= sig3*sig3;
		sig_cut3 = sig_cut*sig_cut*sig_cut;
		sig_cut9 = sig_cut3*sig_cut3*sig_cut3;

		if( cdvdw_sig_repulsion )
			return (4.0/9.0) * pi * pair_ptr->sigrep * sig3 * sig_cut9   /   pbc.volume;
		else if ( polarvdw ) //only repulsion term, if polarvdw is on
			return (16.0/9.0) * pi * pair_ptr->epsilon * sig3 * sig_cut9   /   pbc.volume;
		else //if polarvdw is off, do the usual thing
			return ((16.0/3.0)*pi * pair_ptr->epsilon * sig3)  *  ((1.0/3.0)*sig_cut9 - sig_cut3)  /  pbc.volume;
	}
	else return pair_ptr->lrc; //use stored value

}


double System::lj_lrc_self( Atom * atom_ptr, double cutoff ) 
{
	double sig_cut, sig3, sig_cut3, sig_cut9;

	if (  ((atom_ptr->sigma != 0)  && (atom_ptr->epsilon != 0)) &&  // non-zero parameters
	     ! (atom_ptr->frozen                                  ) &&  // not frozen
	     ! (atom_ptr->spectre)                                      // not spectre 
	) { 
		
		sig_cut  = fabs(atom_ptr->sigma)/cutoff;
		sig3     = fabs(atom_ptr->sigma);
		sig3    *= sig3*sig3;
		sig_cut3 = sig_cut*sig_cut*sig_cut;
		sig_cut9 = sig_cut3*sig_cut3*sig_cut3;

		if( cdvdw_sig_repulsion ) 
			return (1.0/3.0) * pi * hBar/kB * au2invseconds * atom_ptr->omega * atom_ptr->polarizability * atom_ptr->polarizability/sig3 * sig_cut9   /   pbc.volume;
		else if( polarvdw ) //only repulsion term, if polarvdw is on
			return (16.0/9.0) * pi * atom_ptr->epsilon * sig3 * sig_cut9   /   pbc.volume;
		else //if polarvdw is off, do the usual thing
			return ((16.0/3.0) * pi * atom_ptr->epsilon*sig3)  *  ((1.0/3.0) * sig_cut9 - sig_cut3) / pbc.volume;
	}

	return 0;
}



double System::lj_fh_corr( Molecule * molecule_ptr, Pair * pair_ptr, int order, double term12, double term6 )
{
	double reduced_mass;
	double dE, d2E, d3E, d4E; //energy derivatives
	double corr;
	double ir = 1.0/pair_ptr->rimg;
	double ir2 = ir*ir;
	double ir3 = ir2*ir;
	double ir4 = ir3*ir;

	if ( (order != 2) && (order != 4) ) 
		throw invalid_setting; //must be order 2 or 4

	reduced_mass = AMU2KG*molecule_ptr->mass*pair_ptr->molecule->mass /
	               (molecule_ptr->mass+pair_ptr->molecule->mass);

	if( cdvdw_sig_repulsion ) {
		dE = -6.0*pair_ptr->sigrep*(2.0*term12 - term6) * ir;
		d2E = 6.0*pair_ptr->sigrep*(26.0*term12 - 7.0*term6) * ir2;
	} else {
		dE = -24.0*pair_ptr->epsilon*(2.0*term12 - term6) * ir;
		d2E = 24.0*pair_ptr->epsilon*(26.0*term12 - 7.0*term6) * ir2;
	}

	//2nd order correction
	corr = M2A2 *
		(hBar2 / (24.0*kB*temperature*reduced_mass)) *
		(d2E + 2.0*dE/pair_ptr->rimg);

	if(order >= 4) {

		if( cdvdw_sig_repulsion ) {
			d3E = -336.0*pair_ptr->sigrep*(6.0*term12 - term6) * ir3;
			d4E = 3024.0*pair_ptr->sigrep*(10.0*term12 - term6) * ir4;
		} else { 
			d3E = -1344.0*pair_ptr->epsilon*(6.0*term12 - term6) * ir3;
			d4E = 12096.0*pair_ptr->epsilon*(10.0*term12 - term6) * ir4;
		}
	
		//4th order corection
		corr += M2A4 *
			(hBar4/(1152.0*kB2*temperature*temperature*reduced_mass*reduced_mass)) *
			( 15.0*dE*ir3 +	4.0*d3E*ir + d4E );
	}

	return corr;
}



double System::rd_crystal_self( Atom * aptr, double cutoff ) 
{

	double curr_pot, term12, term6;
	double sigma_over_r6, sigma_over_r12, sigma_over_r, r;
	int i[3], p, q;
	double a[3];
	curr_pot = 0;

	if ( aptr->sigma == 0 && aptr->epsilon == 0 ) return 0; //skip if no LJ interaction

	sigma_over_r6 = 0;
	sigma_over_r12 = 0; //need to init these guys

	for ( i[0] = -(rd_crystal_order-1); i[0]<=rd_crystal_order-1; i[0]++ )
		for ( i[1] = -(rd_crystal_order-1); i[1]<=rd_crystal_order-1; i[1]++ )
			for ( i[2] = -(rd_crystal_order-1); i[2]<=rd_crystal_order-1; i[2]++ ) {
				if ( !i[0] && !i[1] && !i[2] ) continue; //no (0,0,0)
				//calculate pair separation (atom with it's image)
				for ( p=0; p<3; p++ ) {
					a[p] = 0;
					for ( q=0; q<3; q++ )
						a[p] += pbc.basis[q][p] * i[q];
				}
				r = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

				if ( r > cutoff ) 
					continue; //too far away! will be included in LRC if enabled
				sigma_over_r = fabs(aptr->sigma)/r;
				sigma_over_r6 += 0.5*pow(sigma_over_r,6); //multiply by 0.5 to get counting correct
				sigma_over_r12 += 0.5*pow(sigma_over_r,12);
			}

	if( spectre ) {
		term6 = 0;
		curr_pot = term12 = sigma_over_r12;
	} 
	else {
		if( polarvdw )
			term6=0; //vdw calc'd by vdw.c
		else
			term6 = sigma_over_r6;

		if(aptr->sigma < 0.0)
			term12 = 0; //attractive only
		else
			term12 = sigma_over_r12;

		if( cdvdw_sig_repulsion )
			curr_pot = 0.75*hBar/kB*au2invseconds*aptr->omega*aptr->polarizability*aptr->polarizability / pow(aptr->sigma,6)*term12; //C6*sig^6/r^12
		else if( polarvdw )
			curr_pot = 4.0*aptr->epsilon*term12;
		else
			curr_pot = 4.0*aptr->epsilon*(term12 - term6);
	}
	return curr_pot;
}



double System::lj_buffered_14_7()
{
	double potential = 0.0, potential_classical;
	double r_over_sigma, first_term, second_term;

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	Pair     * pair_ptr     = nullptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {
					pair_ptr->rd_energy = 0;

					// make sure we're not excluded or beyond the cutoff
					if(!((pair_ptr->rimg > pbc.cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
						r_over_sigma = pair_ptr->rimg/pair_ptr->sigma;
						first_term = pow(1.07/(r_over_sigma+0.07),7);
						second_term = (1.12/(pow(r_over_sigma,7)+0.12)-2);
						potential_classical = pair_ptr->epsilon*first_term*second_term;
						pair_ptr->rd_energy += potential_classical;

						// cavity autoreject
						if(cavity_autoreject)
							if(pair_ptr->rimg < cavity_autoreject_scale * pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
					}
				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}

double System::lj_buffered_14_7_nopbc()
{
	double potential           = 0,
	       potential_classical = 0,
	       r_over_sigma        = 0,
	       first_term          = 0,
	       second_term         = 0;

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	Pair     * pair_ptr     = nullptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {
					pair_ptr->rd_energy = 0;

					// make sure we're not excluded or beyond the cutoff
					if(!((pair_ptr->rimg > pbc.cutoff) || pair_ptr->rd_excluded || pair_ptr->frozen)) {
						r_over_sigma = pair_ptr->rimg/pair_ptr->sigma;
						first_term = pow(1.07/(r_over_sigma+0.07),7);
						second_term = (1.12/(pow(r_over_sigma,7)+0.12)-2);
						potential_classical = pair_ptr->epsilon*first_term*second_term;
						pair_ptr->rd_energy += potential_classical;

						// cavity autoreject
						if( cavity_autoreject )
							if(pair_ptr->rimg < cavity_autoreject_scale*pair_ptr->sigma)
								pair_ptr->rd_energy = MAXVALUE;
					}
				}
				potential += pair_ptr->rd_energy;
			}
		}
	}

	return potential;
}