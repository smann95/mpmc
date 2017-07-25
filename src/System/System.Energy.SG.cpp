// Silvera-Goldman H2 potential

#include <math.h>

#include "Atom.h"
#include "Molecule.h"
#include "Pair.h"
#include "System.h"

// Silvera-Goldman parameters
//http://www.pnas.org/content/99/3/1129.full.pdf
const double ALPHA                      = 1.713;    // unitless
const double BETA                       = 1.5671;   // 1/a.u.
const double GAMMA                      = 0.00993;  // 1/a.u.^2
const double  C6                        = 12.14;    // multipole term1 a.u.^6
const double  C8                        = 215.2;    // multipole term2 a.u.^8
const double  C10                       = 4813.9;   // multipole term3 a.u.^10
const double  C9                        = 143.1;    // 3-body term a.u.^9
const double  RM                        = 8.321;    // position of max well depth (a.u.) times 1.28


double System::sg() {

	Molecule * molecule_ptr              = nullptr;
	Atom     * atom_ptr                  = nullptr;
	Pair     * pair_ptr                  = nullptr;
	double     rimg                      = 0,
	           r6                        = 0,
	           r8                        = 0,
	           r9                        = 0,
	           r10                       = 0,
	           r_rm                      = 0,
	           repulsive_term            = 0,
	           multipole_term            = 0,
	           exponential_term          = 0,
	           first_r_diff_term         = 0,
	           second_r_diff_term        = 0,
	           first_derivative          = 0,
	           second_derivative         = 0,
	           potential_classical       = 0,
	           potential_fh_second_order = 0,
	           potential                 = 0;
	

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					pair_ptr->rd_energy = 0;
					rimg = pair_ptr->rimg;

					if(rimg < pbc.cutoff) {

						// convert units to Bohr radii 
						rimg /= AU2ANGSTROM;

						// classical pairwise part 
						repulsive_term = exp(ALPHA - BETA*rimg - GAMMA*rimg*rimg);

						r6 = pow(rimg, 6);
						r8 = pow(rimg, 8);
						r9 = pow(rimg, 9);
						r10 = pow(rimg, 10);
						multipole_term = C6/r6 + C8/r8 + C10/r10 - C9/r9;


						r_rm = RM/rimg;
						if(rimg < RM)
							exponential_term = exp(-pow((r_rm - 1.0), 2));
						else
							exponential_term = 1.0;

						potential_classical = (repulsive_term - multipole_term*exponential_term);
						pair_ptr->rd_energy += potential_classical;

						if(feynman_hibbs) {

							// FIRST DERIVATIVE 
							first_derivative = (-BETA -2.0*GAMMA*rimg)*repulsive_term;
							first_derivative += (6.0*C6/pow(rimg, 7) + 8.0*C8/pow(rimg, 9) - 9.0*C9/pow(rimg, 10) + 10.0*C10/pow(rimg, 11))*exponential_term;
							first_r_diff_term = (r_rm*r_rm - r_rm)/rimg;
							first_derivative += -2.0*multipole_term*exponential_term*first_r_diff_term;

							// SECOND DERIVATIVE
							second_derivative = (pow((BETA + 2.0*GAMMA*rimg), 2) - 2.0*GAMMA)*repulsive_term;
							second_derivative += (-exponential_term)*(42.0*C6/pow(rimg, 8) + 72.0*C8/pow(rimg, 10) - 90.0*C9/pow(rimg,11) + 110.0*C10/pow(rimg, 10));
							second_derivative += exponential_term*first_r_diff_term*(12.0*C6/pow(rimg, 7) + 16.0*C8/pow(rimg, 9) - 18.0*C9/pow(rimg,10) + 20.0*C10/pow(rimg, 11));
							second_derivative += exponential_term*pow(first_r_diff_term, 2)*4.0*multipole_term;
							second_r_diff_term = (3.0*r_rm*r_rm - 2.0*r_rm)/(rimg*rimg);
							second_derivative += exponential_term*second_r_diff_term*2.0*multipole_term;

							potential_fh_second_order = pow(METER2ANGSTROM, 2)*(hBar*hBar/(24.0*kB*temperature*(AMU2KG*molecule_ptr->mass)))*(second_derivative + 2.0*first_derivative/rimg);
							pair_ptr->rd_energy += potential_fh_second_order;
						}

						// convert units from Hartrees back to Kelvin 
						pair_ptr->rd_energy *= HARTREE2KELVIN;

					}

				} // recalculate
			} // pair
		} // atom
	} // molecule

	potential = 0;
	for( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next )
		for( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next )
			for( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next )
				potential += pair_ptr->rd_energy;

	return potential;

}

// same as above, but no periodic boundary conditions 
double System::sg_nopbc( Molecule * molecules ) {

	Molecule * molecule_ptr   = nullptr;
	Atom     * atom_ptr       = nullptr;
	Pair     * pair_ptr       = nullptr;
	double     potential      = 0,
	           multipole_term = 0,
	           result         = 0,
	           exp_result     = 0,
	           r=0, r6=0, r8=0, r10=0, r9=0, r_rm=0, r_rm_2=0, r_exp=0;
	

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->recalculate_energy) {

					r = pair_ptr->r/AU2ANGSTROM;


					r6 = pow(r, 6);
					r8 = pow(r, 8);
					r10 = pow(r, 10);
					r9 = pow(r, 9);

					multipole_term = C6/r6 + C8/r8 + C10/r10 - C9/r9;

					if(r < RM) {
						r_rm = RM/r;
						r_rm -= 1.0;
						r_rm_2 = pow(r_rm, 2);
						r_rm_2 *= -1.0;
						r_exp = exp(r_rm_2);

						multipole_term *= r_exp;

					}

					result = ALPHA - BETA*r - GAMMA*r*r;

					exp_result = exp(result);
					exp_result -= multipole_term;
					pair_ptr->rd_energy = HARTREE2KELVIN*exp_result;

				} // recalculate 
			} // pair 
		} // atom 
	} // molecule 

	potential = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
				potential += pair_ptr->rd_energy;

	return potential;

}

