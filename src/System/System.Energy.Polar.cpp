#include "constants.h"
#include "Output.h"
#include "Pair.h"
#include "SafeOps.h"
#include "System.h"
#include "UsefulMath.h"



void print_matrix( int N, double **matrix ) {

	int i, j;

	printf("\n");
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			printf("%.3f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


void zero_out( Molecule * m ) {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	int p;

	//zero out the electric field for each site
	for ( mptr = m; mptr; mptr=mptr->next ) {
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for ( p=0; p<3; p++ ) {
				aptr->ef_static[p] = 0.0;
				aptr->ef_static_self[p] = 0.0;
			}
		}
	}

	return;
}


// get the induction energy 
double System::polar() {

	Molecule * molecule_ptr     = nullptr;
	Atom     * atom_ptr         = nullptr;

	int        num_iterations   =  0;
	double     potential        =  0;
	char       linebuf[maxLine] = {0};

	// take measures to let N fluctuate
	if(  (ensemble == ENSEMBLE_UVT || ensemble == ENSEMBLE_REPLAY)    &&    ! polar_zodid )
		thole_resize_matrices();

	// get the A matrix
	if(  ! polar_zodid  ) {
		thole_amatrix();
		if( polarizability_tensor ) {
			Output::out( "POLAR: A matrix:\n" );
			print_matrix(   3 * ((int)checkpoint->thole_N_atom), A_matrix   );
		}
	}

	// find the dipoles

	if( polar_ewald_full ) {
		//do a full-ewald polarization treatment
		ewald_full(); 

	} else if( polar_iterative ) {	
		//solve the self-consistent problem
		thole_field(); //calc e-field
		num_iterations = thole_iterative(); //calc dipoles

		nodestats->polarization_iterations = (double) num_iterations; //statistics
		observables->dipole_rrms = get_dipole_rrms();

		if( iterator_failed ) {
			switch( ensemble ) {
			case ENSEMBLE_UVT:
			case ENSEMBLE_NVT:
			case ENSEMBLE_NVE:
			case ENSEMBLE_NPT:
				sprintf(linebuf,"POLAR: polarization iterative solver convergence failure on mc step %d.\n", step);
				Output::err(linebuf);
			break;
			case ENSEMBLE_REPLAY:
				sprintf(linebuf,"POLAR: polarization iterative solver convergence failure on configuration %d.\n", step);
				Output::err(linebuf);
			break;
			case ENSEMBLE_SURF:
			case ENSEMBLE_SURF_FIT:
			case ENSEMBLE_TE:
			default:
				sprintf(linebuf,"POLAR: polarization iterative solver failed to reach convergence.\n");
				Output::err(linebuf);
			}
		}

	} else {	
		//do matrix inversion
		thole_field(); //calc e-field
		thole_bmatrix(); //matrix inversion
		thole_bmatrix_dipoles(); //get dipoles

		// output the 3x3 molecular polarizability tensor 
		if( polarizability_tensor ) {
			Output::out("POLAR: B matrix:\n");
			print_matrix(3*((int)checkpoint->thole_N_atom), B_matrix);
			thole_polarizability_tensor();
			throw exception_ok;
		}
	}

	// calculate the polarization energy as 1/2 mu*E
	potential = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			potential += UsefulMath::dddotprod(atom_ptr->mu,atom_ptr->ef_static);
			if( polar_palmo )
				potential += UsefulMath::dddotprod(atom_ptr->mu,atom_ptr->ef_induced_change);
		}
	}
	potential *= -0.5;

	#ifdef DEBUG
		fprintf(stderr,"mu MOLECULE ATOM * DIPOLES * STATIC * INDUCED * pot/atom -0.5*mu*E_s\n");
		for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
			for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
				fprintf(stderr,"mu %4d %4d * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %8.5lf %8.5lf %8.5lf * %lf %lf\n", 
					molecule_ptr->id, atom_ptr->id, 
					atom_ptr->mu[0], atom_ptr->mu[1], atom_ptr->mu[2],
					atom_ptr->ef_static[0], atom_ptr->ef_static[1], atom_ptr->ef_static[2],
					atom_ptr->ef_induced[0], atom_ptr->ef_induced[1], atom_ptr->ef_induced[2],
					potential/system->natoms, -0.5*atom_ptr->mu[0]*atom_ptr->ef_static[0]);
			}
		}
	#endif

	return potential;
}


// RRMS of dipoles
double System::get_dipole_rrms() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;

	double     N            = 0,
	           dipole_rrms  = 0;

	

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			if(std::isfinite(atom_ptr->dipole_rrms)) 
				dipole_rrms += atom_ptr->dipole_rrms;
			N++;
		}
	}
	return dipole_rrms / N;
}


// calculate the dipole field tensor 
void System::thole_amatrix() {

	int     ii, jj,
	        N = natoms;
	Pair  * pair_ptr = nullptr;
	double  damp1=0, damp2=0, wdamp1=0, wdamp2=0, v=0, s=0,
	        r=0, r2=0, ir3=0, ir5=0, ir=0,

	        rcut     = pbc.cutoff,
	        rcut2    = rcut  * rcut,
	        rcut3    = rcut2 * rcut,
	
	        l        = polar_damp,
	        l2       = l  * l,
	        l3       = l2 * l,

	        explr    = 0, //exp(-l*r)
	        explrcut = exp(-l*rcut);

	zero_out_amatrix(N);

	// set the diagonal blocks 
	for( int i = 0; i < N; i++) {
		ii = i*3;
		for( int p = 0; p < 3; p++) {
			if(atom_array[i]->polarizability != 0.0)
				A_matrix[ii+p][ii+p] = 1.0/atom_array[i]->polarizability;
			else
				A_matrix[ii+p][ii+p] = MAXVALUE;
		}
	}

	// calculate each Tij tensor component for each dipole pair
	for( int i = 0; i < (N - 1); i++) {
		ii = i*3;
		pair_ptr = atom_array[i]->pairs;
		for( int j = (i + 1); j < N; j++, pair_ptr = pair_ptr->next) {
			jj = j*3;

			r = pair_ptr->rimg;
			r2 = r*r;

			// inverse displacements
			if(pair_ptr->rimg == 0.)
				ir3 = ir5 = MAXVALUE;
			else {
				ir = 1.0/r;
				ir3 = ir*ir*ir;
				ir5 = ir3*ir*ir;
			}

			//evaluate damping factors
			switch( damp_type ) {
				case DAMPING_OFF:
					if ( pair_ptr->es_excluded )
						damp1 = damp2 = wdamp1 = wdamp2 = 0.0; 
					else 
						damp1 = damp2 = wdamp1 = wdamp2 = 1.0; 
					break;
				case DAMPING_LINEAR:
					s = l * pow(atom_array[i]->polarizability*atom_array[j]->polarizability, 1.0/6.0);
					v = r/s;
					if ( r < s ) {
						damp1 = (4.0 - 3.0*v)*v*v*v;
						damp2 = v*v*v*v;
					} else {
						damp1 = damp2 = 1.0;
					}
					break;
				case DAMPING_EXPONENTIAL:
					explr = exp(-l*r);
					damp1 = 1.0 - explr*(0.5*l2*r2 + l*r + 1.0);
					damp2 = damp1 - explr*(l3*r2*r/6.0);
					if( polar_wolf_full ) { //subtract off damped interaction at r_cutoff
						wdamp1 = 1.0 - explrcut*(0.5*l2*rcut2+l*rcut+1.0);
						wdamp2 = wdamp1 - explrcut*(l3*rcut3/6.0);
					}
					break;
				default:
					Output::err("error: something unexpected happened in thole_matrix.c");
			}

			// build the tensor
			for( int p = 0; p < 3; p++) {
				for( int q = 0; q < 3; q++) {

					A_matrix[ii+p][jj+q] = -3.0*pair_ptr->dimg[p]*pair_ptr->dimg[q]*damp2*ir5;
					if( polar_wolf_full ) 
						A_matrix[ii+p][jj+q] -= -3.0*pair_ptr->dimg[p]*pair_ptr->dimg[q]*wdamp2*ir*ir/rcut3;

					// additional diagonal term
					if(p == q) {
						A_matrix[ii+p][jj+q] += damp1*ir3;
						if( polar_wolf_full ) 
							A_matrix[ii+p][jj+q] -= wdamp1/(rcut3);
					}	
				}
			}

			// set the lower half of the tensor component 
			for( int p = 0; p < 3; p++)
				for( int q = 0; q < 3; q++)
					A_matrix[jj+p][ii+q] = A_matrix[ii+p][jj+q];

		} // end j 
	} // end i 

	return;
}


void System::zero_out_amatrix ( int N ) {
	
	// zero out the matrix 
	for(int i = 0; i < 3*N; i++)
		for(int j = 0; j < 3*N; j++)
			A_matrix[i][j] = 0;
	return;
}


//do full polarization calculation using ewald
//see nymand and linse jcp 112 6152 (2000)
void System::ewald_full() {

	const int MAX_ITERATION_COUNT = 128;

	//int max_iter=system->polar_max_iter;  (unused variable)
	int keep_iterating, iteration_counter;

	//calculate static e-field
	zero_out(molecules);
	recip_term();
	real_term();

	//calculate induced e-field
	init_dipoles_ewald();

	keep_iterating = 1;
	iteration_counter = 0;
	while ( keep_iterating ) {

		if ( iteration_counter >= MAX_ITERATION_COUNT && polar_precision ) {
			iterator_failed = 1;
			return;
		}

		//set induced field to zero
		clear_ef_induced();

		//calculate induced field
		induced_real_term();
		induced_recip_term();
		induced_corr_term();

		if ( polar_rrms || polar_precision > 0 )
			calc_dipole_rrms();

		//recalculate dipoles using new induced field
		new_dipoles(iteration_counter);

		keep_iterating = are_we_done_yet(iteration_counter);

		if ( polar_palmo &&	 !keep_iterating ) //if last iteration
			ewald_palmo_contraction();

		iteration_counter++;
	}

	return;
}


//we deviate from drexel's treatment, and instead do a trig identity to get from a pairwise sum to two atomwise rums
//or ignore drexel, and derive this term from eq (29) in nymand and linse
void System::recip_term() {

	Molecule * mptr       = nullptr;
	Atom     * aptr       = nullptr;
	int        l[3]       = {0},
	           kmax       =  ewald_kmax;
	double     ea         =  polar_ewald_alpha, //actually sqrt(ea)
	           k[3]       = {0}, 
	           k2         =  0,
	           kweight[3] = {0},
	           float1     =  0,
	           float2     =  0;


	//k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
	for (l[0] = 0; l[0] <= kmax; l[0]++) {
		for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
			for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {

				// if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
				if ( UsefulMath::iidotprod(l,l) > kmax*kmax ) continue; 

				for( int p = 0; p < 3; p++) {
					k[p] = 0;
					for( int q = 0; q < 3; q++)
						k[p] += 2.0 * pi * pbc.reciprocal_basis[p][q]*l[q];
				}
				k2 = UsefulMath::dddotprod(k,k);

				kweight[0] = k[0]/k2 * exp(-k2/(4.0*ea*ea));
				kweight[1] = k[1]/k2 * exp(-k2/(4.0*ea*ea));
				kweight[2] = k[2]/k2 * exp(-k2/(4.0*ea*ea));

				float1 = float2 = 0;
				for (mptr = molecules; mptr; mptr=mptr->next )
					for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
						float1 += aptr->charge * cos( UsefulMath::dddotprod(k,aptr->pos) );
						float2 += aptr->charge * sin( UsefulMath::dddotprod(k,aptr->pos) );
					}

				for (mptr = molecules; mptr; mptr=mptr->next )
					for (aptr = mptr->atoms; aptr; aptr=aptr->next ) {
						for ( int p=0; p<3; p++ ) {
							aptr->ef_static[p] += kweight[p] * sin( UsefulMath::dddotprod(k,aptr->pos) ) * float1;
							aptr->ef_static[p] -= kweight[p] * cos( UsefulMath::dddotprod(k,aptr->pos) ) * float2;
						}
					}

			} //l2
		} //l1
	} //l0

	for ( mptr=molecules; mptr; mptr=mptr->next ) {
		for ( aptr=mptr->atoms; aptr; aptr=aptr->next ) {
			for( int p=0; p<3; p++ ) {
				//factor of 2 more, since we only summed over hemisphere
				aptr->ef_static[p] *= 8.0*pi / pbc.volume;
			}
		}
	}

	return;
}


void System::real_term() {

	Molecule * mptr   = nullptr;
	Atom     * aptr   = nullptr;
	Pair     * pptr   = nullptr;
	
	double     r      = 0,
	           r2     = 0,
	           factor = 0,
	           a      = polar_ewald_alpha; //some ambiguity between ea and ea^2 across the literature;
	

	for( mptr = molecules; mptr; mptr=mptr->next ) {
		for( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for( pptr = aptr->pairs; pptr; pptr=pptr->next ) { //for each pair
				if (pptr->frozen) continue; //if the pair is frozen (i.e. MOF-MOF interaction) it doesn't contribute to polar
				r = pptr->rimg;
				if ( (r > pbc.cutoff) || (r == 0.0) ) continue; //if outside cutoff sphere (not sure why r==0 ever) -> skip
				r2 = r*r; 
				if (pptr->es_excluded) {
					//need to subtract self-term (interaction between a site and a neighbor's screening charge (on the same molecule)
					factor = (2.0 * a * OneOverSqrtPi * exp(-a*a*r2) * r - erf(a*r))    /    (r*r2);
					for( int p=0; p<3; p++ ) {
						aptr->ef_static[p] += factor*pptr->atom->charge * pptr->dimg[p];
						pptr->atom->ef_static[p] -= factor*aptr->charge * pptr->dimg[p];
					}
				} //excluded
				else { //not excluded

					factor = (2.0 * a * OneOverSqrtPi * exp(-a*a*r2) * r   +   erfc(a*r))    /    (r2*r);
					for ( int p=0; p<3; p++ ) { // for each dim, add e-field contribution for the pair
						aptr->ef_static[p] += factor*pptr->atom->charge * pptr->dimg[p];
						pptr->atom->ef_static[p] -= factor*aptr->charge * pptr->dimg[p];
					}
				} //excluded else
			} //ptr
		} //aptr
	} //mptr

	return;
}


//set zeroth iteration dipoles
void System::init_dipoles_ewald() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;

	for ( mptr=molecules; mptr; mptr=mptr->next )
		for ( aptr=mptr->atoms; aptr; aptr=aptr->next )
			for( int p=0; p<3; p++ ) {
				aptr->old_mu[p] = 0;
				aptr->new_mu[p] = aptr->mu[p] = aptr->polarizability*aptr->ef_static[p];
			}

	return;
}


//reset the ef_induced values
void System::clear_ef_induced() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	

	for( mptr=molecules; mptr; mptr=mptr->next )
		for( aptr=mptr->atoms; aptr; aptr=aptr->next )
			for( int p=0; p<3; p++ )
				aptr->ef_induced[p] = 0;

	return;
}


void System::induced_recip_term() {
	Molecule * mptr     = nullptr;
	Atom     * aptr     = nullptr,
	        ** aarray   = nullptr;
	int        N        =  0,
	           l[3]     = {0},
	           kmax     = ewald_kmax;
	double     Psin     =  0,
	           Pcos     =  0,
	           kweight  =  0,
	           a        = polar_ewald_alpha,
	           dotprod1 =  0,
	           dotprod2 =  0,
	           k[3]     = {0},
	           k2       =  0;

	//make atom array
	for ( mptr=molecules; mptr; mptr=mptr->next ) {
		for ( aptr=mptr->atoms; aptr; aptr=aptr->next ) {
			SafeOps::realloc(aarray, sizeof(Atom *) * (N+1), __LINE__, __FILE__);
			aarray[N] = aptr;
			N++;
		}
	}

	//k-space sum (symmetry for k -> -k, so we sum over hemisphere, avoiding double-counting on the face)
	for ( l[0] = 0; l[0] <= kmax; l[0]++ ) 
		for ( l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++ ) 
			for ( l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++ ) {

				// if |l|^2 > kmax^2, then it doesn't contribute (outside the k-space cutoff)
				if ( UsefulMath::iidotprod(l,l) > kmax*kmax ) continue; 

				for( int p = 0; p < 3; p++) {
					k[p] = 0;
					for( int q = 0; q < 3; q++)
						k[p] += 2.0 * pi * pbc.reciprocal_basis[p][q]*l[q];
				}
				k2 = UsefulMath::dddotprod(k,k);

				for( int p=0; p<3; p++ ) 	
					kweight = 8.0 * pi/pbc.volume * exp(-k2/(4.0*a*a))/k2 * k[p];

				//calculate Pcos, Psin for this k-point
				Pcos = Psin = 0;
				for( int j=0; j<N; j++ ) {
					dotprod1 = UsefulMath::dddotprod(k,aarray[j]->mu);
					dotprod2 = UsefulMath::dddotprod(k,aarray[j]->pos);	
					Pcos += dotprod1 * cos(dotprod2);
					Psin += dotprod1 * sin(dotprod2);
				}

				//calculate ef_induced over atom array
				for( int i=0; i<N; i++ ) {
					//for each cartesian dimension
					for( int p=0; p<3; p++ ) {
				
						dotprod1 = UsefulMath::dddotprod(k,aarray[i]->pos);
						aarray[i]->ef_induced[p] += kweight * ( -sin(dotprod1)*Psin - cos(dotprod1)*Pcos );
	
					} //dim
				} //ef_incuded over atom array
			} //kspace	
		
	free(aarray);

	return;
}


void System::induced_real_term() {

	Molecule * mptr    = nullptr;
	Atom     * aptr    = nullptr;
	Pair     * pptr    = nullptr;
	double     erfcar  = 0, 
	           expa2r2 = 0, 
	           T       = 0, //dipole-interaction tensor component
	           s1      = 0, 
	           s2      = 0, //common term (s_2 via eq 10. JCP 133 243101)
	           a       = polar_ewald_alpha, //ewald damping
	           l       = polar_damp, //polar damping
	           r=0, ir=0, ir3=0, ir5=0;
	
	

	for ( mptr = molecules; mptr; mptr=mptr->next ) {
		for ( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for ( pptr = aptr->pairs; pptr; pptr=pptr->next ) {
				if ( aptr->polarizability == 0 || pptr->atom->polarizability == 0 ) 
					continue; //don't waste CPU time
				if ( pptr->rimg > pbc.cutoff )
					continue; //if outside cutoff sphere skip
				//some things we'll need
				r=pptr->rimg;
				ir=1.0/r; ir3=ir*ir*ir; ir5=ir*ir*ir3;
				erfcar=erfc(a*r);
				expa2r2=exp(-a*a*r*r);

				//E_static_realspace_i = sum(i!=j) d_xi d_xj erfc(a*r)/r u_j 
				s2 = erfcar + 2.0*a*r*OneOverSqrtPi * expa2r2 + 4.0*a*a*a*r*r*r/3.0*OneOverSqrtPi*expa2r2 - damp_factor(l*r,3) ;

				for( int p=0; p<3; p++ ) {
					for( int q=p; q<3; q++ ) {  //it's symmetric!

						if ( p==q )
							s1 = erfcar + 2.0*a*r*OneOverSqrtPi * expa2r2 - damp_factor(l*r,2);
						else
							s1 = 0;

						//real-space dipole interaction tensor
						T = 3.0 * pptr->dimg[p] * pptr->dimg[q] * s2 * ir5 - s1 * ir3; 

						aptr->ef_induced[p] += T*pptr->atom->mu[q];
						pptr->atom->ef_induced[p] += T*aptr->mu[q];

						if ( p!=q ) {
							aptr->ef_induced[q] += T*pptr->atom->mu[p];
							pptr->atom->ef_induced[q] += T*aptr->mu[p];
						}

					} //loop over q dim
				} //loop over p dim
			} //pptr loop
		} //aptr loop
	} //mptr loop
						
	return;
}


//damping term (see e.g. Souaille et al., Comp. Phys. Comm. 180 276-301) below eq (9).
//signs are intentionally reversed (different convention)
double System::damp_factor( double t, int i ) {

	double temp = 1.0 + t + 0.5 * t*t;

	if( i == 3 )
		temp += t*t*t / 6.0;

	return temp * exp(-t);
}


void System::induced_corr_term() {
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	
	double     a          = polar_ewald_alpha,
	           totalmu[3] = {0};

	for( int p=0; p<3; p++ ) 
		totalmu[p] = 0;

	for( mptr = molecules; mptr; mptr=mptr->next ) 
		for( aptr = mptr->atoms; aptr; aptr=aptr->next )
			for( int p=0; p<3; p++ )
				totalmu[p] += aptr->mu[p];

	//other term
	for( mptr = molecules; mptr; mptr=mptr->next ) 
		for( aptr = mptr->atoms; aptr; aptr=aptr->next )
			for( int p=0; p<3; p++ ) 

				aptr->ef_induced[p] += -4.0 * pi/(3.0*pbc.volume) * totalmu[p]   +   4.0*a*a*a/(3.0*SqrtPi) * aptr->mu[p];

	return;
}


void System::calc_dipole_rrms() {
	
	Atom ** aa    = atom_array;
	double  carry = 0;
	
	// get the dipole RRMS 
	for( int i = 0; i < natoms; i++) {
		// mean square difference 
		aa[i]->dipole_rrms = 0;
		for( int p = 0; p < 3; p++) {
			carry = aa[i]->new_mu[p] - aa[i]->old_mu[p];
			aa[i]->dipole_rrms += carry*carry;
		}

		// normalize
		aa[i]->dipole_rrms /= UsefulMath::dddotprod(aa[i]->new_mu,aa[i]->new_mu);
		aa[i]->dipole_rrms = sqrt(aa[i]->dipole_rrms);
		if( !std::isfinite(aa[i]->dipole_rrms) )
			aa[i]->dipole_rrms = 0;
	}

	#ifdef DEBUG
		double totalrms = 0;
		for ( int i=0; i< system->natoms; i++ ) {
			totalrms +=  aa[i]->dipole_rrms;
		}
		fprintf(stderr,"TOTAL DIPOLE RRMS %lf\n", totalrms);
	#endif

	return;
}


void System::new_dipoles(int count) {

	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	
	for( mptr = molecules; mptr; mptr=mptr->next ) {
		for( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			for( int p=0; p<3; p++ ) {

				//set dipoles
				aptr->old_mu[p] = aptr->mu[p];
				if( polar_sor ) {
					aptr->new_mu[p] = aptr->polarizability*(aptr->ef_static[p]+aptr->ef_induced[p]);
					aptr->new_mu[p] = aptr->mu[p] = polar_gamma*aptr->new_mu[p] + (1.0-polar_gamma)*aptr->old_mu[p];
				}
				else if( polar_esor ) {
					aptr->new_mu[p] = aptr->polarizability*(aptr->ef_static[p]+aptr->ef_induced[p]);
					aptr->new_mu[p] = aptr->mu[p] = (1.0 - exp(-polar_gamma*(count+1)))*aptr->new_mu[p] +
						exp(-polar_gamma*(count+1))*aptr->old_mu[p];
				}
				else {
					//if no sor, still need new_mu for polar_palmo
					aptr->mu[p] = aptr->new_mu[p] = aptr->polarizability*(aptr->ef_static[p]+aptr->ef_induced[p]);
				}
			
			}
		}
	}

	return;
}


int System::are_we_done_yet( int iteration_counter ) {
	
	Atom ** aa            = atom_array;
	int     N             = natoms;
	double  allowed_sqerr = 0,
	        error         = 0;
	
	if( polar_precision == 0.0) {	// ... by fixed iteration ...
		if(iteration_counter != polar_max_iter)
			return 1;
	} 

	else { // ... or by dipole precision 
		allowed_sqerr = polar_precision * polar_precision * DEBYE2SKA * DEBYE2SKA;
		for( int i = 0; i < N; i++) { //check the change in each dipole component
			for( int p = 0; p < 3; p++) {
				error = aa[i]->new_mu[p] - aa[i]->old_mu[p];
				if(error*error > allowed_sqerr)
					return 1; //we broke tolerance
			}
		}
	}

	return 0;
}


void System::ewald_palmo_contraction() {
	
	Molecule * mptr = nullptr;
	Atom     * aptr = nullptr;
	
	//set induced field to zero
	clear_ef_induced();

	//calculate induced field
	induced_real_term();
	induced_recip_term();
	induced_corr_term();

	for( mptr = molecules; mptr; mptr=mptr->next ) {
		for( aptr = mptr->atoms; aptr; aptr=aptr->next ) {
			if( aptr->polarizability == 0 ) continue;
			for( int p=0; p<3; p++ ) {
				aptr->ef_induced_change[p] = //current induced - last induced (backed out from dipole values)
					aptr->ef_induced[p] - (aptr->new_mu[p]/aptr->polarizability - aptr->ef_static[p]);
			}
		}
	}

	return;
}


//called from energy/polar.c
// calculate the field with periodic boundaries
void System::thole_field() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;

	// zero the field vectors 
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			for( int p = 0; p < 3; p++) {
				atom_ptr->ef_static[p] = 0;
				atom_ptr->ef_static_self[p] = 0;
			}

		}
	}

	// calculate the electrostatic field 
	if( polar_ewald )
		ewald_estatic();
	else if( polar_wolf || polar_wolf_full )
		thole_field_wolf();
	else
		thole_field_nopbc();

}


// calculate the field without ewald summation/wolf 
void System::thole_field_nopbc() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	Pair     * pair_ptr     = nullptr;
	double     r            = 0;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(pair_ptr->frozen) 
					continue;
				if (molecule_ptr == pair_ptr->molecule) 
					continue; //don't let molecules polarize themselves
				
				r = pair_ptr->rimg;

				//inclusive near the cutoff
				if((r - SMALL_dR < pbc.cutoff ) && ( r != 0. )) {

					for( int p = 0; p < 3; p++ ) {
						atom_ptr->ef_static[p] += pair_ptr->atom->charge*pair_ptr->dimg[p]/(r*r*r);
						pair_ptr->atom->ef_static[p] -= atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
					}

				} // cutoff 

			} // pair
		} // atom
	} // molecule

	return;
}


// calc field using wolf sum (JCP 124 234104 (2006) equation 19
void System::thole_field_wolf() {

	Molecule  * molecule_ptr = nullptr;
	Atom      * atom_ptr     = nullptr;
	Pair      * pair_ptr     = nullptr;
	
	double      r            = 0,
	            rr           = 0, // 1/r (reciprocal of r)
	            R            = pbc.cutoff,
	            rR           = 1./R;
	//used for polar_wolf_alpha (aka polar_wolf_damp)
	double      a            = polar_wolf_alpha,
	            erR          = erfc(a*R), //complementary error functions
	            cutoffterm   = (erR*rR*rR + 2.0*a*OneOverSqrtPi*exp(-a*a*R*R)*rR),
	            bigmess      = 0;

	//init lookup table if needed
	if( polar_wolf_alpha_lookup && !(polar_wolf_alpha_table) )
		polar_wolf_alpha_table = polar_wolf_alpha_lookup_init();

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if ( molecule_ptr == pair_ptr->molecule ) 
					continue; //don't let molecules polarize themselves
				if ( pair_ptr->frozen ) 
					continue; //don't let the MOF polarize itself

				r = pair_ptr->rimg;

				if((r - SMALL_dR < pbc.cutoff) && (r != 0.)) {
					rr = 1./r;

					//we will need this shit if wolf alpha != 0 
					if ( (a != 0) & polar_wolf_alpha_lookup )
						bigmess=polar_wolf_alpha_getval(r);
					else if ( a != 0 ) //no lookup  
						bigmess=(erfc(a*r)*rr*rr+2.0*a*OneOverSqrtPi*exp(-a*a*r*r)*rr);

					for( int p=0; p<3; p++ ) { 
						//see JCP 124 (234104)
						if ( a == 0 ) {
							atom_ptr->ef_static[p] += (pair_ptr->atom->charge)*(rr*rr-rR*rR)*pair_ptr->dimg[p]*rr;
							pair_ptr->atom->ef_static[p] -= (atom_ptr->charge)*(rr*rr-rR*rR)*pair_ptr->dimg[p]*rr;
						} else {
							atom_ptr->ef_static[p] += pair_ptr->atom->charge*(bigmess-cutoffterm)*pair_ptr->dimg[p]*rr;
							pair_ptr->atom->ef_static[p] -= atom_ptr->charge*(bigmess-cutoffterm)*pair_ptr->dimg[p]*rr;
						}
					
					}

				} // no lookup table
			} // pair
		} // atom
	} // molecule

	return;
}


//only calculate the static e-field via ewald
//see http://www.pages.drexel.edu/~cfa22/msim/node50.html
void System::ewald_estatic() {

	//calculate static e-field
	// no need to zero out dipoles; this is done in polar.c
	recip_term();
	real_term();

	return;
}


// the point of this is to store all the polar_wolf_alpha calculations in a big array, and then just look shit up
// that way we don't need to calculate erfc's and exp's over and over and over
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double * System::polar_wolf_alpha_lookup_init () {
	double * rval = nullptr;
	int      i    = 0;
	double   a    = polar_wolf_alpha,
	         r    = 0,
	         rr   = 0;

	polar_wolf_alpha_table_max = (int)( ceil( polar_wolf_alpha_lookup_cutoff )) * 1000;
	SafeOps::malloc( rval, sizeof(double) * polar_wolf_alpha_table_max, __LINE__, __FILE__ );

	for ( i=1; i< polar_wolf_alpha_table_max; i++ ) {
		r = (double)i/1000.;
		rr = 1.0/r;
		rval[i] = erfc(a*r)*rr*rr+2.0*a*OneOverSqrtPi*exp(-a*a*r*r)*rr;
	}
	//store the zeroth component without blowing up
	rval[0]=rval[1];
	

	return rval;
}


double System::polar_wolf_alpha_getval( double r ) {
	
	int i = (int)(r * 1000);
	
	if ( i >= polar_wolf_alpha_table_max ) 
		return 0.0; //answer will be zero if cutoff is large enough

	return polar_wolf_alpha_table[i];
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// iterative solver of the dipole field tensor returns the number of iterations required 
int System::thole_iterative() {

	int     N                 = natoms,
	        iteration_counter = 0, 
	        keep_iterating    = 1;
	Atom ** aa                = atom_array;
	int   * ranked_array      = nullptr;

	
	

	// array for ranking
	SafeOps::calloc( ranked_array, N, sizeof(int), __LINE__, __FILE__ );
	for( int i = 0; i < N; i++)
		ranked_array[i] = i;

	//set all dipoles to alpha*E_static * polar_gamma
	init_dipoles();

	// if ZODID is enabled, then stop here and just return the alpha*E dipoles
	if(polar_zodid) {
		free(ranked_array);
		return(0);
	}

	// iterative solver of the dipole field equations 
	keep_iterating = 1;
	iteration_counter = 0;
	while (keep_iterating) {
		iteration_counter++;

		// divergence detection
		// if we fail to converge, then return dipoles as alpha*E
		if(iteration_counter >= MAX_ITERATION_COUNT && polar_precision) {
			for( int i = 0; i < N; i++ )
				for( int p = 0; p < 3; p++ ) {
					aa[i]->mu[p] = aa[i]->polarizability * (aa[i]->ef_static[p] + aa[i]->ef_static_self[p]);
					aa[i]->ef_induced_change[p] = 0.0; //so we don't break palmo
				}
			//set convergence failure flag
			iterator_failed = 1;
			
			free(ranked_array);
			return(iteration_counter);
		}

		//zero out induced e-field
		for( int i=0; i<natoms; i++ ) 
			for( int p=0; p<3; p++ ) 
				aa[i]->ef_induced[p] = 0;

		//save the current dipole information if we want to calculate precision (or if needed for relaxation)
		if( polar_rrms || polar_precision > 0 || polar_sor || polar_esor ) { 
			for( int i = 0; i < N; i++) 
				for( int p = 0; p < 3; p++) 
					aa[i]->old_mu[p] = aa[i]->mu[p];
		}

		// contract the dipoles with the field tensor (gauss-seidel/gs-ranked optional)
		contract_dipoles(ranked_array);

		if( polar_rrms || polar_precision > 0 )
			calc_dipole_rrms();

		// determine if we are done... 
		keep_iterating = are_we_done_yet(iteration_counter);

		// if we would be finished, we want to contract once more to get the next induced field for palmo
		if( polar_palmo && !keep_iterating ) 
			palmo_contraction(ranked_array);

		//new gs_ranking if needed
		if ( polar_gs_ranked && keep_iterating )
			update_ranking( ranked_array );

		// save the dipoles for the next pass
		for( int i = 0; i < N; i++) {
			for( int p = 0; p < 3; p++) {
				// allow for different successive over-relaxation schemes
				if( polar_sor )
					aa[i]->mu[p] = polar_gamma*aa[i]->new_mu[p] + (1.0 - polar_gamma)*aa[i]->old_mu[p];
				else if( polar_esor )
					aa[i]->mu[p] = (1.0 - exp(-polar_gamma*iteration_counter))*aa[i]->new_mu[p] + exp(-polar_gamma*iteration_counter)*aa[i]->old_mu[p];
				else
					aa[i]->mu[p] = aa[i]->new_mu[p];
			}
		}

	} //end iterate
	free(ranked_array);

	// return the iteration count
	return(iteration_counter);
}


// set them to alpha*E_static
void System::init_dipoles() {

	Atom ** aa = atom_array;

	for( int i=0; i<natoms; i++ ) {
		for( int p=0; p<3; p++ ) {
			aa[i]->mu[p] = aa[i]->polarizability*(aa[i]->ef_static[p]+aa[i]->ef_static_self[p]);
			// should improve convergence since mu's typically grow as induced fields are added in
			if( !polar_sor && !polar_esor )
				aa[i]->mu[p] *= polar_gamma; 
		}
	}
	return;
}


void System::contract_dipoles( int * ranked_array ) {
	int     ii    = 0,
	        jj    = 0,
	        index = 0;
	Atom ** aa    = atom_array;

	for( int i = 0; i < natoms; i++) {
		index = ranked_array[i]; //do them in the order of the ranked index
		ii = index*3;
		if ( aa[index]->polarizability == 0 ) { //if not polar
			//aa[index]->ef_induced[p] is already 0
			aa[index]->new_mu[0] = aa[index]->new_mu[1] = aa[index]->new_mu[2] = 0; //might be redundant?
			aa[index]->mu[0] = aa[index]->mu[1] = aa[index]->mu[2] = 0; //might be redundant?
			continue;
		}	
		for( int j = 0; j < natoms; j++) {
			jj = j*3;
			if(index != j) 
				for( int p = 0; p < 3; p++)
					aa[index]->ef_induced[p] -= UsefulMath::dddotprod( (A_matrix[ii+p]+jj), aa[j]->mu );
		} // end j 

		// dipole is the sum of the static and induced parts
		for( int p = 0; p < 3; p++ ) {
			aa[index]->new_mu[p] = aa[index]->polarizability*(aa[index]->ef_static[p] + aa[index]->ef_static_self[p] + aa[index]->ef_induced[p]);

			// Gauss-Seidel
			if( polar_gs || polar_gs_ranked )
				aa[index]->mu[p] = aa[index]->new_mu[p];
		}

	} // end matrix multiply

	return;
}


void System::palmo_contraction( int * ranked_array ) {

	int     ii    = 0,
	        jj    = 0,
	        index = 0,
	        N     = natoms;
	Atom ** aa    = atom_array; 

	// calculate change in induced field due to this iteration
	for( int i = 0; i < N; i++ ) {
		index = ranked_array[i];
		ii = index*3;

		for( int p=0; p<3; p++ )
			aa[index]->ef_induced_change[p] = -aa[index]->ef_induced[p];

		for( int j = 0; j < N; j++ ) {
			jj = j*3;
			if(index != j) 
				for( int p = 0; p < 3; p++ )
					aa[index]->ef_induced_change[p] -= UsefulMath::dddotprod( A_matrix[ii+p]+jj, aa[j]->mu );
		}
	} 

	return;
}


void System::update_ranking( int * ranked_array ) {

	int     sorted = 0,
	        tmp    = 0,
	        N      = natoms; 
	Atom ** aa     = atom_array;

	// rank the dipoles by bubble sort
	if( polar_gs_ranked ) {
		for( int i = 0; i < N; i++ ) {
			sorted = 1;
			for( int j = 0; j < (N-1); j++ ) {
				if(aa[ranked_array[j]]->rank_metric < aa[ranked_array[j+1]]->rank_metric) {
					sorted = 0;
					tmp = ranked_array[j];
					ranked_array[j] = ranked_array[j+1];
					ranked_array[j+1] = tmp;
				}
			}
			if(sorted) 
				break;
		}
	}

	return;
}


// invert the A matrix
void System::thole_bmatrix() {

	Molecule * molecule_ptr = nullptr;
	Atom     * atom_ptr     = nullptr;
	int        N            = 0;

	// count the number of atoms
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			++N;

	UsefulMath::invert_matrix( 3*N, A_matrix, B_matrix );
}


// get the dipoles by vector matrix multiply
void System::thole_bmatrix_dipoles() {

	int      ii          = 0, 
	         N           = natoms;
	double * mu_array    = nullptr,
	       * field_array = nullptr;

	// allocate working arrays
	SafeOps::calloc( mu_array,    3*N, sizeof(double), __LINE__, __FILE__ );
	SafeOps::calloc( field_array, 3*N, sizeof(double), __LINE__, __FILE__ );
	
	// copy the field in
	for( int i = 0; i < N; i++ ) {
		ii = i*3;
		for( int p = 0; p < 3; p++ )
			field_array[ii+p] = atom_array[i]->ef_static[p] + atom_array[i]->ef_static_self[p];
	}

	// multiply the supervector with the B matrix
	for( int i = 0; i < 3*N; i++ )
		for( int j = 0; j < 3*N; j++ )
			mu_array[i] += B_matrix[i][j]*field_array[j];

	// copy the dipoles out
	for( int i = 0; i < N; i++ ) {
		ii = i*3;
		for( int p = 0; p < 3; p++ )
			atom_array[i]->mu[p] = mu_array[ii+p];
	}

	// free the working arrays
	free(mu_array);
	free(field_array);

}



// calculate the molecular polarizability tensor from the B matrix 
void System::thole_polarizability_tensor() {
	char linebuf[ maxLine ];
	int ii, jj, N;

	double isotropic;

	N = checkpoint->thole_N_atom;

	// clear the polarizability tensor
	for( int p = 0; p < 3; p++)
		for( int q = 0; q < 3; q++)
			C_matrix[p][q] = 0;

	// sum the block terms for the 3x3 molecular tensor
	for( int p = 0; p < 3; p++ ) {
		for( int q = 0; q < 3; q++ ) {
			for( int i = 0; i < N; i++ ) {
				ii = i*3;
				for( int j = 0; j < N; j++ ) {
					jj = j*3;
					C_matrix[p][q] += B_matrix[ii+p][jj+q];
				}
			}
		}
	}

	// get the isotropic term
	isotropic = 0;
	for( int p = 0; p < 3; p++ )
		isotropic += C_matrix[p][p];
	isotropic /= 3.0;

	Output::out("POLARIZATION: polarizability tensor (A^3):\n");
	Output::out("##########################\n");
	for( int p = 0; p < 3; p++ ) {
		for( int q = 0; q < 3; q++ ) {
			sprintf( linebuf, "%.4f ", C_matrix[p][q] );
			Output::out(linebuf);
		}
		Output::out("\n");
	}

	Output::out( "##########################\n" );
	sprintf( linebuf, "isotropic = %.4f\n", isotropic );
	Output::out(linebuf);
	sprintf( linebuf, "XX/ZZ = %.4f\n", C_matrix[0][0]/C_matrix[2][2] );
	Output::out(linebuf);

}