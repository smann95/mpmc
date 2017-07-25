#include <cstring>

#include "Output.h"
#include "Pair.h"
#include "SafeOps.h"
#include "System.h"
#include "UsefulMath.h"


// returns the total potential energy for the system and updates our observables
double System::energy() {

	double potential_energy  = 0, 
	       rd_energy         = 0, 
	       coulombic_energy  = 0,
	       polar_energy      = 0,
	       vdw_energy        = 0,
	       three_body_energy = 0,
	       kinetic_energy    = 0;

	#ifdef POLARTIMING
		static double timing = 0;
		static int count = 0;
	#endif

	natoms = countNatoms();

	// get the pairwise terms necessary for the energy calculation 
	pairs();

	// Only on the first simulation step, make sure that all recalculate flags are set if we made a 
	// volume change (or just reverted from one) set recalculate flags OR if replaying a trajectory 
	// we set last_volume at the end of this function
	if(   last_volume != pbc.volume   ||   ensemble == ENSEMBLE_REPLAY   ||   observables->energy == 0.0  )
		flag_all_pairs();

	// get the electrostatic potential
	if(  ! ( use_sg || rd_only )  ) {

		if( spectre )
			coulombic_energy = coulombic_nopbc( molecules );
		else if( gwp ) {
			coulombic_energy = coulombic_nopbc_gwp();
			kinetic_energy   = coulombic_kinetic_gwp();
			observables->kinetic_energy = kinetic_energy;
		} else
			coulombic_energy = coulombic();

		observables->coulombic_energy = coulombic_energy;

		// get the polarization potential
		if( polarization ) {

			#ifdef POLARTIMING
				// get timing of polarization energy function for cuda comparison 
				Output::GetTimeOfDay( &old_time );
			#endif

			#ifdef CUDA
				if(system->cuda)
					polar_energy = (double)polar_cuda(system);
				else
					polar_energy = polar(system);
			#else
				polar_energy = polar();
			#endif // CUDA 

			#ifdef POLARTIMING
				Output::GetTimeOfDay( &new_time );
				timing = timing * (double)count/((double)count+1.0) 
					+ (double)((new_time.tv_sec-old_time.tv_sec)*1e6+(new_time.tv_usec-old_time.tv_usec)) * 1.0/((double)count +1.0);
				count++;
				if ( system->corrtime ) {
					if ( count % system->corrtime == 0 ) sprintf(linebuf, "OUTPUT: Polarization energy function took %lf us\n", timing);
					output(linebuf);
				}
				else	{
					sprintf(linebuf, "OUTPUT: Polarization energy function took %lf us\n", timing);
					output(linebuf);
				}
			#endif

			observables->polarization_energy = polar_energy;

		}
		if( polarvdw ) {
			#ifdef CUDA
				if (system->cuda) {
					error("error: cuda polarvdw not yet implemented!\n");
					die(-1);
				}
				else
				vdw_energy = vdw(system);
			#else
				vdw_energy = vdw();
			#endif
				observables->vdw_energy = vdw_energy;
		}

	}



	// get the repulsion/dispersion potential
	if( rd_anharmonic )
		rd_energy = anharmonic();
	else if( use_sg )
		rd_energy = sg();
	else if( use_dreiding )
		rd_energy = dreiding();
	else if( using_lj_buffered_14_7 )
		rd_energy = lj_buffered_14_7();
	else if( using_disp_expansion )
		rd_energy = disp_expansion();
	else if( cdvdw_exp_repulsion )
		rd_energy = exp_repulsion();
	else if( !gwp )
		rd_energy = lj();
	observables->rd_energy = rd_energy;
    
	if( using_axilrod_teller )
	{
		three_body_energy = axilrod_teller();
		observables->three_body_energy = three_body_energy;
	}

	// sum the total potential energy 
	potential_energy = rd_energy + coulombic_energy + polar_energy + vdw_energy + three_body_energy;
	
	// not truly potential, but stick it there for convenience of MC 
	//
	//    POSSIBLE BUG: kinetic_energy was uninitialized, and previously only given a value inside the conditional: 
	//
	//        if(!(system->use_sg || system->rd_only)) {}
	//
	//    If this conditional fails, but (system->gwp) is true (idk if this is possible), an un-initialized value would have been
	//    added to potential_energy. Now, 0 is being added, but am not sure if this is the desired behavior. -bt
	 
	
	if( gwp )
		potential_energy += kinetic_energy;
	observables->energy = potential_energy;

	countN();
	observables->spin_ratio /= observables->N;

	// for NVE
	if(ensemble == ENSEMBLE_NVE) {
		observables->kinetic_energy = total_energy - potential_energy;
		observables->temperature = (2.0/3.0) * observables->kinetic_energy/observables->N;
	}

	// need this for the isosteric heat 
	observables->NU = observables->N * observables->energy;

	// set last known volume
	last_volume = pbc.volume;

	if( cavity_autoreject_absolute )
		potential_energy += cavity_absolute_check();

	return potential_energy;
}


//returns interaction VDW energy
double System::vdw() {

	int       N;              // number of atoms
	double    e_total,        // total energy
	          e_iso;          // isolation energy (atoms @ infinity)
	double  * sqrtKinv;       // matrix K^(-1/2); cholesky decomposition of K
	double ** Am = A_matrix;  // A_matrix
	mtx_t   * Cm;             // C_matrix (we use single pointer due to LAPACK requirements)
	double  * eigvals;        // eigenvales
	double    fh_corr, lr_corr;

	N = natoms;

	//allocate arrays. sqrtKinv is a diagonal matrix. d,e are used for matrix diag.
	sqrtKinv = getsqrtKinv( N );

	//calculate energy vdw of isolated molecules
	e_iso = sum_eiso_vdw ( sqrtKinv );

	//Build the C_Matrix
	Cm = build_M (3*N, 0, Am, sqrtKinv);

	//setup and use lapack diagonalization routine dsyev_()
	eigvals = lapack_diag (Cm, polarvdw ); //eigenvectors if system->polarvdw == 2
	if( polarvdw == 2 )
		printevects(Cm);

	//return energy in inverse time (a.u.) units
	e_total = eigen2energy( eigvals, Cm->dim, temperature );
	e_total *= au2invseconds * half_hBar; //convert a.u. -> s^-1 -> K

	//vdw energy comparison
	if ( polarvdw == 3 )
		printf("VDW Two-Body | Many Body = %lf | %lf\n", twobody(), e_total-e_iso );

	if( feynman_hibbs ) {
		if( vdw_fh_2be ) fh_corr = fh_vdw_corr_2be(); //2be method
		else fh_corr = fh_vdw_corr(); //mpfd
	}
	else fh_corr=0;

	if ( rd_lrc ) lr_corr = lr_vdw_corr();
	else lr_corr=0;

	//cleanup and return
	free(sqrtKinv);
	free(eigvals);
	free_mtx(Cm);

	return e_total - e_iso + fh_corr + lr_corr;

}


//build the matrix K^(-1/2) -- see the PDF
double * System::getsqrtKinv( int N ) {
	double   * sqrtKinv;
	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	int i = 0;

	//malloc 3*N wastes an insignificant amount of memory, but saves us a lot of index management
	 SafeOps::malloc( sqrtKinv, 3 * N * sizeof(double), __LINE__, __FILE__ );
	
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			//Seek through atoms, calculate sqrtKinv*
			sqrtKinv[i] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			sqrtKinv[i+1] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			sqrtKinv[i+2] = sqrt(atom_ptr->polarizability)*atom_ptr->omega;
			i+=3;
		}
	}
	
	return sqrtKinv;
}


//go through each molecule and determine the VDW energy associated with each isolated molecule
double System::sum_eiso_vdw ( double * sqrtKinv ) {

	char       linebuf[maxLine];
	double     e_iso = 0;
	Molecule * mp;
	vdw_t    * vp;
	vdw_t    * vpscan;

	//loop through molecules. if not known, calculate, store and count. otherwise just count.
	for ( mp = molecules; mp; mp=mp->next ) {
		for ( vp = vdw_eiso_info; vp != NULL; vp=vp->next ) { //loop through all vp's
			if ( strncmp(vp->mtype,mp->moleculetype, maxLine ) == 0 ) {
					e_iso += vp->energy; //count energy
					break; //break out of vp loop. the current molecule is accounted for now. go to the next molecule
			}
			else continue; //not a match, check the next vp
		} //vp loop

		if ( vp == NULL ) { //if the molecule was unmatched, we need to grow the list
			// end of vp list and we haven't matched yet -> grow vdw_eiso_info
			// scan to the last non-NULL element
			if ( vdw_eiso_info == NULL ) {
				SafeOps::calloc( vdw_eiso_info, 1, sizeof(vdw_t), __LINE__, __FILE__ );
				vpscan = vdw_eiso_info; //set scan pointer
			} else {
				for ( vpscan = vdw_eiso_info; vpscan->next != NULL; vpscan=vpscan->next );
				SafeOps::calloc( vpscan->next, 1, sizeof(vdw_t), __LINE__, __FILE__ ); //allocate space
				vpscan = vpscan->next;
			} //done scanning and malloc'ing
		
			//set values
			strncpy( vpscan->mtype, mp->moleculetype, maxLine ); //assign moleculetype
			vpscan->energy = calc_e_iso( sqrtKinv, mp ); //assign energy
			if ( std::isfinite(vpscan->energy) == 0 ) { //if nan, then calc_e_iso failed
				sprintf(linebuf,"VDW: Problem in calc_e_iso.\n");
				Output::out( linebuf );
				throw infinite_energy_calc;
			}
			//otherwise count the energy and move to the next molecule
			e_iso += vpscan->energy;

		} //vp==NULL
	} //mp loop	

	////all of this logic is actually really bad if we're doing surface fitting, since omega will change... :-(
	//free everything so we can recalc next step
	if( ensemble == ENSEMBLE_SURF_FIT )  {
		free_vdw_eiso( vdw_eiso_info );
		vdw_eiso_info = NULL;
	}
	
	return e_iso;
}

//calculate energies for isolated molecules
//if we don't know it, calculate it and save the value
double System::calc_e_iso ( double * sqrtKinv, Molecule * mptr ) {

	int        nstart, nsize; 
	double     e_iso;         //total vdw energy of isolated molecules
	mtx_t    * Cm_iso;        //matrix Cm_isolated
	double   * eigvals;       //eigenvalues of Cm_cm
	Molecule * molecule_ptr;
	Atom     * atom_ptr;

	nstart=nsize=0; //loop through each individual molecule
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr=molecule_ptr->next ) {
		if ( molecule_ptr != mptr ) {  //count atoms then skip to next molecule
			for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) nstart++;
			continue;
		}

		//now that we've found the molecule of interest, count natoms, and calc energy
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) nsize++;

		//build matrix for calculation of vdw energy of isolated molecule
		Cm_iso  = build_M(3*(nsize), 3*nstart, A_matrix, sqrtKinv);
		//diagonalize M and extract eigenvales -> calculate energy
		eigvals = lapack_diag( Cm_iso, 1 ); //no eigenvectors
		e_iso   = eigen2energy( eigvals, Cm_iso->dim, temperature );

		//free memory
		free(eigvals);
		free_mtx(Cm_iso);

		//convert a.u. -> s^-1 -> K
		return e_iso * au2invseconds * half_hBar;
	}	

	//unmatched molecule
	return NAN; //we should never get here
}


//free vdw pointer which keeps track of e_iso energies
void System::free_vdw_eiso(vdw_t * vdw_eiso_info) {

	vdw_t * vp;
	vdw_t ** varray = NULL;
	int i=0;

	for ( vp = vdw_eiso_info; vp; vp=vp->next ) {
		SafeOps::realloc( varray, sizeof(vdw_t *)*(i+1), __LINE__, __FILE__ );
		varray[i]=vp;
		i++;
	}

	while ( i>0 ) {
		i--;
		free(varray[i]);
	}

	free(varray);

	return;
}


//build C matrix for a given molecule/system, with atom indicies (offset)/3..(offset+dim)/3
System::mtx_t * System::build_M ( int dim, int offset, double ** Am, double * sqrtKinv ) {
	int i; //dummy
	int iA, jA; //Am indicies
	int iC, jC; //Cm indicies
	int nonzero; //non-zero col/rows in Am
	mtx_t * Cm; //return matrix Cm

	//count non-zero elements
	nonzero=0;
	for ( i=offset; i<dim+offset; i++ )
		if ( sqrtKinv[i] != 0 ) nonzero++;

	//allocate
	Cm = alloc_mtx(nonzero);

	// build lapack compatible matrix from Am[offset..dim, offset..dim]
	iC=jC=-1; //C index
	for ( iA=offset; iA<dim+offset; iA++ ) {
		if ( sqrtKinv[iA] == 0 ) continue; //suppress rows/cols full of zeros
		iC++; jC=-1;
		for ( jA=offset; jA<=iA; jA++ ) {
			if ( sqrtKinv[jA] == 0 ) continue; //suppress
			jC++;
			(Cm->val)[iC+jC*(Cm->dim)]=
				Am[iA][jA]*sqrtKinv[iA]*sqrtKinv[jA];
		}
	}

	return Cm;
}


System::mtx_t * System::alloc_mtx( int dim ) {
	//alloc matrix variable and set dim
	mtx_t * M = nullptr;
	SafeOps::malloc( M, sizeof(mtx_t), __LINE__, __FILE__ );
	M->dim=dim;
	//alloc matrix storage space
	SafeOps::calloc( M->val, dim*dim, sizeof(double), __LINE__, __FILE__ );
	
	return M;
}


void System::free_mtx( mtx_t * M ) {
	free(M->val);
	free(M);
	return;
}


void System::printevects( mtx_t * M ) {
	int r,c;

	printf("%%vdw === Begin Eigenvectors ===\n");
	for ( r=0; r < (M->dim); r++ ) {
		for ( c=0; c < (M->dim); c++ ) {
			printf("%.2le ", (M->val)[r+c*M->dim]);
		}
		printf("\n");
	}
	printf("%%vdw=== End Eigenvectors ===\n");
}

// not needed unless T >> 300 
//double wtanh ( double w, double T ) {
//	if ( w < 1.0E-10 ) TWOoverHBAR*T/au2invsec; //from Taylor expansion
//	if ( T == 0 ) return w;
//	return w/tanh(halfHBAR*w*au2invsec/T);
//}
double System::eigen2energy( double * eigvals, int dim, double temperature ) {
	int i;
	double rval=0;
	
	if ( eigvals == NULL ) return 0;

	for ( i=0; i<dim; i++ ) {
		if ( eigvals[i] < 0 ) eigvals[i]=0;
		//rval += wtanh(sqrt(eigvals[i]), temperature);
		rval += sqrt(eigvals[i]);
	}
	return rval;
}

//with damping
double System::twobody() {
	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;
	double     energy = 0;

	//for each pair
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > pbc.cutoff ) continue;
				//check if fh is non-zero
				if ( atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 || 
					atom_ptr->omega == 0 || pair_ptr->atom->omega == 0 ) continue; //no vdw energy

				//calculate two-body energies
				energy += e2body(atom_ptr,pair_ptr,pair_ptr->rimg);
			}
		}
	}

	return energy;
}

//calculate T matrix element for a particular separation
double System::e2body( Atom * atom, Pair * pair, double r ) {
	double energy;
	double lr  = polar_damp * r;
	double lr2 = lr*lr;
	double lr3 = lr*lr2;
	double Txx = pow(r,-3)*(-2.0+(0.5*lr3+lr2+2*lr+2)*exp(-lr));
	double Tyy = pow(r,-3)*(1-(0.5*lr2+lr+1)*exp(-lr));
	double * eigvals;
	mtx_t * M = alloc_mtx(6);
	
	//only the sub-diagonals are non-zero
	M->val[1]=M->val[2]=M->val[4]=M->val[5]=M->val[6]=M->val[8]=M->val[9]=M->val[11]=0;
	M->val[12]=M->val[13]=M->val[15]=M->val[16]=M->val[19]=M->val[20]=M->val[22]=M->val[23]=0;
	M->val[24]=M->val[26]=M->val[27]=M->val[29]=M->val[30]=M->val[31]=M->val[33]=M->val[34]=0;

	//true diagonals
	M->val[0]=M->val[7]=M->val[14]=(atom->omega)*(atom->omega);
	M->val[21]=M->val[28]=M->val[35]=(pair->atom->omega)*(pair->atom->omega);

	//sub-diagonals
	M->val[3]=M->val[18]=
		(atom->omega)*(pair->atom->omega)*sqrt(atom->polarizability*pair->atom->polarizability)*Txx;
	M->val[10]=M->val[17]=M->val[25]=M->val[32]=
		(atom->omega)*(pair->atom->omega)*sqrt(atom->polarizability*pair->atom->polarizability)*Tyy;

	eigvals=lapack_diag(M,1);
	energy = eigen2energy(eigvals, 6, temperature);

	//subtract energy of atoms at infinity
//	energy -= 3*wtanh(atom->omega, system->temperature);
	energy -= 3*atom->omega;
//	energy -= 3*wtanh(pair->atom->omega, system->temperature);
	energy -= 3*pair->atom->omega;

	free(eigvals);
	free_mtx(M);

  return energy * au2invseconds * half_hBar;
}

//    LAPACK using 1D arrays for storing matricies.
//    / 0  3  6 \
//    | 1  4  7 |    =    [ 0 1 2 3 4 5 6 7 8 ]
//    \ 2  5  8 /									
double * System::lapack_diag ( mtx_t * M, int jobtype ) {
	
	char     job;      //job type
	char     uplo='L'; //operate on lower triagle
	double * work;     //working space for dsyev
	int      lwork;    //size of work array
	int      rval=0;   //returned from dsyev_
	double * eigvals;
	char     linebuf[maxLine];

	//eigenvectors or no?
	if ( jobtype == 2 )
		job='V';
	else
		job = 'N';

	if ( M->dim == 0 ) return NULL;

	//allocate eigenvalues array
	SafeOps::malloc( eigvals, M->dim*sizeof(double), __LINE__, __FILE__ );
	
	//optimize the size of work array
	lwork = -1;
	SafeOps::malloc( work, sizeof(double), __LINE__, __FILE__ );
	//dsyev_(&job, &uplo, &(M->dim), M->val, &(M->dim), eigvals, work, &lwork, &rval); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//now optimize work array size is stored as work[0]
	lwork = (int)work[0];
	SafeOps::realloc( work, lwork*sizeof(double), __LINE__, __FILE__ );
	//diagonalize
	//dsyev_(&job, &uplo, &(M->dim), M->val, &(M->dim), eigvals, work, &lwork, &rval); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if ( rval != 0 ) {
		sprintf(linebuf,"error: LAPACK: dsyev returned error: %d\n", rval);
		Output::err(linebuf);
		throw lapack_error;
	}

	free(work);

	return eigvals;
}

// long-range correction
double System::lr_vdw_corr() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double w1, w2;   //omegas
	double a1, a2;   //alphas
	double cC;       //leading coefficient to r^-6
	double corr = 0; //correction to the energy

	//skip if PBC isn't set-up
	if ( pbc.volume == 0 ) {
		Output::err( "VDW: PBC not set-up. Did you define your basis? Skipping LRC.\n" );
		return 0;
	}

	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
					//skip if frozen
					if ( pair_ptr->frozen ) continue;
					// skip if same molecule  // don't do this... this DOES contribute to LRC
					// if ( molecule_ptr == pair_ptr->molecule ) continue;
					// fetch alphas and omegas
					a1 = atom_ptr->polarizability;
					a2 = pair_ptr->atom->polarizability;
					w1 = atom_ptr->omega;
					w2 = pair_ptr->atom->omega;
					if ( w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0 ) continue; //no vdw energy
					// 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
					cC=1.5 * c_hBar * w1*w2/(w1+w2) * au2invseconds * a1 * a2;

					// long-range correction
					corr += -4.0/3.0 * pi * cC * pow(pbc.cutoff,-3) / pbc.volume;
			}
		}
	}

	return corr;
}

// feynman-hibbs correction - molecular pair finite differencing method
double System::fh_vdw_corr() {

	const double FINITE_DIFF = 0.01; //too small -> vdw calc noises becomes a problem

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;
	double     rm;                //reduced mass
	double     E[5];              //energy at five points, used for finite differencing
	double     dv, d2v, d3v, d4v; //derivatives
	double     corr = 0;          //correction to the energy
	double     corr_single;       //single vdw interaction energy
	double     h = FINITE_DIFF ;  //small dr used for finite differencing //too small -> vdw calculation noise becomes a problem

	//for each pair
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) {
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) {
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) {
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > pbc.cutoff ) continue;
				//check if fh is non-zero
				if ( atom_ptr->polarizability == 0 || pair_ptr->atom->polarizability == 0 || 
					atom_ptr->omega == 0 || pair_ptr->atom->omega == 0 ) continue; //no vdw energy

				//calculate two-body energies
				E[0]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg-h-h); //smaller r
				E[1]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg-h); 
				E[2]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg); //current r
				E[3]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg+h); //larger r
				E[4]=e2body(atom_ptr,pair_ptr,pair_ptr->rimg+h+h);

				//derivatives (Numerical Methods Using Matlab 4E 2004 Mathews/Fink 6.2)
				dv = (E[3]-E[1])/(2.0*h);
				d2v = (E[3]-2.0*E[2]+E[1])/(h*h);
				d3v = (E[4]-2*E[3]+2*E[1]-E[0])/(2*pow(h,3));
				d4v = (E[4]-4*E[3]+6*E[2]-4*E[1]+E[0])/pow(h,4);
				
				// reduced mass
				rm=AMU2KG*(molecule_ptr->mass)*(pair_ptr->molecule->mass)/
					((molecule_ptr->mass)+(pair_ptr->molecule->mass));

				//2nd order correction
				corr_single = pow(METER2ANGSTROM, 2) * (hBar*hBar/(24.0*kB*temperature*rm)) * (d2v + 2.0*dv/pair_ptr->rimg);
				//4th order correction
				if( feynman_hibbs_order >= 4 )
					corr_single += pow(METER2ANGSTROM, 4)*(pow(hBar, 4) /
						(1152.0*pow(kB*temperature*rm, 2))) *
						(15.0*dv/pow(pair_ptr->rimg, 3) + 4.0*d3v/pair_ptr->rimg + d4v);

				corr += corr_single;
			}
		}
	}

	return corr;
}

// feynman-hibbs using 2BE (shitty)
double System::fh_vdw_corr_2be() {

	Molecule * molecule_ptr;
	Atom     * atom_ptr;
	Pair     * pair_ptr;

	double rm; //reduced mass
	double w1, w2; //omegas
	double a1, a2; //alphas
	double cC; //leading coefficient to r^-6
	double dv, d2v, d3v, d4v; //derivatives
	double corr = 0; //correction to the energy
	double corr_single; //single vdw interaction energy

	//for each pair
	for ( molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next ) { 
		for ( atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next ) { 
			for ( pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next ) { 
				//skip if frozen
				if ( pair_ptr->frozen ) continue;
				//skip if they belong to the same molecule
				if ( molecule_ptr == pair_ptr->molecule ) continue;
				//skip if distance is greater than cutoff
				if ( pair_ptr->rimg > pbc.cutoff ) continue;
				//fetch alphas and omegas
				a1=atom_ptr->polarizability;
				a2=pair_ptr->atom->polarizability;
				w1=atom_ptr->omega;
				w2=pair_ptr->atom->omega;

				if ( w1 == 0 || w2 == 0 || a1 == 0 || a2 == 0 ) continue; //no vdw energy
				// 3/4 hbar/k_B(Ks) omega(s^-1)  Ang^6
				cC=1.5 * c_hBar * w1*w2/(w1+w2) * au2invseconds * a1 * a2;
				// reduced mass
				rm = AMU2KG * (molecule_ptr->mass) * (pair_ptr->molecule->mass)/
					((molecule_ptr->mass)+(pair_ptr->molecule->mass));

				//derivatives 
				dv = 6.0*cC*pow(pair_ptr->rimg,-7);
				d2v= dv * (-7.0)/pair_ptr->rimg;
				if( feynman_hibbs_order >= 4 ) {
					d3v= d2v* (-8.0)/pair_ptr->rimg;
					d4v= d3v* (-9.0)/pair_ptr->rimg;
				}

				//2nd order correction
				corr_single = pow(METER2ANGSTROM, 2)*(hBar*hBar/(24.0*kB*temperature*rm))*(d2v + 2.0*dv/pair_ptr->rimg);
				//4th order correction
				if ( feynman_hibbs_order >= 4 )
				corr_single += pow(METER2ANGSTROM, 4)*(pow(hBar, 4) /
				(1152.0*pow(kB*temperature*rm, 2))) *
				(15.0*dv/pow(pair_ptr->rimg, 3) + 4.0*d3v/pair_ptr->rimg + d4v);

				corr += corr_single;
			}
		}
	}

	return corr;
}

// energy of molecule in a 1D anharmonic well
double System::anharmonic() {

	double     k, g, x;
	double     energy;
	Molecule * molecule_ptr;
	Atom     * atom_ptr;

	k = rd_anharmonic_k;
	g = rd_anharmonic_g;

	energy = 0;
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			x = atom_ptr->pos[0];
			energy += anharmonic_energy(k, g, x);

			if( feynman_hibbs ) {

				if( feynman_kleinert ) {

					energy = anharmonic_fk( temperature, atom_ptr->mass, k, g, x );

				} else {

					if( feynman_hibbs_order == 2 )
						energy += anharmonic_fh_second_order( temperature, atom_ptr->mass, k, g, x );
					else if( feynman_hibbs_order == 4 )
						energy += anharmonic_fh_fourth_order( temperature, atom_ptr->mass, k, g, x );
				}
			}
		}
	}

	return(energy);
}

// anharmonic potential
double System::anharmonic_energy(double k, double g, double x) {

	double potential;

	potential =   0.5 * k * pow(x, 2)    +    0.25 * g * pow(x, 4);

	return potential;
}

// Feynman-Kleinert iterative method of effective potential */
double System::anharmonic_fk(double temperature, double mass, double k, double g, double x) {

	int    keep_iterating;
	double a_sq;               // width a^2 (A^2) 
	double omega, omega_sq;    // spacing Omega^2 (K/A^2) 
	double prev_a_sq;          // last a_sq 
	double tolerance;          // iterative tolerance 
	double V_a;                // V_a^2 (K) 
	double potential;          // W_1 (K) 
	double conversion_factor;  // hbar^2*(m2A)^2/(k*m)

	// convert the mass to kg 
	mass *= AMU2KG;

	conversion_factor = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(kB*mass);

	// initial guess a^2 = beta/12 
	a_sq = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(12.0*kB*temperature*mass);

	// solve self-consistently 
	keep_iterating = 1;
	while(keep_iterating) {

		// save the last a_sq for tolerance 
		prev_a_sq = a_sq;

		omega_sq = conversion_factor*(k + 3.0*g*a_sq + 3.0*g*pow(x, 2)); omega = sqrt(omega_sq);
		a_sq = conversion_factor*(temperature/omega_sq)*((omega/(2.0*temperature))*(1.0/tanh(omega/(2.0*temperature))) - 1.0);

		tolerance = fabs(prev_a_sq - a_sq);
		if(tolerance < FEYNMAN_KLEINERT_TOLERANCE)
			keep_iterating = 0;

	}

	V_a = 0.5*a_sq*k + 0.75*g*pow(a_sq, 2) + 0.5*(k + 3.0*g*a_sq)*pow(x, 2) + 0.25*g*pow(x, 4);

	potential = temperature*log(sinh(omega/(2.0*temperature))/(omega/(2.0*temperature))) - 0.5*omega_sq*a_sq/conversion_factor + V_a;

	return potential;
}

// up to FH h^2 effective potential term 
double System::anharmonic_fh_second_order(double temperature, double mass, double k, double g, double x) {

	double first_derivative;
	double second_derivative;
	double potential;

	mass *= AMU2KG;

	first_derivative = k*x + g*pow(x, 3);
	second_derivative = k + 3.0*g*pow(x, 2);

	potential = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(24.0*kB*temperature*mass)*(second_derivative + 2.0*first_derivative/x);

	return potential;
}

// up to FH h^4 effective potential term 
double System::anharmonic_fh_fourth_order(double temperature, double mass, double k, double g, double x) {

	double first_derivative, second_derivative;
	double other_derivatives;
	double potential;

	mass *= AMU2KG;

	first_derivative = k*x + g*pow(x, 3);
	second_derivative = k + 3.0*g*pow(x, 2);
	other_derivatives = 15.0*k/pow(x, 2) + 45.0*g;

	potential = pow(METER2ANGSTROM, 2)*pow(hBar, 2)/(24.0*kB*temperature*mass)*(second_derivative + 2.0*first_derivative/x);
	potential += pow(METER2ANGSTROM, 4)*pow(hBar, 4)/(1152.0*pow(kB*temperature*mass, 2))*other_derivatives;

	return potential;
}