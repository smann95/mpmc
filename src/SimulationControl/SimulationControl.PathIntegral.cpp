#include <time.h>
#include "Output.h"
#include "SimulationControl.h"
#include "SafeOps.h"
extern int rank, size;



bool SimulationControl::check_pi_nvt_options() {

	if( !sys.beads ) {
		if( size < 3 ) {
			Output::out("Path Integrals require at least 3 MPI nodes to run.\n");
			throw insufficient_MPI_size_for_PI;
		}
	} 

	if(sys.beads < 0) {
		Output::out( "Local beads for Path Integral computations must be positive.\n");
		throw invalid_bead_count;
	}

	return ok;
}

void SimulationControl::initialize_PI_NVT_Systems() {

	char linebuf[maxLine] = {0};

	System **temp_sys;
	SafeOps::calloc(temp_sys, size, sizeof(System*), __LINE__, __FILE__);
	for( int i=0; i<size; i++ ) {
		temp_sys[i] = new System(sys);
		systems.push_back(temp_sys[i]);
		systems[i]->setup_simulation_box();
		sprintf( linebuf, "SIM_CONTROL: simulation box %d configured.\n", i );
		Output::out(linebuf);
	}
		
	systems[rank]->allocateStatisticsMem();

	// Set up the pairs list and other considerations only for the system,
	// systems[rank], that this controller will be directly manipulating. 
		
	systems[rank]->allocate_pair_lists();
	Output::out( "SIM_CONTROL: finished allocating pair lists\n" );
				
	systems[rank]->pairs();          // get all of the pairwise interactions, exclusions, etc.
	systems[rank]->flag_all_pairs(); // set all pairs to initially have their energies calculated 
	Output::out( "SIM_CONTROL: finished calculating pairwise interactions\n" );

	if( systems[rank]->cavity_bias )
		systems[rank]->setup_cavity_grid();

	// if polarization active, allocate the necessary matrices 
	if(  systems[rank]->polarization   &&   !systems[rank]->cuda   &&   !systems[rank]->polar_zodid  )
		systems[rank]->thole_resize_matrices();

	// if histogram calculation flag is set, allocate grid
	if( systems[rank]->calc_hist ) {
		systems[rank]->setup_histogram();
	}

	if( systems[rank]->preset_seed )
		systems[rank]->mt_rand.seed( sys.preset_seed );
	else 
		systems[rank]->mt_rand.seed( (unsigned int)time(0) );

		
	if(  ! systems[rank]->use_sg   ||   systems[rank]->rd_only  )
	{
			
		sprintf(linebuf, "SIM_CONTROL: SYS[%d] Ewald gaussian width = %f A\n", rank, systems[rank]->ewald_alpha);
		Output::out(linebuf);
		sprintf(linebuf, "SIM_CONTROL: SYS[%d] Ewald kmax = %d\n", rank, systems[rank]->ewald_kmax);
		Output::out(linebuf);
	}

	// ensure that all SPECTRE charges lie within the restricted domain 
	if( systems[rank]->spectre )
		systems[rank]->spectre_wrapall();
}

bool SimulationControl::pi_nvt_mc() {

	System test[10];

	double  initial_energy  = 0, 
	        final_energy    = 0,
	        current_energy  = 0;


	System::mpiData mpi;
	test[0].setup_mpi_dataStructs( mpi );
/*
	if( sorbateCount > 1 ) 
		SafeOps::calloc( sorbateGlobal, sorbateCount, sizeof(sorbateAverages_t), __LINE__, __FILE__ );
/*
	if( cavity_bias )
		cavity_update_grid(); // update the grid for the first time 

	observables->volume = pbc.volume; // set volume observable

	initial_energy      = energy();   // get the initial energy of the system 


	#ifdef QM_ROTATION
		// solve for the rotational energy levels 
		if(system->quantum_rotation) quantum_system_rotational_energies(system);
	#endif // QM_ROTATION 


	// be a bit forgiving of the initial state 
	if( !std::isfinite(initial_energy) ) 
		initial_energy = observables->energy = MAXVALUE;

	
	// write initial observables to stdout and logs
	if( ! rank ) {

		open_files(); // open output files 

		calc_system_mass();
		// average in the initial values once  (we don't want to double-count the initial state when using MPI)
		update_root_averages( observables, avg_observables );
		// average in the initial sorbate values
		if ( sorbateCount > 1 ) {
			update_sorbate_info(); //local update
			update_root_sorb_averages( sorbateInfo ); //global update
		}
		// write initial observables exactly once
		if( fp_energy ) 
			write_observables( fp_energy, observables, temperature );
		if( fp_energy_csv ) 
			write_observables_csv( fp_energy_csv, observables, temperature );
		Output::out( "MC: initial values:\n" );
		write_averages();
	}

	// save the initial state 
	do_checkpoint();

	
	// main MC loop 
	for( step = 1; step <= numsteps; step++ ) {

		// restore the last accepted energy 
		initial_energy = observables->energy;

		// perturb the system 
		make_move();

		// calculate the energy change 
		final_energy = energy();

		#ifdef QM_ROTATION
			// solve for the rotational energy levels 
			if(system->quantum_rotation && (system->checkpoint->movetype == MOVETYPE_SPINFLIP))
				quantum_system_rotational_energies(system);
		#endif // QM_ROTATION 

		// treat a bad contact as a reject 
		if( !std::isfinite(final_energy) ) {
			observables->energy = MAXVALUE;
			nodestats->boltzmann_factor = 0;
		} 
		else 
			boltzmann_factor( initial_energy, final_energy );

		// Metropolis function 
		if(  (get_rand() < nodestats->boltzmann_factor)   &&   (iterator_failed == 0)  ) {	
			/////////// ACCEPT

			current_energy = final_energy;

			// checkpoint 
			do_checkpoint();
			register_accept();

			// SA 
			if( simulated_annealing )
			{
				if( simulated_annealing_linear==1 )
				{
					temperature = temperature + (simulated_annealing_target - temperature)/(numsteps - step);
					if(  (numsteps - step)   ==   0   )
						temperature = simulated_annealing_target;
				}
				else
					temperature = simulated_annealing_target + (temperature - simulated_annealing_target)*simulated_annealing_schedule;
			}

		} else { 
			// REJECT ///////////////////////////////////////////////////////////////
			//reset the polar iterative failure flag and restore from last checkpoint
			current_energy = initial_energy; //used in parallel tempering
			iterator_failed = 0;
			restore();
			register_reject();
		}

		// perform parallel_tempering
		if ( (parallel_tempering) && ((step % ptemp_freq) == 0) )
			temper_system( current_energy );

		// track the acceptance_rate 
		track_ar( nodestats );

		// each node calculates its stats 
		update_nodestats( nodestats, avg_nodestats );


		// do this every correlation time, and at the very end 
		if(  !(step % corrtime)  ||  (step == numsteps)  ) {

			// copy observables and avgs to the mpi send buffer 
			// histogram array is at the end of the message 
			if( calc_hist ) {
				zero_grid( grids->histogram->grid );
				population_histogram();
			}

			// update frozen and total system mass
			calc_system_mass();

			// update sorbate info on each node
			if ( sorbateCount > 1 )
				update_sorbate_info();

			//write trajectory files for each node -> one at a time to avoid disk congestion
			#ifdef MPI
				for ( int j=0; j<size; j++ ) {
					MPI_Barrier(MPI_COMM_WORLD);
					if( j == rank )
						write_states();
				}
			#else
				write_states();
			#endif

			// restart files for each node -> one at a time to avoid disk congestion
			if( write_molecules_wrapper(pqr_restart) < 0 ) {
				Output::err("MC: could not write restart state to disk\n");
				throw unknown_file_error;
			}

			// dipole/field data for each node -> one at a time to avoid disk congestion
			#ifdef MPI
				if ( polarization ) {
					for ( int j=0; j<size; j++ ) {
						MPI_Barrier(MPI_COMM_WORLD);
						if ( j == rank ) {
							write_dipole();
							write_field();
						}
					}
				}
			#else
				if ( polarization ) {
					write_dipole();
					write_field();
				}
			#endif

			// zero the send buffer 
			memset( mpi.snd_strct, 0, mpi.msgsize );
			memcpy( mpi.snd_strct, observables, sizeof(observables_t) );
			memcpy( mpi.snd_strct + sizeof(observables_t), avg_nodestats, sizeof(avg_nodestats_t) );
			if( calc_hist )
				mpi_copy_histogram_to_sendbuffer(
					mpi.snd_strct + sizeof(observables_t) + sizeof(avg_nodestats_t), 
					grids->histogram->grid
				);
			if ( sorbateCount > 1 )
				memcpy(
					mpi.snd_strct   +   sizeof(observables_t)   +   sizeof(avg_nodestats_t)   +   calc_hist * n_histogram_bins * sizeof(int), //compensate for the size of hist data, if neccessary
					sorbateInfo,
					sorbateCount * sizeof( sorbateInfo_t )
				);

			if( ! rank )
				memset(mpi.rcv_strct, 0, size * mpi.msgsize);

			#ifdef MPI
				MPI_Gather( snd_strct, 1, msgtype, rcv_strct, 1, msgtype, 0, MPI_COMM_WORLD );
				MPI_Gather( &temperature, 1, MPI_DOUBLE, temperature_mpi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
				// need to gather data for sorbate stats also
			#else
				memcpy( mpi.rcv_strct, mpi.snd_strct, mpi.msgsize );
				mpi.temperature[0] = temperature;
			#endif // MPI 

			// head node collects all observables and averages 
			if( !rank ) {
				// clear avg_nodestats to avoid double-counting 
				clear_avg_nodestats();
				//loop for each core -> shift data into variable_mpi, then average into avg_observables
				for( int j = 0; j < size; j++ ) { 
					// copy from the mpi buffer 
					memcpy( mpi.observables,   mpi.rcv_strct + j*mpi.msgsize,                         sizeof(observables_t  ));
					memcpy( mpi.avg_nodestats, mpi.rcv_strct + j*mpi.msgsize + sizeof(observables_t), sizeof(avg_nodestats_t));
					if( calc_hist )
						mpi_copy_rcv_histogram_to_data( mpi.rcv_strct + j*mpi.msgsize + sizeof(observables_t) + sizeof(avg_nodestats_t), grids->histogram->grid );
					if( sorbateCount > 1 )
						memcpy( mpi.sinfo, 
							mpi.rcv_strct   +   j*mpi.msgsize   +   sizeof(observables_t)   +   sizeof(avg_nodestats_t)   +   calc_hist * n_histogram_bins * sizeof(int), //compensate for the size of hist data, if neccessary
							sorbateCount * sizeof(sorbateInfo_t)
						);

					// write observables 
					if( fp_energy ) 
						write_observables( fp_energy, mpi.observables, mpi.temperature[j] );
					if( fp_energy_csv ) 
						write_observables_csv( fp_energy_csv, mpi.observables, mpi.temperature[j] );
					// collect the averages 
					// if parallel tempering, we will collect obserables from the coldest bath. this can't be done for
					// nodestats though, since nodestats are averaged over each corrtime, rather than based on a single 
					// taken at the corrtime 
					update_root_nodestats( mpi.avg_nodestats, avg_observables );
					if ( ! parallel_tempering ) {
						update_root_averages( mpi.observables, avg_observables );
						if( calc_hist )
							update_root_histogram();
						if( sorbateCount > 1 ) 
							update_root_sorb_averages( mpi.sinfo );
					}
					else if ( ptemp->index[j] == 0 ) {
						update_root_averages( mpi.observables, avg_observables );
						if( calc_hist )
							update_root_histogram();
						if( sorbateCount > 1 )
							update_root_sorb_averages( mpi.sinfo );
					}
				}

				output_file_data();

			} // !rank 

		} // corrtime 

	} // main loop 
	


	// restart files for each node
	if(  write_molecules_wrapper(pqr_output) < 0  ) {
		Output::err( "MC: could not write final state to disk\n" );
		throw unknown_file_error;
	}

	if(!rank) {
		free( mpi.rcv_strct   );
		free( mpi.temperature );
	}

	if( sorbateCount > 1 ) 
		free( mpi.sinfo );
	
	free( mpi.snd_strct     );
	free( mpi.observables   );
	free( mpi.avg_nodestats );
*/

	
	return ok;
}