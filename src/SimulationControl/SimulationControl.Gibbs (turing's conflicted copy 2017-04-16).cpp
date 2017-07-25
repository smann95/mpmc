

#include "SafeOps.h"
#include "SimulationControl.h"
#include "Output.h"

extern int rank, size;


bool SimulationControl::check_Gibbs_options() {

	// Right now, Gibbs is going through the normal check_mc_options() validity check.
	// Place checks for any additional constraints here.
	// If constraints are violated:
	// return fail;
	
	// Checks on probabilities are done in setup_Gibbs_Systems()---if volume_probability is
	// not explicitly set, it is set to 1/N and the system is not set up (to be able to
	// compute N) until setup_Gibbs_Systems()

	return ok;
}



void SimulationControl::setup_Gibbs_Systems() {

	char linebuf[maxLine] = {0};

	// THIS WILL SET UP TWO SYSTEMS, systems[0] AND systems[1] IN THE USUAL MPMC METHOD
	System **temp_sys;
	SafeOps::calloc(temp_sys, 2, sizeof(System*), __LINE__, __FILE__);
	temp_sys[0] = new System(sys);
	temp_sys[1] = new System(sys);
	systems.push_back(temp_sys[0]);
	systems.push_back(temp_sys[1]);
	
	for( int i=0; i<2; i++ ) {
		
		systems[i]->setup_simulation_box();
			
		sprintf( linebuf, "SIM_CONTROL, SYS %d: simulation box configured.\n", i );
		Output::out(linebuf);
		
		systems[i]->allocateStatisticsMem();

		systems[i]->allocate_pair_lists();
		sprintf( linebuf, "SIM_CONTROL, SYS %d: finished allocating pair lists\n", i );
		Output::out(linebuf);
				
		systems[i]->pairs();          // get all of the pairwise interactions, exclusions, etc.
		systems[i]->flag_all_pairs(); // set all pairs to initially have their energies calculated 
		sprintf( linebuf, "SIM_CONTROL, SYS %d: finished calculating pairwise interactions\n", i );
		Output::out(linebuf);

		if( systems[i]->cavity_bias )
			systems[i]->setup_cavity_grid();

		// if polarization active, allocate the necessary matrices 
		if(  systems[i]->polarization   &&   !systems[i]->cuda   &&   !systems[i]->polar_zodid  )
			systems[i]->thole_resize_matrices();

		// if histogram calculation flag is set, allocate grid
		if( systems[i]->calc_hist ) {
			systems[i]->setup_histogram();
		}

		// Seed RNG?

		if(  ! systems[i]->use_sg   ||   systems[i]->rd_only  )
		{
			char linebuf[maxLine] = {0};
			sprintf(linebuf, "SIM_CONTROL, SYS %d: Ewald gaussian width = %f A\n", i, systems[i]->ewald_alpha);
			Output::out(linebuf);
			sprintf(linebuf, "SIM_CONTROL, SYS %d: Ewald kmax = %d\n", i, systems[i]->ewald_kmax);
			Output::out(linebuf);
		}

		// ensure that all SPECTRE charges lie within the restricted domain 
		if( systems[i]->spectre )
			systems[i]->spectre_wrapall();
	}
	

	// Set the probability of the volume move to 1/N if it wasn't explicitly set
	if( sys.volume_probability == 0.0 ) {
		sys.volume_probability =
		systems[0]->volume_probability =
		systems[1]->volume_probability = 1.0/(double)( systems[0]->observables->N + systems[1]->observables->N );
	}
	sprintf(linebuf, "SIM_CONTROL: volume change probability is %lf.\n", systems[0]->volume_probability);
	Output::out(linebuf);


	// If quantum rotation is turned off, make spinflip probability 0.0, just in case
	if(sys.quantum_rotation) {
		sprintf(linebuf, "SIM_CONTROL:      spinflip probability is %lf.\n", sys.spinflip_probability);
		Output::out(linebuf);
	} else
		sys.spinflip_probability = 0.0;


	if(sys.transfer_probability == 0.0) {
		Output::err( "(ERROR) SIM_CONTROL: transfer move probability was either not set, or set to 0.0 in a Gibbs NVT simulation. Set with keyword \"transfer_probability\" in input file." );
		throw missing_setting;
	} else {
		sprintf( linebuf, "SIM_CONTROL:      transfer probability is %lf.\n", sys.transfer_probability );
		Output::out(linebuf);
	}

	sprintf(linebuf, "SIM_CONTROL:      displace probability is %lf.\n", 1.0 - (sys.spinflip_probability + sys.transfer_probability + sys.volume_probability) );
	Output::out(linebuf);

	// Make sure probability of all moves adds up to 1.0 (displacement is the
	// default move so all the other moves must some to a probability < 1.0
	if (sys.spinflip_probability + sys.volume_probability + sys.transfer_probability >= 1.0) {
		// To allow displacement probabilities of 0.0, change the above check to   ... >    1.0
		Output::err("(ERROR) SIM_CONTROL: Invalid probabilities set. The summed frequencies for spinflip, volume, transfer, and displacement moves may not exceed 1.0.");
	}
}



bool SimulationControl::Gibbs_mc()  {

	double  initial_energy[2]  = {0}, 
	        final_energy  [2]  = {0},
	        current_energy[2]  = {0};
	System::mpiData mpi   [2];
	
	for( int i=0; i<2; i++ ) {
		if (systems[i]->sorbateCount > 1)
			SafeOps::calloc(systems[i]->sorbateGlobal, systems[i]->sorbateCount, sizeof(System::sorbateAverages_t), __LINE__, __FILE__);
		if( systems[i]->cavity_bias      )
			systems[i]->cavity_update_grid(); // update the grid for the first time 
		systems[i]->observables->volume = systems[i]->pbc.volume; // set volume observable
		initial_energy[i] = systems[i]->mc_initial_energy();
		mpi[i] = systems[i]->setup_mpi();
	}

	// save the initial state 
	System::do_checkpoint_NVT_Gibbs( systems );
	// backup observables
	// pick a move
	//   -change move if impossible
	// pick an atom
	//   -make arrays
	//   -mark molecule targets 
	//   -clean up


	int s=1; 
	int max_step = systems[0]->numsteps;

	// main MC loop 
	for( s = 1; s <= max_step; s++ ) {

		systems[0]->step = systems[1]->step = s;
		
		// restore the last accepted energy 
		initial_energy[0] = systems[0]->observables->energy;
		initial_energy[1] = systems[1]->observables->energy;
		
		// perturb the system 
		// systems[0]->make_move();
		// systems[1]->make_move();


		// calculate the energy change 
		final_energy[0] = systems[0]->energy();
		final_energy[1] = systems[1]->energy();

		#ifdef QM_ROTATION
			// solve for the rotational energy levels 
			if(system->quantum_rotation && (system->checkpoint->movetype == MOVETYPE_SPINFLIP))
				quantum_system_rotational_energies(system);
		#endif // QM_ROTATION 

		// treat a bad contact as a reject 
		if( !std::isfinite(final_energy[0]) ) {
			systems[0]->observables->energy           = MAXVALUE;
			systems[0]->nodestats  ->boltzmann_factor = 0;
		} 
		else
			System::boltzmann_factor_NVT_Gibbs( *systems[0], initial_energy[0], final_energy[0], *systems[0], initial_energy[0], final_energy[0]);

		if( !std::isfinite(final_energy[1]) ) {
			systems[1]->observables->energy           = MAXVALUE;
			systems[1]->nodestats  ->boltzmann_factor = 0;
		} 
		else
			systems[1]->boltzmann_factor( initial_energy[1], final_energy[1] );


		// Metropolis function (System 0)
		if(  (get_rand() < systems[0]->nodestats->boltzmann_factor)   &&   (systems[0]->iter_success == 0)  ) {
			/////////// ACCEPT

			current_energy[0] = final_energy[0];

			// checkpoint 
			systems[0]->do_checkpoint();
			systems[0]->register_accept();

			// SA 
			if( systems[0]->simulated_annealing )
			{
				if( systems[0]->simulated_annealing_linear == 1 )
				{
					systems[0]->temperature = systems[0]->temperature + (systems[0]->simulated_annealing_target - systems[0]->temperature)/(max_step - s);
					if(  (max_step - s)   ==   0   )
						systems[0]->temperature = systems[0]->simulated_annealing_target;
				}
				else
					systems[0]->temperature = systems[0]->simulated_annealing_target + (systems[0]->temperature - systems[0]->simulated_annealing_target)*systems[0]->simulated_annealing_schedule;
			}

		} else { 
			// REJECT ///////////////////////////////////////////////////////////////
			//reset the polar iterative failure flag and restore from last checkpoint
			current_energy[0] = initial_energy[0]; //used in parallel tempering
			systems[0]->iter_success = 0;
			systems[0]->restore();
			systems[0]->register_reject();
		}


		// Metropolis function (System 1)
		if(  (get_rand() < systems[1]->nodestats->boltzmann_factor)   &&   (systems[1]->iter_success == 0)  ) {
			/////////// ACCEPT

			current_energy[1] = final_energy[1];

			// checkpoint 
			systems[1]->do_checkpoint();
			systems[1]->register_accept();

			// SA 
			if( systems[1]->simulated_annealing )
			{
				if( systems[1]->simulated_annealing_linear == 1 )
				{
					systems[1]->temperature = systems[1]->temperature + (systems[1]->simulated_annealing_target - systems[1]->temperature)/(max_step - s);
					if(  (max_step - s)   ==   0   )
						systems[1]->temperature = systems[1]->simulated_annealing_target;
				}
				else
					systems[1]->temperature = systems[1]->simulated_annealing_target + (systems[1]->temperature - systems[1]->simulated_annealing_target)*systems[1]->simulated_annealing_schedule;
			}

		} else { 
			// REJECT ///////////////////////////////////////////////////////////////
			//reset the polar iterative failure flag and restore from last checkpoint
			current_energy[1] = initial_energy[1]; //used in parallel tempering
			systems[1]->iter_success = 0;
			systems[1]->restore();
			systems[1]->register_reject();
		}



		// perform parallel_tempering
		if ( (systems[0]->parallel_tempering) && ((s % systems[0]->ptemp_freq) == 0) )
			systems[0]->temper_system( current_energy[0] );
		if ( (systems[1]->parallel_tempering) && ((s % systems[1]->ptemp_freq) == 0) )
			systems[1]->temper_system( current_energy[1] );

		// track the acceptance_rate 
		systems[0]->track_ar( systems[0]->nodestats );
		systems[1]->track_ar( systems[1]->nodestats );


		// each node calculates its stats 
		systems[0]->update_nodestats( systems[0]->nodestats, systems[0]->avg_nodestats );
		systems[1]->update_nodestats( systems[1]->nodestats, systems[1]->avg_nodestats );


		// do this every correlation time, and at the very end 
		if(  !(s % systems[0]->corrtime)  ||  (s == max_step)  )
			systems[0]->do_corrtime_bookkeeping(mpi[0]);
		if(  !(s % systems[1]->corrtime)  ||  (s == max_step)  )
			systems[1]->do_corrtime_bookkeeping(mpi[1]);
			

	} // main loop 
	


	// output restart files for each node
	for( int i=0; i<2; i++ ) {
		if(  systems[i]->write_molecules_wrapper( systems[i]->pqr_output ) < 0  ) {
			Output::err( "MC: could not write final state to disk\n" );
			throw unknown_file_error;
		}

		if(!rank) {
			free( mpi[i].rcv_strct   );
			free( mpi[i].temperature );
		}

		if( systems[i]->sorbateCount > 1 )
			free( mpi[i].sinfo );
	
		free( mpi[i].snd_strct     );
		free( mpi[i].observables   );
		free( mpi[i].avg_nodestats );
	}

	return ok;
}