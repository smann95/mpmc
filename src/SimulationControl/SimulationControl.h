#pragma once

#include <stdint.h>
#include <stdio.h>
#include <vector>

#include "constants.h"
#include "System.h"




class SimulationControl
{

public:
	SimulationControl( char * inFilename );
	~SimulationControl();
	bool runSimulation();
	void initializeSimulationObjects();

private:
	
	System   sys;
	std::vector<System*> systems;

	
	
	
	
	

private:
	
	

	void read_config(char *inFilename );
	// Parses a simulation input file and populates the System Controller with
	// the data found therein. Data read in is largely unvalidated.

	bool process_command( char token[maxTokens][maxLine] );
	// Read each input line from the input file and set individual system flags accordingly.
	// Minimal error checking performed, mostly in the form of checking for malformed
	// command sequences. 

	// The following functions validate the system, checking for missing, invalid or incompatible options & settings.
	bool check_system();
	bool check_mc_options();
	bool check_spectre_options();
	bool check_io_files_options();
	bool check_feynman_hibbs_options();
	bool check_simulated_annealing_options();
	bool check_hist_options();
	bool check_polarization_options();
	bool check_qrot_options();
	double get_rand();

	// Allocate data structures and initialize system data for multi-system ensembles
	void initialize_PI_NVT_Systems();
	void setup_Gibbs_Systems();

	// SimulationControl.PathIntegral.cpp
	bool pi_nvt_mc();
	bool check_pi_nvt_options();

	// SimulationControl.Gibbs.cpp
	bool Gibbs_mc();
	bool check_Gibbs_options();
	void setup_Gibbs_systems();

};


