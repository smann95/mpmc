// Space Research Group
// Department of Chemistry
// University of South Florida
#define VERSION 0


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constants.h"
#include "Output.h"
#include "SimulationControl.h"
#ifdef MPI
	#include <mpi.h>
#endif
int rank = 0;  
int size = 1;


void die(int code);  // Kill any MPI Processes and stop execution
void usage(char *progname); // Print usage notes, then stop execution


// Interrupt signal handler 
#ifdef __unix__
	#include <unistd.h>
	#include <signal.h>
	SimulationControl *sc;
	// on SIGTERM, cleanup and exit 
	void signal_handler( int sigtype ) {
		Output::out(" ************ SIGTERM, SIGUSR1 or SIGUSR2 received, exiting *************\n");
		//sc->close_files();
		//sc->cleanup();
		#ifdef MPI
		//	if(!rank)
		//		sc->close_files();
		#else
			die( EXIT_FAILURE );
		#endif // MPI
	}
#endif








int main(int argc, char * argv[]) {

	time_t t = time(nullptr);
	struct tm tm = *localtime(&t);
	char linebuf[maxLine] = {'\0'};

	
	#ifdef MPI
		// Start up the MPI chain 
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
	#endif // MPI



	if (!rank) {
		sprintf(linebuf, "MPMC (Massively Parallel Monte Carlo) r%d++ - 2012-2016 GNU Public License\n", VERSION);
		Output::out(linebuf);
		sprintf(linebuf, "MAIN: processes started on %d cores @ %d-%d-%d %d:%d:%d\n", size, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		Output::out(linebuf);
		// Check for validity of command line argument count
		if (argc != 2) 
			usage(argv[0]);
	}

	
	
	try {

		SimulationControl simController(argv[1]);  // create the system
		Output::out( "MAIN: Simulation parameters established.\n" );

		simController.initializeSimulationObjects();
		Output::out( "MAIN: System(s) data structures allocated and initialized.\n" );
		
		#if defined( _POSIX_VERSION )
			// install the signal handler to catch SIGTERM cleanly
			signal(SIGTERM, signal_handler ); signal(SIGUSR1, signal_handler );	signal(SIGUSR2, signal_handler );
			sc = &simController;     ////////////////////////////////////////////////////////////////////////////////////////////////////////////// ADD SUPPORT FOR CLEAN UP AT THE TOP OF THIS FILE
		#endif
		
		simController.runSimulation();   

	}
	catch (int e) { 
		sprintf(linebuf, "Program exiting with MPMC error code %d\n", e);
		Output::err(linebuf);
		die(EXIT_FAILURE); 
	}
	
	die(EXIT_SUCCESS);
}








// Print usage info and quit
void usage(char *progname) {
	if (!rank) {
		fprintf( stderr, "usage: %s <config>\n", progname);
		fprintf( stderr, "See: https://github.com/mpmccode/mpmc\n");
	}
	die(EXIT_FAILURE);
}



// Kill MPI before quitting, when neccessary
void die(int code) {
	#ifdef MPI
		MPI_Finalize();
	#endif
	exit(code);
}