// @2007, Jonathan Belof
// Space Research Group
// Department of Chemistry
// University of South Florida


#include "Atom.h"
#include "Molecule.h"
#include "Pair.h"
#include "SafeOps.h"
#include "System.h"



// Allocate the pair lists 
void System::allocate_pair_lists() {

	Pair *pair_ptr, *prev_pair_ptr;

	natoms = countNatoms();

	// build atom and molecule arrays
	rebuild_arrays();
	int n = natoms;

	// setup the pairs, lower triangular
	for(int i = 0; i < (n - 1); i++) {

		SafeOps::calloc( atom_array[i]->pairs, 1, sizeof(Pair), __LINE__, __FILE__ );
		pair_ptr = atom_array[i]->pairs;
		prev_pair_ptr = pair_ptr;

		for(int j = (i + 1); j < n; j++) {
			SafeOps::calloc(pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__ );
			prev_pair_ptr = pair_ptr;
			pair_ptr = pair_ptr->next;
		}

		prev_pair_ptr->next = nullptr;
		free(pair_ptr);
	}
}




// add new pairs for when a new molecule is created */
void System::update_pairs_insert() {

	int         n = 0;
	Molecule  * molecule_ptr;
	Atom      * atom_ptr;
	Pair      * pair_ptr;


	// count the number of atoms per molecule
	for( atom_ptr = checkpoint->molecule_altered->atoms; atom_ptr; atom_ptr = atom_ptr->next, n++);

	// add n number of pairs to altered and all molecules ahead of it in the list
	for(molecule_ptr = molecules; molecule_ptr != checkpoint->tail; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			// go to the end of the pair list
			if(atom_ptr->pairs) {

				// go to the end of the pair list
				for(pair_ptr = atom_ptr->pairs; pair_ptr->next; pair_ptr = pair_ptr->next);

				// tag on the extra pairs
				for(int i = 0; i < n; i++) {
					SafeOps::calloc( pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__);
					pair_ptr = pair_ptr->next;
				}

			} else {

				// needs a new list
				SafeOps::calloc( atom_ptr->pairs, 1, sizeof(Pair), __LINE__, __FILE__ );
				pair_ptr = atom_ptr->pairs;
				for( int i = 0; i < (n - 1); i++) {
					SafeOps::calloc( pair_ptr->next, 1, sizeof(Pair), __LINE__, __FILE__ );
					pair_ptr = pair_ptr->next;
				}

			}

		} // for atom 
	} // for molecule
}




// remove pairs when a molecule is deleted
void System::update_pairs_remove() {

	int         n            = 0,
		        m            = 0;
	Molecule  * molecule_ptr = nullptr;
	Atom      * atom_ptr     = nullptr;
	Pair      * pair_ptr     = nullptr,
	         ** pair_array   = nullptr;


	// count the number of atoms per molecule
	for( atom_ptr = checkpoint->molecule_backup->atoms, n=0; atom_ptr; atom_ptr = atom_ptr->next, n++);

	// remove n number of pairs for all molecules ahead of the removal point
	SafeOps::calloc( pair_array, 1, sizeof(Pair *), __LINE__, __FILE__ );
	for( molecule_ptr = molecules; molecule_ptr != checkpoint->tail; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			// build the pair pointer array
			for( pair_ptr = atom_ptr->pairs, m = 0; pair_ptr; pair_ptr = pair_ptr->next, m++) {

				SafeOps::realloc( pair_array, sizeof(Pair *)*(m + 1), __LINE__, __FILE__ );
				pair_array[m] = pair_ptr;

			}

			for(int i = (m - n); i < m; i++)
				free(pair_array[i]);

			// handle the end of the list
			if((m - n) > 0)
				pair_array[(m - n - 1)]->next = nullptr;
			else
				atom_ptr->pairs = nullptr;

		}
	}

	// free our temporary array
	free(pair_array);

}