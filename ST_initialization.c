/***********************************************************************/
// ST_initialization.c 
//      Contains definitions of all functions related to initialization.
//      The current initialization methods are _run_seed_initialization
//      and _run_spinup. Note that _run_seed_initialization is non-functional 
//      as of 12/11/19, but the logic is in place to call the function. This 
//      file uses the "_" prefix in function names to denote private 
//      functions that should NEVER be called outside of this file.
//
// TO ADD AN INITIALIZATION METHOD:
//      Adding a method is simple. runInitialization() takes care of memory 
//      management, iterating, and loading cells. All you need to do in your 
//      new function is define what the program should do for ONE year and 
//      for ONE cell. Once you have your function written, add an entry for 
//      it to the InitializationMethod enumerator in ST_initialization.h, 
//      add your method to _read_grid_setup() of ST_grid.c, then add your 
//      function to the switch statement in runInitialization(). To see
//      an example initialization function check out _run_spinup().
/***********************************************************************/

// ST_initialization.h contains declarations for runInitialization and loadInitializationConditions 
#include "ST_initialization.h" 
#include "ST_grid.h"
#include "ST_stats.h"
#include "ST_globals.h"
#include "ST_defines.h"
#include "sw_src/pcg/pcg_basic.h"
#include "sw_src/rands.h"
#include "sxw_funcs.h"
#include "myMemory.h"
#include "filefuncs.h"
#include "ST_progressBar.h"
#include "ST_stats.h"

/********** Local functions. These should all be treated as private. *************/
static void _run_spinup(void);
static void _run_seed_initialization(void);
static void _beginInitialization(void);
static void _endInitialization(void);
static void _saveAsInitializationConditions(void);

/***************************** Externed variables **********************************/
/* Note that in an ideal world we wouldn't need to extern any variables because 
   every module would declare them in a header file. Hopefully we can get this
   cleaned up soon! -CH */

// these are SOILWAT variables that we need...
extern SW_SOILWAT SW_Soilwat;
extern SW_SITE SW_Site;
extern SW_VEGPROD SW_VegProd;
extern SW_WEATHER SW_Weather;

extern pcg32_random_t grid_rng;         // Gridded mode's unique RNG.

/* We need to seed these RNGs when using the gridded mode but do not use them in this file. */
extern pcg32_random_t environs_rng;     // Used exclusively in ST_environs.c
extern pcg32_random_t mortality_rng;    // Used exclusively in ST_mortality.c
extern pcg32_random_t resgroups_rng;    // Used exclusively in ST_resgroups.c
extern pcg32_random_t species_rng;      // Used exclusively in ST_species.c
extern pcg32_random_t markov_rng;       // Used exclusively in SW_Markov.c

extern Bool UseProgressBar;             // From ST_main.c
extern Bool* _SomeKillage;              // From ST_mortality.c

/********* Modular functions defined elsewhere ************/
/* Again, we should clean this up eventually. -CH */

void rgroup_Establish(void);
void rgroup_PartResources(void);
void rgroup_Grow(void); 
void mort_Main(Bool* killed);
void rgroup_IncrAges(void);
void _kill_annuals(void);
void _kill_maxage(void);
void proportion_Recovery(void); 
void _kill_extra_growth(void); 	
void parm_Initialize(void);
void Plot_Initialize(void);
void copy_rgroup(const GroupType* src, GroupType* dest);
void copy_species(const SpeciesType* src, SpeciesType* dest);

/* Initializes the plot with whichever method you have specified with initializationMethod. 
   This function takes care of EVERYTHING involved with initialization.
   After calling this function you can load in the initialization information by calling 
   loadInitializationConditions(). */
void runInitialization(void){
	_beginInitialization();

    /* Dummy accumulators to ensure we do not collect statistics */
	StatType *dummy_Dist, *dummy_Ppt, *dummy_Temp,
  		*dummy_Grp, *dummy_Gsize, *dummy_Gpr, *dummy_Gmort, *dummy_Gestab,
  		*dummy_Spp, *dummy_Indv, *dummy_Smort, *dummy_Sestab, *dummy_Sreceived;
	FireStatsType *dummy_Gwf;

	/* For iterating over gridCells */
	int i, j;

	/* For iterating over years */
	IntS year;

    /* Initialization is technically an iteration so we need to seed the RNGs. */
    RandSeed(SuperGlobals.randseed, &environs_rng);
    RandSeed(SuperGlobals.randseed, &mortality_rng);
    RandSeed(SuperGlobals.randseed, &resgroups_rng);
    RandSeed(SuperGlobals.randseed, &species_rng);
    RandSeed(SuperGlobals.randseed, &grid_rng);
    RandSeed(SuperGlobals.randseed, &markov_rng);

    // Initialize the plot for each grid cell
    for (i = 0; i < grid_Rows; i++){
        for (j = 0; j < grid_Cols; j++){
            load_cell(i, j);
            Plot_Initialize();
            Globals->currIter = 1; // Iteration doesn't really matter, I set it to 1 here just in case.
        }
    }
    unload_cell(); // Reset the global variables

    /* Iterate through the number of years requested in inputs. */
    for (year = 1; year <= SuperGlobals.runInitializationYears; year++)
    {
        if(UseProgressBar){
            logProgress(0, year, INITIALIZATION); // iter = 0 because we are not actually in an iterations loop.
        }
        for (i = 0; i < grid_Rows; ++i)
        { // for each row
            for(j = 0; j < grid_Cols; ++j)
            { // for each column
                // If we should run spinup on this cell
                if(gridCells[i][j].mySpeciesInit.useInitialization){
                    // Load up a cell
                    load_cell(i, j);
                } else {
                    continue; // No spinup requested. Move on to next cell.
                }

                /* This step is important. load_cell loaded in the actual accumulators, but we do not want
                    to accumulate stats while in spinup. We need to load in dummy accumulators to ensure
                    we ignore everything that happens in spinup. */
                stat_Copy_Accumulators(dummy_Dist, dummy_Ppt, dummy_Temp, dummy_Grp, dummy_Gsize, dummy_Gpr, dummy_Gmort, dummy_Gestab,
                                        dummy_Spp, dummy_Indv, dummy_Smort, dummy_Sestab, dummy_Sreceived, dummy_Gwf, TRUE);

                Globals->currYear = year;

                switch (initializationMethod){
		            case INIT_WITH_SPINUP:
			            _run_spinup();
			            break;
		            case INIT_WITH_SEEDS:
			            _run_seed_initialization();
			            break;
		            default:
			            break;
	            }

            } /* end column */
        } /* end row */
        unload_cell(); // Reset the global variables
    } /* end model run for this year*/

    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    SXW_Reset(gridCells[0][0].mySXW->f_watin);
    //TODO: This is a shortcut. swc history is not used and shouldn't be until this is fixed.
    Mem_Free(SW_Soilwat.hist.file_prefix);
    SW_Soilwat.hist.file_prefix = NULL;
    ChDir("..");

	_endInitialization();
}

/* Prepares for initialization by turning on species that have requested initialization and turning off 
   species that have not requested initialization. This function should be accompanied by a call to 
   _endInitialization. */
static void _beginInitialization(void){
	int i, j; 				/* For iterating over cells */
	SppIndex sp;			/* For iterating over species */
	Bool temporary_storage;	/* For swapping variables */

	DuringInitialization = TRUE;

	/* Swap Species[sp]->use_me and mySpeciesInit.shouldBeInitialized[sp]. shouldBeInitialized is an array 
	   of booleans that represent whether the given species should be used in initialization. use_me is a 
	   boolean that represents whether the given species should be used in production. By swaping them we 
	   save space, but we have to remember to swap them back before the production run. */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			if(!gridCells[i][j].mySpeciesInit.useInitialization){
				continue;
			}
			load_cell(i, j);	/* We could do this without loading the cell, but there would be no guarantee
									that ForEachGroup would iterate correctly */

			/* Begin swaping variables */
			ForEachSpecies(sp){
				// Temporarily store use_me
				temporary_storage = Species[sp]->use_me;
				// Swap use_me
				Species[sp]->use_me = gridCells[i][j].mySpeciesInit.shouldBeInitialized[sp]; 
				// Swap shouldBeInitialized[sp]
				gridCells[i][j].mySpeciesInit.shouldBeInitialized[sp] = temporary_storage;
			} /* End for each species */
		} /* End for each column */
	} /* End for each row */
	unload_cell(); // Reset the global variables
}

/* Return the program to the state it needs to be in for the main simulation. This should only be called if you
   have called _beginInitialization. */
static void _endInitialization(void){
	// Calling this function a second time will swap the variables back to their original state.
	_beginInitialization();
	// Save the state of the program as our initialization conditions.
	_saveAsInitializationConditions();
	// We have now exited initialization.
	DuringInitialization = FALSE;
}

/* Save the current state of the program as spinup conditions.
 *  
 * This is low level function. If you have already called
 * _endInitialization() there is no need to call this function. */
static void _saveAsInitializationConditions(){
	// Save gridCells as initializationCells
	initializationCells = gridCells;
	// Nullify gridCells. This ensures we can no longer modify initializationCells on accident.
	gridCells = NULL;

	// The easiest way to reallocate gridCells is reread the files.
	rereadInputs();
}

/* Load the state of the program right after initialization. */
void loadInitializationConditions(){
	int row, col;
	GrpIndex rg;
	SppIndex sp;

	for(row = 0; row < grid_Rows; ++row){
		for(col = 0; col < grid_Cols; ++col){
			load_cell(row, col);

			ForEachSpecies(sp){
				copy_species(initializationCells[row][col].mySpecies[sp], Species[sp]);
			}

			ForEachGroup(rg){
				copy_rgroup(initializationCells[row][col].myGroup[rg], RGroup[rg]);
			}

			copy_environment(&initializationCells[row][col].myEnvironment, Env);
			copy_plot(&initializationCells[row][col].myPlot, Plot);
			copy_succulent(&initializationCells[row][col].mySucculent, Succulent);

			unload_cell();
		}
	}

}

/* "Spinup" the model by running without stat collection, fire, or grazing.
 * 
 * Fire and grazing will potentially be added to spinup as a future feature. */
static void _run_spinup(void)
{
    // We don't want seed dispersal durring spinup, so we'll store it here,
    // turn it off, then turn it back on afterwards.
    Bool myUseSeedDispersal = UseSeedDispersal;
    UseSeedDispersal = FALSE;

	Bool killedany;             // killedany for mortality functions

    rgroup_Establish(); 		// Establish individuals. Excludes annuals.
    Env_Generate();				// Generated the SOILWAT environment
    rgroup_PartResources();		// Distribute resources
    rgroup_Grow(); 				// Grow
    mort_Main(&killedany); 		// Mortality that occurs during the growing season
    rgroup_IncrAges(); 			// Increment ages of all plants
    _kill_annuals(); 			// Kill annuals
    _kill_maxage();             // Kill plants that reach max age
    _kill_extra_growth(); 		// Kill superfluous growth	

    // Turn UseSeedDispersal back on (if it was ever on)
    UseSeedDispersal = myUseSeedDispersal;		
}

/* TODO: This is a dummy method. It needs to be implemented once seed dispersal is fully planned.
 *
 * Note on implementing: The logic behind calling this function already exists. You can view the 
 *     function call inside runInitialization, but the important thing to note is that this 
 *     function is called inside the following structure:
 *     
 *     for each year
 *         for each grid cell
 *             _run_seed_initialization
 * 
 *     You can assume that prior to this function being called the correct grid cell has been loaded
 *     into the global variables (RGroup, Species, Env, ect.) and that this function must
 *     operate on a yearly time step. */
static void _run_seed_initialization(void){
    if(Globals->currYear == 1){
	    printf("\nYou have attempted to initialize with seed dispersal.\n" 
	            "This option is currently in development and will be availible soon.\n");
    }
    initializationMethod = INIT_WITH_NOTHING;
	return;
}

/* Free memory allocated to initializationCells. This function should only be called once per simulation. */
void freeInitializationMemory(void)
{
	// Remember where gridCells pointed.
	CellType** locationOfGridCells = gridCells;

	// If initializationCells is allocated.
	if(initializationCells){
		// move gridCells to point to initializationCells. This allows us to deallocate using free_grid_memory();
		gridCells = initializationCells;
		// Since we already have a function to deallocate gridCells we can use it.
		free_grid_memory();
		// And we need to reset gridCells in case it hasn't been deallocated.
		gridCells = locationOfGridCells;
	}
}