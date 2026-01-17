/**
 * \file ST_spinup.c
 * \brief Definitions of all [spinup](\ref SPINUP) functions.
 *
 * Contains definitions of all functions related to spinup. The
 * underscore prefix in function names to denote private functions that should
 * NEVER be called outside of this file.
 *
 * \author Chandler Haukap
 * \date September 2019
 * \ingroup SPINUP_PRIVATE
 */

#include <stdlib.h>

// ST_spinup.h contains declarations for runSpinup and
// loadSpinupConditions
#include "ST_spinup.h"
#include "ST_stats.h"
#include "ST_globals.h" // externs UseProgressBar
#include "ST_mortality.h"
#include "sw_src/include/rands.h"
#include "sxw_funcs.h"
#include "sw_src/include/myMemory.h"
#include "sw_src/include/filefuncs.h"
#include "ST_progressBar.h"
#include "ST_grid.h" // externs grid_rng
#include "sw_src/include/SW_Output_outarray.h"
#include "ST_functions.h" // externs `environs_rng`, `resgroups_rng`, `species_rng`


/* =================================================== */
/*                  Global Variables                   */
/* --------------------------------------------------- */
/**
 * \brief Stores the state of the cells following spinup.
 *
 * This information must be stored so it can be loaded in every iteration of
 * the spinup.
 *
 * \ingroup SPINUP
 */
CellType** spinupCells;

/**
 * \brief TRUE if spinup should be run.
 * \ingroup SPINUP
 */
Bool shouldSpinup;

/**
 * \brief TRUE if the program is currently in spinup.
 * \ingroup SPINUP
 */
Bool DuringSpinup;


/******** Local functions. These should all be treated as private. ***********/
static void _run_spinup(void);
static void _beginSpinup(void);
static void _endSpinup(void);
static void _saveAsSpinupConditions(void);


/********* Modular functions defined elsewhere ************/
/* Again, we should clean this up eventually. -CH */

void rgroup_Establish(void);
void rgroup_PartResources(void);
void rgroup_Grow(void);
void rgroup_IncrAges(void);
void parm_Initialize(void);
void Plot_Initialize(void);
void copy_rgroup(const GroupType* src, GroupType* dest);
void copy_species(const SpeciesType* src, SpeciesType* dest);

/**
 * \brief Spin up \ref gridCells.
 *
 * This function takes care of EVERYTHING involved with spinup. After calling
 * this function you can load in the spinup information by calling
 * \ref loadSpinupConditions.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP
 */
void runSpinup(void){
	_beginSpinup();

    /* Dummy accumulators to ensure we do not collect statistics */
	StatType *dummy_Dist, *dummy_Ppt, *dummy_Temp,
  		*dummy_Grp, *dummy_Gsize, *dummy_Gpr, *dummy_Gmort, *dummy_Gestab,
  		*dummy_Spp, *dummy_Indv, *dummy_Smort, *dummy_Sestab, *dummy_Sreceived, *dummy_Grazed;
	FireStatsType *dummy_Gwf;

	/* For iterating over gridCells */
	int i, j;

	/* For iterating over years */
	IntS year;


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
    for (year = 1; year <= SuperGlobals.runSpinupYears; year++)
    {
        if(UseProgressBar){
            logProgress(0, year, SPINUP); // iter = 0 because we are not actually in an iterations loop.
        }
        for (i = 0; i < grid_Rows; ++i)
        { // for each row
            for(j = 0; j < grid_Cols; ++j)
            { // for each column
                // If we should run spinup on this cell
                if(gridCells[i][j].mySpeciesInit.useSpinup){
                    // Load up a cell
                    load_cell(i, j);
                } else {
                    continue; // No spinup requested. Move on to next cell.
                }

                set_all_rngs(SuperGlobals.randseed, 0, year, j * grid_Rows + i);

                /* This step is important. load_cell loaded in the actual
                   accumulators, but we do not want to accumulate stats while
                   in spinup. We need to load in dummy accumulators to ensure
                   we ignore everything that happens in spinup. */
                stat_Copy_Accumulators(dummy_Dist, dummy_Ppt, dummy_Temp,
                                       dummy_Grp, dummy_Gsize, dummy_Gpr,
                                       dummy_Gmort, dummy_Gestab, dummy_Spp,
                                       dummy_Indv, dummy_Smort, dummy_Sestab,
                                       dummy_Sreceived, dummy_Grazed, dummy_Gwf, TRUE);

                Globals->currYear = year;
			    _run_spinup();

            } /* end column */
        } /* end row */
        unload_cell(); // Reset the global variables
    } /* end model run for this year*/

    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    SXW_Reset(gridCells[0][0].mySXW->f_watin, TRUE);
    //TODO: This is a shortcut. swc history is not used and shouldn't be until this is fixed.
    free(SoilWatRun.SoilWatIn.hist.file_prefix);
    SoilWatRun.SoilWatIn.hist.file_prefix = NULL;
    ChDir("..");

	_endSpinup();
}

/**
 * \brief Prepares for spinup by turning on species that have requested spinup
 *        and turning off species that have not requested spinup.
 *
 * This function should be accompanied by a call to \ref _endSpinup.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP_PRIVATE
 */
static void _beginSpinup(void){
	int i, j; 				/* For iterating over cells */
	SppIndex sp;			/* For iterating over species */
	Bool temporary_storage;	/* For swapping variables */

	DuringSpinup = TRUE;

	/* Swap Species[sp]->use_me and mySpeciesInit.shouldSpinup[sp]. shouldSpinup is an array
	   of booleans that represent whether the given species should be used in spinup. use_me is a
	   boolean that represents whether the given species should be used in production. By swaping them we
	   save space, but we have to remember to swap them back before the production run. */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			if(!gridCells[i][j].mySpeciesInit.useSpinup){
				continue;
			}
			load_cell(i, j);	/* We could do this without loading the cell, but there would be no guarantee
									that ForEachGroup would iterate correctly */

			/* Begin swaping variables */
			ForEachSpecies(sp){
				// Temporarily store use_me
				temporary_storage = Species[sp]->use_me;
				// Swap use_me
				Species[sp]->use_me = gridCells[i][j].mySpeciesInit.shouldSpinup[sp];
				// Swap shouldBeSpinup[sp]
				gridCells[i][j].mySpeciesInit.shouldSpinup[sp] = temporary_storage;
			} /* End for each species */
		} /* End for each column */
	} /* End for each row */
	unload_cell(); // Reset the global variables
}

/**
 * \brief Return the program to the state it needs to be in for the main
 *        simulation.
 *
 * \ref gridCells from gridded mode was modified heavily during spinup. Calling
 * this function reverts it to it's original state.
 *
 * This should only be called if you have called \ref _beginSpinup.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP_PRIVATE
 */
static void _endSpinup(void){
	// Calling this function a second time will swap the variables back to their original state.
	_beginSpinup();
	// Save the state of the program as our spinup conditions.
	_saveAsSpinupConditions();
	// We have now exited spinup.
	DuringSpinup = FALSE;
}

/**
 * \brief Save the current state of the program as spinup conditions.
 *
 * This is low level function. If you have already called \ref _endSpinup()
 * there is no need to call this function.
 *
 * \sideeffect
 *     This function will reread the input files.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP_PRIVATE
 */
static void _saveAsSpinupConditions(){
	// Save gridCells as spinupCells
	spinupCells = gridCells;
	// Nullify gridCells. This ensures we can no longer modify spinupCells on accident.
	gridCells = NULL;

	// The easiest way to reallocate gridCells is reread the files.
	rereadInputs();
}

/**
 * \brief Load the state of the program following spinup into \ref gridCells.
 *
 * This function must be called before every iteration of
 * [gridded mode](\ref GRID). Of course, this function will only work after a
 * call to \ref runSpinup.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP
 */
void loadSpinupConditions(){
	int row, col;
	GrpIndex rg;
	SppIndex sp;

	for(row = 0; row < grid_Rows; ++row){
		for(col = 0; col < grid_Cols; ++col){
			load_cell(row, col);

			ForEachSpecies(sp){
				copy_species(spinupCells[row][col].mySpecies[sp], Species[sp]);
			}

			ForEachGroup(rg){
				copy_rgroup(spinupCells[row][col].myGroup[rg], RGroup[rg]);
			}

			copy_environment(&spinupCells[row][col].myEnvironment, Env);
			copy_plot(&spinupCells[row][col].myPlot, Plot);
			copy_succulent(&spinupCells[row][col].mySucculent, Succulent);

			unload_cell();
		}
	}

}

/**
 * \brief "Spinup" the model by running without stat collection, fire, or
 *        grazing.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP_PRIVATE
 */
static void _run_spinup(void)
{
    // We don't want seed dispersal during spinup, so we'll store it here,
    // turn it off, then turn it back on afterwards.
    Bool myUseSeedDispersal = UseSeedDispersal;
    UseSeedDispersal = FALSE;

    // This is a flag from SOILWAT2. It MUST be FALSE during spinup.
    SuperGlobals.prepare_IterationSummary = FALSE;

	Bool killedany;             // killedany for mortality functions

    rgroup_Establish(); 		// Establish individuals. Excludes annuals.
    Env_Generate();				// Generated the SOILWAT environment
    rgroup_PartResources();		// Distribute resources
    rgroup_Grow(); 				// Grow
    mort_Main(&killedany); 		// Mortality that occurs during the growing season
    rgroup_IncrAges(); 			// Increment ages of all plants
    killAnnuals(); 			// Kill annuals
    killMaxage();             // Kill plants that reach max age
    killExtraGrowth(); 		// Kill superfluous growth

    // Turn UseSeedDispersal back on (if it was ever on)
    UseSeedDispersal = myUseSeedDispersal;
}

/**
 * \brief Free memory allocated to spinupCells.
 *
 * This function should only be called once per simulation.
 *
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP
 */
void freeSpinupMemory(void)
{
	// Remember where gridCells pointed.
	CellType** locationOfGridCells = gridCells;

	// If spinupCells is allocated.
	if(spinupCells){
		// move gridCells to point to spinupCells. This allows us to deallocate using free_grid_memory();
		gridCells = spinupCells;
		// Since we already have a function to deallocate gridCells we can use it.
		free_grid_memory();
		// And we need to reset gridCells in case it hasn't been deallocated.
		gridCells = locationOfGridCells;
	}
}
