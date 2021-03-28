/** 
 * \file ST_grid.c
 * \brief Function definitions for the gridded mode.
 *
 * This module handles the gridded mode of STEPWAT2. To accomplish this we use
 * a grid of cells represented by the CellType struct. The entire grid of
 * cells can be referenced by the gridCells variable which is a 2d array of
 * CellTypes. To allow this module to use the same functions as non-gridded
 * mode the CellType structs must be loaded into the global variables using
 * the load_cell function. As long as a cell is loaded in you can be sure that
 * all functions will work as expected.
 *
 * In addition to all of the functionality of non-gridded mode, gridded mode
 * has three additional features: [spinup](\ref SPINUP), 
 * [seed dispersal](\ref SEED_DISPERSAL) and [colonization](\ref COLONIZATION).
 * Spinup allows vegetation to establish before the simulation
 * experiments begin. Seed dispersal allows each cell to disperse seeds to
 * nearby cells.
 * 
 * See issue #262 and pull request #375 on GitHub for a discussion of the 
 * overhaul of this module.
 * 
 * \author DLM (initial programming)
 * \author Fredrick Pierson
 * \author Chandler Haukap
 * \date March through July 2019
 * \ingroup GRID_PRIVATE
 */

/*******************************************************/
/* -------------- INCLUDES / DEFINES ----------------- */
/*******************************************************/
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sw_src/rands.h"
#include "ST_grid.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sxw_funcs.h"
#include "ST_spinup.h"
#include "ST_progressBar.h"
#include "ST_colonization.h"
#include "sw_src/SW_Output.h"
#include "sw_src/SW_Output_outtext.h"
#include "sw_src/SW_Output_outarray.h"
#include "ST_mortality.h"

int grid_Cells;
Bool UseDisturbances, UseSoils;

/***************************** Externed variables **********************************/
/* Note that in an ideal world we wouldn't need to extern any variables because 
   every module would declare them in a header file. Hopefully we can get this
   cleaned up soon! -CH */

// these are SOILWAT variables that we need...
extern SW_SOILWAT SW_Soilwat;
extern SW_WEATHER SW_Weather;

extern pcg32_random_t grid_rng;         // Gridded mode's unique RNG.

/* We need to seed these RNGs when using the gridded mode but do not use them in this file. */
extern pcg32_random_t environs_rng;     // Used exclusively in ST_environs.c
extern pcg32_random_t resgroups_rng;    // Used exclusively in ST_resgroups.c
extern pcg32_random_t species_rng;      // Used exclusively in ST_species.c
extern pcg32_random_t markov_rng;       // Used exclusively in SW_Markov.c

extern Bool UseProgressBar;             // From ST_main.c
extern Bool *_SomeKillage;				// From ST_mortality.c 

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
//from ST_species.c
void save_annual_species_relsize(void);
void copy_species(const SpeciesType* src, SpeciesType* dest);

//from ST_resgroups.c
void rgroup_Grow(void);
void rgroup_Establish(void);
void rgroup_IncrAges(void);
void rgroup_PartResources(void);
void copy_rgroup(const GroupType* src, GroupType* dest);

//functions from ST_params.c
void parm_Initialize(void);
void parm_SetFirstName(char *s);
void parm_SetName(char *s, int which);
void parm_free_memory(void);
void files_init(void);
void maxrgroupspecies_init(void);

//from ST_main.c
void Plot_Initialize(void);
void deallocate_Globals(Bool isGriddedMode);

/* Functions from sxw.c */
SXW_t* getSXW(void);
SXW_resourceType* getSXWResources(void);
transp_t* getTranspWindow(void);
void copy_sxw_variables(SXW_t* newSXW, SXW_resourceType* newSXWResources, transp_t* newTransp_window);

extern Bool prepare_IterationSummary; // defined in `SOILWAT2/SW_Output.c`
extern Bool print_IterationSummary;   // defined in `SOILWAT2/SW_Output_outtext.c`
extern Bool storeAllIterations; // defined in `SOILWAT2/SW_Output.c`

/***********************************************************/
/* --------- Locally Used Function Declarations ---------- */
/***********************************************************/
void _copy_soils(SoilType* src, SoilType* dest);
static void printGeneralInfo(void);
static void _init_grid_files(void);
static void _init_SXW_inputs(Bool init_SW, char *f_roots);
static void _allocate_gridCells(int rows, int cols);
static void _Output_AllCellAvgBmass(const char* filename);
static void _Output_AllCellAvgMort(const char* filename);
static void _allocate_accumulators(void);
static void _read_disturbances_in(void);
static void _read_soils_in(void);
static int _read_soil_line(char* buf, SoilType* destination, int layer);
static void _read_spinup_species(void);
static void _read_maxrgroupspecies(void);
static void _read_grid_setup(void);
static void _read_files(void);
static void _init_soilwat_outputs(char* fileName);
static void _init_stepwat_inputs(void);
static void _init_grid_inputs(void);
static void _separateSOILWAT2Output(void);
static void _separateSOILWAT2DailyOutput(char* fileName, int* cellNumbers);
static void _separateSOILWAT2MonthlyOutput(char* fileName, int* cellNumbers);
static void _separateSOILWAT2YearlyOutput(char* fileName, int* cellNumbers);
static int _getNumberSOILWAT2OutputCells(void);
static int* _getSOILWAT2OutputCells(void);

/******************** Begin Model Code *********************/
/***********************************************************/

/**
 * \brief Print information about the simulation to stdout.
 * 
 * \author Chandler Haukap
 * \date August 2019 
 * \ingroup GRID 
 */
static void printGeneralInfo(void){
	/* ------------------- Print some general information to stdout ----------------------- */
    printf("Number of iterations: %d\n", SuperGlobals.runModelIterations);
    printf("Number of years: %d\n", SuperGlobals.runModelYears);
	printf("Number of cells: %d\n\n", grid_Cells);
	if(UseDisturbances) printf("Using grid disturbances file\n");
	if(UseSoils) printf("Using grid soils file\n");
	if(shouldSpinup){
		printf("Number of spinup years: %d\n", SuperGlobals.runSpinupYears);
	}
	if(UseSeedDispersal){
		printf("Dispersing seeds between cells\n");
	}
	printf("\n");
	/* --------------------------- END printing general info -------------------------------- */
}

/**
 * \brief Run gridded mode.
 * 
 * This function is responsible for initializing the [cells](\ref CellType), 
 * running the simulation on all of them, then printing output.
 * 
 * \sideeffect
 *     In theory this function will have no side effects, but memory leaks are
 *     possible due to the massive amount of allocations associated with
 *     gridded mode.
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup GRID
 */
void runGrid(void)
{
	int i, j;
	Bool killedany;
	IntS year, iter;

  if(storeAllIterations) {
    printf("WARNING: SOILWAT2 iteration output unavailable for gridded mode.\n");
    storeAllIterations = FALSE;
  }

	_init_grid_files();				// reads in files.in file
	_read_maxrgroupspecies();       // reads in maxrgroupspecies.in file
	_read_grid_setup();             // reads in grid_setup.in file
    _read_files();                  // reads in Stepwat_Inputs/files.in file
    _init_stepwat_inputs();			// reads the stepwat inputs in
	_init_grid_inputs();			// reads the grid inputs in & initializes the global grid variables
	initColonization(grid_files[GRID_FILE_COLONIZATION]);
	//SWC hist file prefix needs to be cleared
	Mem_Free(SW_Soilwat.hist.file_prefix);
	SW_Soilwat.hist.file_prefix = NULL;

	printGeneralInfo();

	if(shouldSpinup){
		runSpinup();
	} else {
		/* SXW expects to be run from the testing.sagebrush.master/Stepwat_Inputs directory.
        However, we are running from the testing.sagebrush.master directory. To find the 
        location of the SOILWAT input files we need to manually set SXW->f_watin. */
		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
		SXW_Reset(gridCells[0][0].mySXW->f_watin);
		ChDir("..");
	}

	// SOILWAT resets SW_Weather.name_prefix every iteration. This is not the behavior we want 
	// so the name is stored here.
	char SW_prefix_permanent[2048];
	sprintf(SW_prefix_permanent, "%s/%s", grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS], SW_Weather.name_prefix);


  _init_soilwat_outputs(grid_files[GRID_FILE_SOILWAT2_OUTPUT]);
    SW_OUT_set_ncol(); // set number of output columns
	SW_OUT_set_colnames(); // set column names for output files
    if (_getNumberSOILWAT2OutputCells() > 0) {
		SW_OUT_create_summary_files();
		// allocate `p_OUT` and `p_OUTsd` arrays to aggregate SOILWAT2 output across iterations
		setGlobalSTEPWAT2_OutputVariables();
    }

	for (iter = 1; iter <= SuperGlobals.runModelIterations; iter++)
	{ //for each iteration
		/*
		 * 06/15/2016 (akt) Added resetting correct historical weather file path,
		 * as it was resetting to original path value (that was not correct for grid version)from input file after every iteration
		 */
		sprintf(SW_Weather.name_prefix, "%s", SW_prefix_permanent); //updates the directory of the weather files so SOILWAT2 can find them

		// Initialize the plot for each grid cell
		for (i = 0; i < grid_Rows; i++){
			for (j = 0; j < grid_Cols; j++){
				load_cell(i, j);
				Plot_Initialize();
				Globals->currIter = iter;

                if (_getNumberSOILWAT2OutputCells() > 0) {
			        print_IterationSummary = (Bool) (Globals->currIter == SuperGlobals.runModelIterations);
                }
			}
		}
		unload_cell(); // Reset the global variables

		// If we used spinup we need to reset to the state of the program right after spinup.
		if (shouldSpinup){
			loadSpinupConditions();
		}

		RandSeed(SuperGlobals.randseed, &environs_rng);
		RandSeed(SuperGlobals.randseed, &mortality_rng);
		RandSeed(SuperGlobals.randseed, &resgroups_rng);
		RandSeed(SuperGlobals.randseed, &species_rng);
		RandSeed(SuperGlobals.randseed, &grid_rng);
		RandSeed(SuperGlobals.randseed, &markov_rng);

		for (year = 1; year <= SuperGlobals.runModelYears; year++)
		{ //for each year
			if(UseProgressBar){
				logProgress(iter, year, SIMULATION);
			}

            if (UseSeedDispersal){
				disperseSeeds(year);
            }

			// Allow the colonization module to run. This function MUST be
			// Called after disperseSeeds() because it modifies the
			// seedsPresent variable of each Species.
			colonize(year);
            
			for (i = 0; i < grid_Rows; i++){
				for (j = 0; j < grid_Cols; j++)
				{ //for each cell
                
                    /* Ensure that all global variables reference the specific cell */
					load_cell(i, j);

					printf("------------------------Repetition/year = %d / %d\n", iter, year);


					Globals->currYear = year;

					/* The following functions mimic ST_main.c. */

					rgroup_Establish(); 		// Establish individuals.

					Env_Generate();				// Run SOILWAT2 to generate resources.

					rgroup_PartResources();		// Distribute resources
					rgroup_Grow(); 				// Implement plant growth

					mort_Main(&killedany); 		// Mortality that occurs during the growing season

					rgroup_IncrAges(); 			// Increment ages of all plants

					grazing_EndOfYear(); 		// Livestock grazing
					
					save_annual_species_relsize(); // Save annuals before we kill them

					mort_EndOfYear(); 			// End of year mortality.

					stat_Collect(year); 		// Update the accumulators

				    killAnnuals(); 			// Kill annuals
				    killMaxage();             // Kill plants that reach max age
					proportion_Recovery(); 		// Recover from any disturbances
					killExtraGrowth(); 		// Kill superfluous growth

				} /* end model run for this cell*/
			} /* end model run for this row */
			
			unload_cell(); // Reset the global variables
		}/* end model run for this year*/

		// collects the data for the mort output,
        // i.e. fills the accumulators in ST_stats.c.
		if (MortFlags.summary){
			for (i = 0; i < grid_Rows; i++){
				for (j = 0; j < grid_Cols; j++)
				{
					load_cell(i, j);
					stat_Collect_GMort();
					stat_Collect_SMort();
				}
			}
		}
		unload_cell(); 
		//reset soilwat to initial condition
		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
		for(i = 0; i < grid_Rows; ++i){
			for(j = 0; j < grid_Cols; ++j){
				load_cell(i, j);
                int realYears = SuperGlobals.runModelYears;
                SuperGlobals.runModelYears *= _getNumberSOILWAT2OutputCells();
				SXW_Reset(gridCells[i][j].mySXW->f_watin);
				unload_cell();
                SuperGlobals.runModelYears = realYears;
			}
		}
		Mem_Free(SW_Soilwat.hist.file_prefix);
		SW_Soilwat.hist.file_prefix = NULL;
		ChDir("..");

	} /* end iterations */

	if(UseProgressBar){
		logProgress(0, 0, OUTPUT);
	}

	// Output all of the mort and BMass files for each cell.
	for (i = 0; i < grid_Rows; i++){
		for (j = 0; j < grid_Cols; j++)
		{
			int cell = j + (i * grid_Cols);
			load_cell(i, j);
			char fileMort[1024], fileBMass[1024];

			sprintf(fileMort, "%s%d.csv", grid_files[GRID_FILE_PREFIX_MORTAVG], cell);
			sprintf(fileBMass, "%s%d.csv", grid_files[GRID_FILE_PREFIX_BMASSAVG], cell);
			parm_SetName(fileMort, F_MortAvg);
			parm_SetName(fileBMass, F_BMassAvg);

			if (MortFlags.summary && writeIndividualFiles){
				stat_Output_AllMorts();
			}
			if (BmassFlags.summary && writeIndividualFiles){
				stat_Output_AllBmass();
			}
		}
	}
	unload_cell(); // Reset the global variables

    if (recordDispersalEvents){
        outputDispersalEvents(grid_files[GRID_FILE_PREFIX_DISPERSALEVENTS]);
    }

	// Output the Bmass and Mort average statistics (if requested).
	char fileBMassCellAvg[1024], fileMortCellAvg[1024];
	if (BmassFlags.summary){
        sprintf(fileBMassCellAvg, "%s.csv", grid_files[GRID_FILE_PREFIX_BMASSCELLAVG]);
		_Output_AllCellAvgBmass(fileBMassCellAvg);
	}
    if (MortFlags.summary){
        sprintf(fileMortCellAvg, "%s.csv", grid_files[GRID_FILE_PREFIX_MORTCELLAVG]);
        _Output_AllCellAvgMort(fileMortCellAvg);
    }

  // If the user requested SOILWAT2 output this function will separate the
  // output into cell-specific files. If they didn't request SOILWAT2 output
  // this function will do nothing.
  _separateSOILWAT2Output();

	free_grid_memory();	// Free our allocated memory since we do not need it anymore
	parm_free_memory();		// Free memory allocated to the _files array in ST_params.c
	freeColonizationMemory(); // Free memory allocated to the Colonization module.
    freeDispersalMemory(); // Free memory allocated to the Seed Dispersal module.
	if(shouldSpinup) {
		freeSpinupMemory();
	}
	logProgress(0, 0, DONE);
}

/**
 * \brief Read the files.in file.
 * 
 * The files.in file specifies the locations of the other input files.
 * 
 * \sideeffect
 *     This function saves the file names it reads to \ref grid_files and 
 *     \ref grid_directories. 
 * 
 * \ingroup GRID_PRIVATE
 */
static void _init_grid_files(void)
{
	// reads the files.in file
	FILE *f;
	char buf[1024];
	int i;

	f = OpenFile(Parm_name(F_First), "r");

	for (i = 0; i < N_GRID_DIRECTORIES; i++)
	{ //0 is stepwat directory
		if (!GetALine(f, buf))
			break;
		grid_directories[i] = Str_Dup(Str_TrimLeftQ(buf));
	}
	if (i != N_GRID_DIRECTORIES)
		LogError(stderr, LOGFATAL, "Invalid files.in");

	for (i = 0; i < N_GRID_FILES; i++)
	{
		if (!GetALine(f, buf))
			break;
		grid_files[i] = Str_Dup(Str_TrimLeftQ(buf));
	}
	if (i != N_GRID_FILES)
		LogError(stderr, LOGFATAL, "Invalid files.in");

	// opens the log file...
	if (!strcmp("stdout", grid_files[GRID_FILE_LOGFILE])){
        logfp = stdout;
    }
	else {
        logfp = OpenFile(grid_files[GRID_FILE_LOGFILE], "w");
    }

	CloseFile(&f);
}

/**
 * \brief Read all gridded mode files excluding grid_setup.in. 
 * 
 * All associated fields will be populated with the values read.
 * 
 * \sideeffect
 *     This function overrides values in the \ref RGroup and \ref Species
 *     arrays so make sure you have read the non-gridded mode files before
 *     calling this function.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _init_grid_inputs(void)
{
    int i, j;

	if (UseDisturbances){
		_read_disturbances_in();
	}
	if(shouldSpinup){
		_read_spinup_species();
	}
	if (UseSoils) {
		_read_soils_in();
	}

    for(i = 0; i < grid_Rows; ++i) {
        for(j = 0; j < grid_Cols; ++j) {
            gridCells[i][j].DuringSpinup = FALSE;
        }
    }
}

/**
 * \brief Read \ref SXW input files
 * 
 * \ingroup GRID
 */
static void _init_SXW_inputs(Bool init_SW, char *f_roots)
{
    int realYears = SuperGlobals.runModelYears;
    SuperGlobals.runModelYears *= (grid_Cols * grid_Rows);
	SXW_Init(init_SW, f_roots);	// initializes soilwat
    SuperGlobals.runModelYears = realYears;
	if (init_SW)
	{
		char aString[2048];
		sprintf(aString, "%s/%s", grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS], SW_Weather.name_prefix);
		sprintf(SW_Weather.name_prefix, "%s", aString); //updates the directory correctly for the weather files so soilwat can find them
	}
}

static void _init_soilwat_outputs(char* fileName) {
  if(!writeSOILWAT2Output) {
    return;
  }

  char buffer[2048];
  int cell;
  FILE* inFile = fopen(fileName, "r");
  if(!inFile) {
    LogError(logfp, LOGFATAL, "Could not open SOILWAT2 output cells file."
             "\n\tFile named \"%s\"", fileName);
  }

  if(!GetALine(inFile, buffer)) {
    fclose(inFile);
    return;
  }

  while(sscanf(buffer, "%d,%s", &cell, buffer) == 2) {
    if(cell < 0 || cell > (grid_Cols * grid_Rows)) {
      LogError(logfp, LOGWARN, "Invalid cell (%d) specified in file \"%s\"",
               cell, fileName);
    }
    gridCells[cell / grid_Cols][cell % grid_Cols].generateSWOutput = TRUE;
  }

  if(sscanf(buffer, "%d", &cell) == 1) {
    if(cell < 0 || cell > (grid_Cols * grid_Rows)) {
      LogError(logfp, LOGWARN, "Invalid cell (%d) specified in file \"%s\"",
               cell, fileName);
    }
    gridCells[cell / grid_Cols][cell % grid_Cols].generateSWOutput = TRUE;
  }

  fclose(inFile);
}

/**
 * \brief Read in the STEPWAT2 files and populate the grid. This only needs to be called once.
 * 
 * Note that \ref gridCells must be allocated first.
 * 
 * \sideeffect
 *     Many fields in the \ref gridCells array will be assigned values.
 * 
 * \ingroup GRID_PRIVATE 
 */
static void _init_stepwat_inputs(void)
{
	int i, j; 							// Used as indices in gridCells

	ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);			// Change to folder with STEPWAT files
	parm_SetFirstName(grid_files[GRID_FILE_FILES]);	// Set the name of the STEPWAT "files.in" file

	/* Loop through all gridCells. */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			load_cell(i, j); 				     // Load this cell into the global variables
			parm_Initialize();				     // Initialize the STEPWAT variables
			gridCells[i][j].myGroup = RGroup;    // This is necessary because load_cell only points RGroup to our cell-specific
			                                     // resource group array, not to the variable that points to that array.
			gridCells[i][j].mySpecies = Species; // This is necessary because load_cell only points Species to our cell-specific
			                                     // species array, not to the variable that points to that array.
			_init_SXW_inputs(TRUE, NULL);	     // Initialize the SXW and SOILWAT variables

			// Set mySXW to the location of the newly allocated SXW
			gridCells[i][j].mySXW = getSXW();	
			// Set myTranspWindow to the location of the newly allocated transp window
			gridCells[i][j].myTranspWindow = getTranspWindow(); 
			// Set mySXWResources to the location of the newly allocated SXW resources
			gridCells[i][j].mySXWResources = getSXWResources();
		} /* End for each column */
	} /* End for each row */
	unload_cell(); // Reset the global variables

	/* Since the accumulators used in ST_stats.c are local, we need to allocate our accumulators in ST_grid.
	   The other option is to create get functions for every accumulator, which is what we do for SXW variables.
	   The reason we do not do this for accumulators is the sheer number of accumulators, which would require
	   14 get functions (as of 4/5/19). */
	_allocate_accumulators();

	ChDir("..");						// go back to the folder we started in
}

/**
 * \brief Reread input files.
 * 
 * This is useful when transitioning from [spinup](\ref SPINUP)
 * to the regular simulation. It is a quick way to reset everything
 * 
 * Be careful because this function reallocates [the grid](\ref gridCells).
 * Make sure you call \ref free_grid_memory before calling this function. 
 * 
 * \sideeffect
 *     The grid will be reallocated and all fields will be assigned the values
 *     from inputs.
 * 
 * \ingroup GRID
 */
void rereadInputs(void){
    _read_grid_setup();
    _read_files();
    _init_stepwat_inputs();
    _init_grid_inputs();
}

/**
 * \brief Allocates memory for the grid cells.
 * 
 * This function only needs to be called once.
 * 
 * \sideeffect 
 *     \ref grid_Rows * \ref grid_Cols worth of [cells](\ref CellType)
 *     will be allocated.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _allocate_gridCells(int rows, int cols){
	int i, j;
	gridCells = (CellType**) Mem_Calloc(rows, sizeof(CellType*), "_allocate_gridCells: rows");
	for(i = 0; i < rows; ++i){
		gridCells[i] = (CellType*) Mem_Calloc(cols, sizeof(CellType), "_allocate_gridCells: columns");
	}

	/* Allocate all fields specific to gridded mode. This is not necessary for fields like mySpecies 
	   since they are allocated elsewhere in the code.*/
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			// shouldBeInitialized is a dynamically allocated array
			gridCells[i][j].mySpeciesInit.shouldSpinup = (int*)
				Mem_Calloc(MAX_SPECIES, sizeof(int), "_allocate_gridCells: mySpeciesInit");

			gridCells[i][j].someKillage = (Bool*) Mem_Calloc(1, sizeof(Bool), "_allocate_gridCells: someKillage");
			
			// Allocate the cheatgrassPrecip variable for the Mortality module
			setCheatgrassPrecip(0);
			initCheatgrassPrecip();
			gridCells[i][j].myCheatgrassPrecip = getCheatgrassPrecip();
		}
	}
}

/**
 * \brief Initialize each gridCell's accumulators. 
 * 
 * Must be called after STEPWAT inputs have been read so the program knows
 * which accumulators are necessary.
 * 
 * \sideeffect 
 *     Multiple [accumulators](\ref StatType) will be allocated to each 
 *     [cell](\ref CellType). Of course, this means \ref gridCells must be
 *     allocated first.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _allocate_accumulators(void){
	int i, j;
	SppIndex sp;
	GrpIndex rg;

	/* Iterate across all cells */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			/* load_cell is not necessary for the actual accumulators, but it is necessary for
			   the ForEach loops. We still have to refer to the accumulators as
			   gridCells[i][j].<accumulator> because the ST_stats accumulators are local. */
			load_cell(i,j);

  			if (BmassFlags.dist) {
    			gridCells[i][j]._Dist = (StatType*) Mem_Calloc(1, sizeof(StatType), "_allocate_accumulators(Dist)");
		    	gridCells[i][j]._Dist->s = (struct accumulators_st *)
		               		Mem_Calloc( SuperGlobals.runModelYears,
		                           		sizeof(struct accumulators_st),
 		                         		"_allocate_accumulators(Dist)");
		  	}
 		 	if (BmassFlags.ppt) {
		    	gridCells[i][j]._Ppt = (StatType*) Mem_Calloc(1, sizeof(StatType), "_allocate_accumulators(PPT");
		    	gridCells[i][j]._Ppt->s  = (struct accumulators_st *)
		               		Mem_Calloc( SuperGlobals.runModelYears,
		                           		sizeof(struct accumulators_st),
		                          		"_allocate_accumulators(PPT)");
 		 	}
		  	if (BmassFlags.tmp) {
		    	gridCells[i][j]._Temp = (StatType*) Mem_Calloc(1, sizeof(StatType), "_allocate_accumulators(Temp)");
		    	gridCells[i][j]._Temp->s = (struct accumulators_st *)
		               		Mem_Calloc( SuperGlobals.runModelYears,
 		                          		sizeof(struct accumulators_st),
 		                         		"_allocate_accumulators(Temp)");
		  	}
		  	if (BmassFlags.grpb) {
 			   	gridCells[i][j]._Grp = (struct stat_st *)
           				Mem_Calloc( Globals->grpCount,
                		       sizeof(struct stat_st),
                		      "_allocate_accumulators(Grp)");
    			ForEachGroup(rg){
      				gridCells[i][j]._Grp[rg].s = (struct accumulators_st *)
             			Mem_Calloc( SuperGlobals.runModelYears,
                         			sizeof(struct accumulators_st),
                        			"_allocate_accumulators(Grp[rg].s)");
				}

    			if (BmassFlags.size) {
      				gridCells[i][j]._Gsize = (struct stat_st *)
             				Mem_Calloc( Globals->grpCount,
                         				sizeof(struct stat_st),
                        				"_allocate_accumulators(GSize)");
      				ForEachGroup(rg){
          				gridCells[i][j]._Gsize[rg].s = (struct accumulators_st *)
             							Mem_Calloc( SuperGlobals.runModelYears,
                         							sizeof(struct accumulators_st),
                        							"_allocate_accumulators(GSize[rg].s)");
					}
    			}

    			if (BmassFlags.pr) {
      				gridCells[i][j]._Gpr = (struct stat_st *)
             				Mem_Calloc( Globals->grpCount,
                		         		sizeof(struct stat_st),
                		        		"_allocate_accumulators(Gpr)");
      				ForEachGroup(rg){
          				gridCells[i][j]._Gpr[rg].s = (struct accumulators_st *)
             							Mem_Calloc( SuperGlobals.runModelYears,
                    		     					sizeof(struct accumulators_st),
                        							"_allocate_accumulators(Gpr[rg].s)");
					}
    			}

    			if (BmassFlags.wildfire || BmassFlags.prescribedfire) {
      				gridCells[i][j]._Gwf = (struct fire_st *)
             				Mem_Calloc( 1, sizeof(struct fire_st),
                		        		"_allocate_accumulators(Gwf)");

      				gridCells[i][j]._Gwf->wildfire = (int *) Mem_Calloc( 1,
                    		  		 sizeof(int) * SuperGlobals.runModelYears,
                    		  		 "_allocate_accumulators(Gwf->wildfire)");
      
      				gridCells[i][j]._Gwf->prescribedFire = (int **) Mem_Calloc( 1,
                      						sizeof(int **) * SuperGlobals.max_rgroups,
                       						"_allocate_accumulators(Gwf->prescribedfire");

      				ForEachGroup(rg){
        				gridCells[i][j]._Gwf->prescribedFire[rg] = (int *)
          										   Mem_Calloc( SuperGlobals.runModelYears,
                      										   sizeof(int) * SuperGlobals.runModelYears,
                      										   "_allocate_accumulators(Gwf->prescribedFire)");
      				}
    			}
    		}

  			if (MortFlags.group) {
    			gridCells[i][j]._Gestab = (struct stat_st *)
             	  		  Mem_Calloc( Globals->grpCount,
                         			sizeof(struct stat_st),
                         			"_allocate_accumulators(Gestab)");

				gridCells[i][j]._Gmort = (struct stat_st *)
           		 		 Mem_Calloc( Globals->grpCount,
                       				sizeof(struct stat_st),
                      				"_allocate_accumulators(Gmort)");
    			ForEachGroup(rg){
      				gridCells[i][j]._Gestab[rg].s = (struct accumulators_st *)
                     				Mem_Calloc( 1, sizeof(struct accumulators_st),
                                				"_allocate_accumulators(Gestab[rg].s)");
					gridCells[i][j]._Gmort[rg].s = (struct accumulators_st *)
           							Mem_Calloc( GrpMaxAge(rg),
                       							sizeof(struct accumulators_st),
                      							"_allocate_accumulators(Gmort[rg].s)");
				}
  			}

  			if (BmassFlags.sppb) {
      			gridCells[i][j]._Spp = (struct stat_st *)
               		   Mem_Calloc( Globals->sppCount,
                           		   sizeof(struct stat_st),
                          		   "_allocate_accumulators(Spp)");
      			ForEachSpecies(sp){
        			gridCells[i][j]._Spp[sp].s = (struct accumulators_st *)
               			 		  Mem_Calloc( SuperGlobals.runModelYears,
                           					sizeof(struct accumulators_st),
                          					"_allocate_accumulators(Spp[sp].s)");
				}

      			if (BmassFlags.indv) {
        			gridCells[i][j]._Indv = (struct stat_st *)
               		 		Mem_Calloc( Globals->sppCount,
                           		 		sizeof(struct stat_st),
                          		 		"_allocate_accumulators(Indv)");
        			ForEachSpecies(sp){
          				gridCells[i][j]._Indv[sp].s = (struct accumulators_st *)
               				  		Mem_Calloc( SuperGlobals.runModelYears,
                           						sizeof(struct accumulators_st),
                          						"_allocate_accumulators(Indv[sp].s)");
					}
    			}
    		}
  			if (MortFlags.species) {
    			gridCells[i][j]._Sestab = (struct stat_st *)
           					Mem_Calloc( Globals->sppCount,
                       					sizeof(struct stat_st),
                      					"_allocate_accumulators(Sestab)");

				gridCells[i][j]._Smort = (struct stat_st *)
           		 		Mem_Calloc( Globals->sppCount,
                       				sizeof(struct stat_st),
                      				"_allocate_accumulators(Smort)");

    			ForEachSpecies(sp){
      				gridCells[i][j]._Sestab[sp].s = (struct accumulators_st *)
                    				Mem_Calloc( 1, sizeof(struct accumulators_st),
                                				"_allocate_accumulators(Sestab[sp].s)");
					gridCells[i][j]._Smort[sp].s = (struct accumulators_st *)
            		        		Mem_Calloc( SppMaxAge(sp),
            		                    	   sizeof(struct accumulators_st),
             		                   	   	   "_allocate_accumulators(Smort[sp].s)");
				}
  			}

  			if (UseSeedDispersal) {
	  			gridCells[i][j]._Sreceived = Mem_Calloc( Globals->sppCount, sizeof(struct stat_st), "_allocate_accumulators(Sreceived)");

	  			ForEachSpecies(sp) {
		  			gridCells[i][j]._Sreceived[sp].s = (struct accumulators_st *)
					  									Mem_Calloc( SuperGlobals.runModelYears,
														  		   sizeof(struct accumulators_st), 
																   "_allocate_accumulators(Sreceived[sp].s)");
		  			gridCells[i][j]._Sreceived[sp].name = &Species[sp]->name[0];
	  			}
  			}

  			/* "appoint" names of columns*/
  			if (BmassFlags.grpb) {
    			ForEachGroup(rg)
      				gridCells[i][j]._Grp[rg].name = &RGroup[rg]->name[0];
  				}
  			if (MortFlags.group) {
    			ForEachGroup(rg)
      				gridCells[i][j]._Gmort[rg].name = &RGroup[rg]->name[0];
  			}
  			if (BmassFlags.sppb) {
    			ForEachSpecies(sp)
      				gridCells[i][j]._Spp[sp].name = &Species[sp]->name[0];
  			}
			if (MortFlags.species) {
				ForEachSpecies(sp)
    				gridCells[i][j]._Smort[sp].name = &Species[sp]->name[0];
			}
		} /* End for each column */
	} /* End for each row */
	unload_cell(); // Unload the cell to protect the last cell from unintended modification.
}

/**
 * \brief Free all memory allocated to the gridded mode.
 * 
 * This includes not only the [cells](\ref gridCells) but also the fields
 * inside those cells.
 * 
 * \sideeffect 
 *     \ref gridCells will be completely deallocated.
 * 
 * \ingroup GRID
 */
void free_grid_memory(void)
{
	//frees all the memory allocated in this file ST_Grid.c (most of it is dynamically allocated in _init_grid_globals() & _load_grid_globals() functions)
	int i, j;

	/* Free memory that we have allocated in ST_grid.c */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			/* Use deallocate_Globals from ST_main to deallocate global variables,
			   and free_all_sxw_memory from sxw to deallocate SXW variables. */
			load_cell(i,j);
            #ifndef STDEBUG
			deallocate_Globals(TRUE);
            #endif
			freeMortalityMemory();
			free_all_sxw_memory();
			stat_free_mem();
			unload_cell();

			Mem_Free(gridCells[i][j].mySpeciesInit.shouldSpinup);
			Mem_Free(gridCells[i][j].someKillage);

			Mem_Free(gridCells[i][j].mySoils.depth);
            Mem_Free(gridCells[i][j].mySoils.evco);
            Mem_Free(gridCells[i][j].mySoils.gravel);
            Mem_Free(gridCells[i][j].mySoils.imperm);
            Mem_Free(gridCells[i][j].mySoils.matricd);
            Mem_Free(gridCells[i][j].mySoils.pclay);
            Mem_Free(gridCells[i][j].mySoils.psand);
            Mem_Free(gridCells[i][j].mySoils.soiltemp);
            Mem_Free(gridCells[i][j].mySoils.trco_forb);
            Mem_Free(gridCells[i][j].mySoils.trco_grass);
            Mem_Free(gridCells[i][j].mySoils.trco_shrub);
            Mem_Free(gridCells[i][j].mySoils.trco_tree);
		}
	}

	for(i = 0; i < grid_Rows; ++i){
		Mem_Free(gridCells[i]);
	}
	Mem_Free(gridCells);
}

/**
 * \brief Load \ref gridCells[row][col] into the 
 *        [global variables](\ref ST_globals.h).
 * 
 * \param row the row in the \ref gridCells 2d array.
 * \param col the column in the \ref gridCells 2d array.
 * 
 * After calling this function all variables like \ref Species, \ref RGroup,
 * and \ref Env will point to a specific cell. It will also pass along the 
 * cell-specific structs like \ref SXWResources to their specific modules.
 * 
 * Any call to this function should have an accompanying call to 
 * \ref unload_cell().
 * 
 * \sideeffect
 *     All global variables will point to the specified cell.
 * 
 * \ingroup GRID
 */
void load_cell(int row, int col){
    /* RGroup for this cell */
	RGroup = gridCells[row][col].myGroup;

	/*Species for this cell */
	Species = gridCells[row][col].mySpecies;

	/* Succulents corresponding to this cell */
	Succulent = &gridCells[row][col].mySucculent;

	/* This cell's environment. We expect each cell to
	 * have slightly different weather each year */
	Env = &gridCells[row][col].myEnvironment;

	/* Cell's plot data */
	Plot = &gridCells[row][col].myPlot;

	/* Global variables corresponding to this cell */ 
	Globals = &gridCells[row][col].myGlobals;

	/* TRUE if this cell is in spinup mode */
	DuringSpinup = gridCells[row][col].DuringSpinup;

	/* This cell's cheatgrass-wildfire parameters */
	setCheatgrassPrecip(gridCells[row][col].myCheatgrassPrecip);

	_SomeKillage = gridCells[row][col].someKillage;
	UseCheatgrassWildfire = gridCells[row][col].UseCheatgrassWildfire;

	/* Copy this cell's accumulators into the local accumulators in ST_stats.c */
	stat_Copy_Accumulators(gridCells[row][col]._Dist, gridCells[row][col]._Ppt, gridCells[row][col]._Temp,
	                       gridCells[row][col]._Grp, gridCells[row][col]._Gsize, gridCells[row][col]._Gpr,
						   gridCells[row][col]._Gmort, gridCells[row][col]._Gestab, gridCells[row][col]._Spp,
						   gridCells[row][col]._Indv, gridCells[row][col]._Smort, gridCells[row][col]._Sestab,
						   gridCells[row][col]._Sreceived, gridCells[row][col]._Gwf, gridCells[row][col].stats_init);

	/* Copy this cell's SXW variables into the local variables in sxw.c */
	copy_sxw_variables(gridCells[row][col].mySXW, gridCells[row][col].mySXWResources, gridCells[row][col].myTranspWindow);

    prepare_IterationSummary = gridCells[row][col].generateSWOutput;

	// If we have read in the soil information num_layers will be > 0.
	// Otherwise we haven't read the file so there is no point wasting time on this.
	if(gridCells[row][col].mySoils.num_layers > 0){
		RealD soilRegionsLowerBounds[3] = { 30, 70, 100 };
		set_soillayers(gridCells[row][col].mySoils.num_layers, gridCells[row][col].mySoils.depth, gridCells[row][col].mySoils.matricd,
	    	           gridCells[row][col].mySoils.gravel, gridCells[row][col].mySoils.evco, gridCells[row][col].mySoils.trco_grass,
					   gridCells[row][col].mySoils.trco_shrub, gridCells[row][col].mySoils.trco_tree, gridCells[row][col].mySoils.trco_forb,
					   gridCells[row][col].mySoils.psand, gridCells[row][col].mySoils.pclay, gridCells[row][col].mySoils.imperm,
				   	   gridCells[row][col].mySoils.soiltemp, 3, soilRegionsLowerBounds);
	}
}

/**
 * \brief Nullify all global variables. Prevents accidental modification of
 *        [cells](\ref CellType).
 * 
 * Any call to \ref load_cell should be followed by a call to this function.
 * 
 * \sideeffect
 *     All global variables will be nullified, i.e. point to 0.
 * 
 * \ingroup GRID
 */
void unload_cell(){
	Species = NULL;
	RGroup = NULL;
	Succulent = NULL;
	Env = NULL;
	Plot = NULL;
	Globals = NULL;
	// Nullify the accumulators
	stat_Copy_Accumulators(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,FALSE);
	// Nullify sxw
	copy_sxw_variables(NULL,NULL,NULL);
}

/** 
 * \brief Similar to the \ref GetALine function in \ref filefuncs.c, except 
 *        this one checks for carriage return characters and doesn't deal with
 *        whitespace.
 * 
 * It treats '\\r', '\\n', and '\\r\\n' all like they are valid line feeds.
 * 
 * \ingroup GRID_PRIVATE
 */
static Bool GetALine2(FILE *f, char buf[], int limit)
{
	//only read limit characters
	int i = 0, aChar;
	aChar = getc(f);
	if (aChar == EOF)
		return FALSE;
	while ((i < (limit - 1)) && aChar != EOF && aChar != '\r' && aChar != '\n')
	{
		buf[i++] = (char) aChar;
		aChar = getc(f);
	}
	if (aChar == '\r') //this part handles the '\r\n' case
		if (getc(f) != '\n')
			fseek(f, -1, SEEK_CUR); //back up one character in the file because we didn't find a new-line character

	buf[i] = '\0';
	return TRUE;
}

/**
 * \brief Reads the grid disturbance CSV. 
 * 
 * This function will override disturbance inputs from non-gridded mode.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _read_disturbances_in(void)
{
	FILE *f;
	char buf[1024];
	int i, row, col, cell, num = 0;
    GrpIndex rg;

	f = OpenFile(grid_files[GRID_FILE_DISTURBANCES], "r");

	GetALine2(f, buf, 1024); // gets rid of the first line (since it just defines the columns)
	for (i = 0; i < grid_Cells; i++)
	{
	    row = i / grid_Cols;
	    col = i % grid_Cols;

	    load_cell(row, col);
		
		if (!GetALine2(f, buf, 1024))
			break;

        ForEachGroup(rg) {
            num = sscanf(buf, "%d,%u,%u,%u,%hu,%hu,%f,%hu,%f,%hu,%u", &cell,
                &Globals->pat.use, &Globals->mound.use, &Globals->burrow.use,
				&RGroup[rg]->killyr, &RGroup[rg]->killfreq_startyr, 
				&RGroup[rg]->killfreq, &RGroup[rg]->extirp, 
				&RGroup[rg]->grazingfrq, &RGroup[rg]->grazingfreq_startyr, 
				&gridCells[row][col].UseCheatgrassWildfire);
		}

		if (num != 11)
			LogError(logfp, LOGFATAL, "Invalid %s file line %d wrong",
					grid_files[GRID_FILE_DISTURBANCES], i + 2);
	}
	if (i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file wrong number of cells",
				grid_files[GRID_FILE_DISTURBANCES]);

    unload_cell();
	CloseFile(&f);
}

/**
 * \brief Iterates through s until it find nSeperators worth of seperators
 * 
 * Used to do most of the parsing in the \ref _read_soils_in() function
 * 
 * \param s is a char* array.
 * \param seperator is the character used as a separator, for example tab or space.
 * \param nSeperators is the number of separators to read.
 * 
 * \return index of the character following the last separator. 
 * 
 * \ingroup GRID_PRIVATE
 */
static int _get_value_index(char* s, char seperator, int nSeperators)
{
	int i = 0, sep = 0;
	while (*s)
	{
		i++;
		if (*s++ == seperator) //checks if the current char equals the seperator... and then increments the char pointer
			if (++sep == nSeperators) //needs ++sep, not sep++
				break;
	}
	return i;
}

/* Read the grid_soils.csv file and assign values to all gridCells.mySoils variables. */
static void _read_soils_in(void){
	int i, j, row, col, lineReadReturnValue;
	char buf[4096];

	/* tempSoil is allocated the maximum amount of memory that a SoilType could need.
	   It will serve to read in parameters. */
	SoilType tempSoil;
	tempSoil.depth = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.evco = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.gravel = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.imperm = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.matricd = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.pclay = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.psand = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.soiltemp = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.trco_forb = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.trco_grass = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.trco_shrub = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");
	tempSoil.trco_tree = Mem_Calloc(MAX_LAYERS, sizeof(RealF), "_read_soils_in: tempSoil");

	FILE* f = OpenFile(grid_files[GRID_FILE_SOILS], "r");
	if(!GetALine(f, buf)){ // Throw out the header line.
		LogError(logfp, LOGFATAL, "%s file empty.", grid_files[GRID_FILE_SOILS]);
	}

	for(i = 0; i < grid_Cells; ++i){
		row = i / grid_Cols;
	    col = i % grid_Cols;
		load_cell(row, col);

		if(!GetALine(f, buf)){
			LogError(logfp, LOGFATAL, "Too few lines in %s", grid_files[GRID_FILE_SOILS]);
		}
		lineReadReturnValue = _read_soil_line(buf, &tempSoil, 0);
		if (lineReadReturnValue == SOIL_READ_FAILURE){
			LogError(logfp, LOGFATAL, "Error reading %s file.", grid_files[GRID_FILE_SOILS]);
		} 
		/* If _read_soil_line didnt return SUCCESS or FAILURE, 
		   it returned a cell number to copy. */
		else if (lineReadReturnValue != SOIL_READ_SUCCESS){
			if(lineReadReturnValue > i){
				LogError(logfp, LOGFATAL, "%s: Attempted to copy values that have not been read yet.\n"
				                          "\tIf you want to copy a soil make sure you define the layers"
										  "the FIRST time you use it.", grid_files[GRID_FILE_SOILS]);
			}
			_copy_soils(&gridCells[lineReadReturnValue / grid_Cols][lineReadReturnValue % grid_Cols].mySoils, &gridCells[row][col].mySoils);
		}
		/* If we get here we have successfully populated the first layer of soil. 
		   Now we must populate the rest. */
		else {
			for(j = 1; j < tempSoil.num_layers; ++j){
				if(!GetALine(f, buf)){
					LogError(logfp, LOGFATAL, "Too few lines in %s", grid_files[GRID_FILE_SOILS]);
				}
				lineReadReturnValue = _read_soil_line(buf, &tempSoil, j);
				if(lineReadReturnValue != SOIL_READ_SUCCESS){
					LogError(logfp, LOGFATAL, "Different behavior is specified between layers %d and %d of"
					                          " cell %d in file %s. (Perhaps you specified a cell to copy in one"
											  " but not the other?)", j, j+1, i, grid_files[GRID_FILE_SOILS]);
				}
			}
			/* And finally copy the temporary soil into the grid. */
			_copy_soils(&tempSoil, &gridCells[row][col].mySoils);
		}
	}

	unload_cell();
	CloseFile(&f);

	Mem_Free(tempSoil.soiltemp);
	Mem_Free(tempSoil.trco_forb);
	Mem_Free(tempSoil.trco_grass);
	Mem_Free(tempSoil.trco_shrub);
	Mem_Free(tempSoil.trco_tree);
	Mem_Free(tempSoil.depth);
	Mem_Free(tempSoil.evco);
	Mem_Free(tempSoil.gravel);
	Mem_Free(tempSoil.imperm);
	Mem_Free(tempSoil.matricd);
	Mem_Free(tempSoil.pclay);
	Mem_Free(tempSoil.psand);
}

/**
 * \brief Reads a line of soil input from buf into destination.
 * 
 * \param buf line of soil input to parse.
 * \param destination \ref SoilType struct to fill.
 * \param layer layer number of destination to fill (0 indexed).
 * 
 * \return SOIL_READ_FAILURE if buf is incorrectly formatted.
 * \return SOIL_READ_SUCCESS if destination is populated correctly.
 * \return Otherwise returns the cell number that destination should copy.
 * 
 * \ingroup GRID_PRIVATE
 */ 
static int _read_soil_line(char* buf, SoilType* destination, int layer){
	int entriesRead, cellToCopy, cellNum, layerRead;
	entriesRead = sscanf(buf, "%d,,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s",
			                    &cellNum, &destination->num_layers, &layerRead, &destination->depth[layer], 
								&destination->matricd[layer], &destination->gravel[layer], 
								&destination->evco[layer], &destination->trco_grass[layer], 
								&destination->trco_shrub[layer], &destination->trco_tree[layer],
								&destination->trco_forb[layer], &destination->psand[layer], 
								&destination->pclay[layer], &destination->imperm[layer], 
								&destination->soiltemp[layer], destination->rootsFile);

	if(cellNum > grid_Cells){
		LogError(
			logfp, LOGFATAL,
			"%s: cell number (id=%d) is larger than number of cells in grid (n=%d).",
			grid_files[GRID_FILE_SOILS], cellNum, grid_Cells
		);
	}


	/* If the user specified a cell to copy we will perform the copy, regardless
	   of whether or not they also entered parameters. */
	if(entriesRead == 1){
		entriesRead = sscanf(buf, "%d,%d", &cellNum, &cellToCopy);

		if(cellToCopy > grid_Cells){
			LogError(
				logfp, LOGFATAL,
				"%s: cell to copy (id=%d) not present in grid (n=%d).",
				grid_files[GRID_FILE_SOILS], cellToCopy, grid_Cells
			);
		}

		if(entriesRead == 2){
			return cellToCopy;
		} else {
			return SOIL_READ_FAILURE;
		}
	}

	if(entriesRead == 16) {
		if(layerRead > destination->num_layers){
			LogError(
				logfp, LOGFATAL,
				"%s: cell %d has too many soil layers (%d for max=%d).",
				grid_files[GRID_FILE_SOILS], cellNum, layerRead, destination->num_layers
			);
		}

		return SOIL_READ_SUCCESS;
	} else {
		return SOIL_READ_FAILURE;
	}
}

/**
 * \brief Copy one SoilType variable to another.
 * 
 * \param src: An allocated and assigned SoilType
 * \param dest: An unallocated SoilType
 * 
 * \sideeffect
 *     dest will be allocated memory, so do not call this function if dest is
 *     already allocated. 
 * 
 * \ingroup GRID_PRIVATE
 */
void _copy_soils(SoilType* src, SoilType* dest){
	int i;

	dest->num_layers = src->num_layers;
	dest->gravel = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: gravel");
	dest->depth = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: depth");
	dest->matricd = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: matricd");
	dest->evco = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: evco");
	dest->trco_grass = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: trco_grass");
	dest->trco_shrub = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: trco_shrub");
	dest->trco_tree = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: trco_tree");
	dest->trco_forb = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: trco_forb");
	dest->psand = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: psand");
	dest->pclay = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: pclay");
	dest->imperm = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: imperm");
	dest->soiltemp = Mem_Calloc(src->num_layers, sizeof(RealF), "_copy_Soils: soiltemp");

	for(i = 0; i < src->num_layers; ++i){
		dest->gravel[i] = src->gravel[i];
		dest->evco[i] = src->evco[i];
		dest->depth[i] = src->depth[i];
		dest->matricd[i] = src->matricd[i];
		dest->trco_shrub[i] = src->trco_shrub[i];
		dest->trco_grass[i] = src->trco_grass[i];
		dest->trco_forb[i] = src->trco_forb[i];
		dest->trco_tree[i] = src->trco_tree[i];
		dest->psand[i] = src->psand[i];
		dest->pclay[i] = src->pclay[i];
		dest->imperm[i] = src->imperm[i];
		dest->soiltemp[i] = src->soiltemp[i];
	}
}

/**
 * \brief Read the species spinup CSV. 
 * 
 * This function only needs to be called if the user requests spinup.
 * 
 * \sideeffect
 *     The species spinup information of each [cell](\ref gridCells)
 *     will be populated.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _read_spinup_species(void)
{
	FILE *f;
	int i, j, num, cell, do_copy, copy_cell, useSpinup, seeds_Avail,
	    row, col, copy_cell_row, copy_cell_col;
	Bool isAnyCellOnForSpinup = FALSE;
	char buf[4096];

	//open the file/do the reading
	f = OpenFile(grid_files[GRID_FILE_SPINUP_SPECIES], "r");

	GetALine2(f, buf, 4096); // gets rid of the first line (since it just defines the columns)... it's only there for user readability
	for (i = 0; i < grid_Cells; i++)
	{
		useSpinup = FALSE;
	    row = i / grid_Cols;
	    col = i % grid_Cols;

	    load_cell(row, col);

		if (!GetALine2(f, buf, 4096))
			break;

		num = sscanf(buf, "%d,%d,%d", &cell, &do_copy, &copy_cell);

		copy_cell_row = copy_cell / grid_Cols;
		copy_cell_col = copy_cell % grid_Cols;

		if (num != 3)
			LogError(logfp, LOGFATAL, "Invalid %s file", grid_files[GRID_FILE_SPINUP_SPECIES]);

		int stringIndex = _get_value_index(buf, ',', 3); //gets us the index of the string that is right after what we just parsed in

		if (do_copy == 1 && copy_cell > -1 && copy_cell < grid_Cells
				&& cell != 0 && copy_cell < cell)
		{ //copy this cells values from a previous cell's
			for (j = 0; j < Globals->sppCount; j++)
				gridCells[row][col].mySpeciesInit.shouldSpinup[j] =
				    gridCells[copy_cell_row][copy_cell_col].mySpeciesInit.shouldSpinup[j];
			gridCells[row][col].mySpeciesInit.useSpinup =
			    gridCells[copy_cell_row][copy_cell_col].mySpeciesInit.useSpinup;
			continue;
		}
		else if (do_copy == 1)
			LogError(logfp, LOGFATAL,
					"Invalid %s file line %d invalid copy_cell attempt",
					grid_files[GRID_FILE_SPINUP_SPECIES], i + 2);

		//going through each species
		SppIndex s;
		ForEachSpecies(s)
		{
			num = sscanf(&buf[stringIndex], "%d,", &seeds_Avail);
			if (num != 1){
				LogError(logfp, LOGFATAL, "Invalid %s file line %d invalid species input",
						 grid_files[GRID_FILE_SPINUP_SPECIES], i + 2);
			}

			if(seeds_Avail){
				useSpinup = TRUE;
				isAnyCellOnForSpinup = TRUE;
			}

			gridCells[row][col].mySpeciesInit.shouldSpinup[s] = seeds_Avail;
			stringIndex += _get_value_index(&buf[stringIndex], ',', 1);
		}

		gridCells[row][col].mySpeciesInit.useSpinup = useSpinup;
	}

	if (i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file, not enough cells",
				grid_files[GRID_FILE_SPINUP_SPECIES]);

	if(!isAnyCellOnForSpinup){
		LogError(logfp, LOGWARN, "Spinup is on, but no species are turned on for spinup inside %s.",
				grid_files[GRID_FILE_SPINUP_SPECIES]);
	}

    unload_cell();
	CloseFile(&f);
}

/**
 * \brief Read the maxrgroupspecies file.
 * 
 * This will populate the max [RGroup](\ref RGROUP) and [Species](\ref SPECIES)
 * information in the \ref SuperGlobals struct.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _read_maxrgroupspecies(void)
{
    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    parm_SetName(grid_files[GRID_FILE_MAXRGROUPSPECIES], F_MaxRGroupSpecies);
    maxrgroupspecies_init();
    ChDir("..");
}

/**
 * \brief Read the non-gridded mode files.in file.
 * 
 * \ingroup GRID_PRIVATE
 */
static void _read_files(void)
{
    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    files_init();
    ChDir("..");
}

/**
 * \brief Reads the grid setup file and allocates \ref gridCells. 
 * 
 * \ingroup GRID_PRIVATE
 */
static void _read_grid_setup(void)
{
    FILE *f;
    char buf[1024];
    int i, j;

    f = OpenFile(grid_files[GRID_FILE_SETUP], "r");

    GetALine(f, buf);
    i = sscanf(buf, "%d %d", &grid_Rows, &grid_Cols);
    if (i != 2)
        LogError(logfp, LOGFATAL,
                 "Invalid grid setup file (rows/cols line wrong)");

    grid_Cells = grid_Cols * grid_Rows;

    if (grid_Cells > MAX_CELLS)
        LogError(logfp, LOGFATAL,
                 "Number of cells in grid exceeds MAX_CELLS defined in ST_defines.h");

    /* Allocate the 2d array of cells now that we know how many we need */
    _allocate_gridCells(grid_Rows, grid_Cols);

    GetALine(f, buf);
    i = sscanf(buf, "%u", &UseDisturbances);
    if (i != 1)
        LogError(logfp, LOGFATAL,
                 "Invalid grid setup file (disturbances line wrong)");

    GetALine(f, buf);
    i = sscanf(buf, "%u", &UseSoils);
    if (i != 1)
        LogError(logfp, LOGFATAL, "Invalid grid setup file (soils line wrong)");

    GetALine(f, buf);
    i = sscanf(buf, "%d", &j);
    if (i != 1)
        LogError(logfp, LOGFATAL,
                 "Invalid grid setup file (seed dispersal line wrong)");
    UseSeedDispersal = itob(j);

	GetALine(f, buf);
	i = sscanf(buf, "%d", &shouldSpinup);
	if(i < 1){
		LogError(logfp, LOGFATAL, "Invalid grid setup file (Spinup line wrong)");
	}

	GetALine(f, buf);
	i = sscanf(buf, "%hd", &SuperGlobals.runSpinupYears);
	if(i < 1){
		LogError(logfp, LOGFATAL, "Invalid grid setup file (Spinup years line wrong)");
	}

	GetALine(f, buf);
	i = sscanf(buf, "%u", &writeIndividualFiles);
	if(i < 1){
		LogError(logfp, LOGFATAL, "Invalid grid setup file (Individual output line wrong)");
	}

	GetALine(f, buf);
	if (sscanf(buf, "%u", &recordDispersalEvents) != 1) {
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal events output line\n", grid_files[GRID_FILE_SETUP]);
    }

    CloseFile(&f);
}

/**
 * \brief Output a master .csv file containing the averages across all cells. 
 * 
 * \ingroup GRID_PRIVATE
 */
void _Output_AllCellAvgBmass(const char * filename){
	int i, j, year, nobs = 0;	//for iterating
	GrpIndex rg;	//for iterating
	SppIndex sp;	//for iterating

	/* One accumulator for every accumulator in ST_stats.c */
	float ppt, pptstd, pptsos, temp, tempstd, tempsos, dist, wildfire, grp[SuperGlobals.max_rgroups], grpstd[SuperGlobals.max_rgroups], 
		  grpsos[SuperGlobals.max_rgroups], gsize[SuperGlobals.max_rgroups], gpr[SuperGlobals.max_rgroups], 
		  gprsos[SuperGlobals.max_rgroups], gprstd[SuperGlobals.max_rgroups], prescribedfire[SuperGlobals.max_rgroups],
		  spp[SuperGlobals.max_spp_per_grp * SuperGlobals.max_rgroups],
		  indv[SuperGlobals.max_spp_per_grp * SuperGlobals.max_rgroups];

	char buf[2048], tbuf[2048];	// Two buffers: one for accumulating and one for formatting.
	char sep = BmassFlags.sep;	// Separator specified in inputs
	size_t len_buf;

	FILE* file;
	file = fopen(filename, "w");

	buf[0]='\0';

	load_cell(0, 0);

	if(BmassFlags.header){
		make_header_with_std(buf);
		fprintf(file, "%s", buf);
	}

	for(year = 0; year < SuperGlobals.runModelYears; ++year){
		*buf = '\0';	//Initialize buffer as empty string

		/* -------- Initialize all accumulators to 0 ------- */
		nobs = 0;
		ppt = 0;
		pptstd = 0;
		pptsos = 0;		
		temp = 0;
		tempstd = 0;
		tempsos = 0;
		dist = 0;
		wildfire = 0;
		ForEachGroup(rg){
			grp[rg] = 0;
			grpsos[rg] = 0;
			grpstd[rg] = 0;
			gsize[rg] = 0;
			gpr[rg] = 0;
			gprsos[rg] = 0;
			gprstd[rg] = 0; 
			prescribedfire[rg] = 0;
		}
		ForEachSpecies(sp){
			spp[sp] = 0;
			indv[sp] = 0;
		}
		/* ------ End Initialize all accumulators to 0 ------ */

		for(i = 0; i < grid_Rows; ++i){ // For each row
			for(j = 0; j < grid_Cols; ++j){ // For each column
				/* ------------- Accumulate requested output ----------------- */
				nobs++;
				if(BmassFlags.ppt) {
					
					float old_ppt_ave = ppt;
					ppt = get_running_mean(nobs, ppt, gridCells[i][j]._Ppt->s[year].ave);
					pptsos += get_running_sqr(old_ppt_ave, ppt, gridCells[i][j]._Ppt->s[year].ave);
					pptstd = final_running_sd(nobs, pptsos);
				}
				if(BmassFlags.tmp) {
					float old_temp_ave = temp;
					temp = get_running_mean(nobs, temp, gridCells[i][j]._Temp->s[year].ave);
					tempsos += get_running_sqr(old_temp_ave, temp, gridCells[i][j]._Temp->s[year].ave);
					tempstd = final_running_sd(nobs, tempsos);
				}
				if(BmassFlags.dist) dist += gridCells[i][j]._Dist->s[year].nobs;
				if(BmassFlags.grpb) {
					if(BmassFlags.wildfire){
						wildfire += gridCells[i][j]._Gwf->wildfire[year];
					}
					ForEachGroup(rg){
						float old_grp_ave = grp[rg];
						grp[rg] = get_running_mean(nobs, grp[rg], gridCells[i][j]._Grp[rg].s[year].ave);
						grpsos[rg] += get_running_sqr(old_grp_ave, grp[rg], gridCells[i][j]._Grp[rg].s[year].ave);
						grpstd[rg] = final_running_sd(nobs, grpsos[rg]);
						if(BmassFlags.size){
							gsize[rg] += gridCells[i][j]._Gsize[rg].s[year].ave;
						}
						if(BmassFlags.pr){
							float old_gpr_ave = gpr[rg];
							gpr[rg] = get_running_mean(nobs, gpr[rg], gridCells[i][j]._Gpr[rg].s[year].ave);
							gprsos[rg] += get_running_sqr(old_gpr_ave, gpr[rg], gridCells[i][j]._Gpr[rg].s[year].ave);
							gprstd[rg] = final_running_sd(nobs, gprsos[rg]);
						}
						if(BmassFlags.prescribedfire){
							prescribedfire[rg] += gridCells[i][j]._Gwf->prescribedFire[rg][year];
						}
					} // End ForEachGroup
				} // End grpb
				if(BmassFlags.sppb){
					ForEachSpecies(sp){
						spp[sp] += gridCells[i][j]._Spp[sp].s[year].ave;
						if(BmassFlags.indv){
							indv[sp] += gridCells[i][j]._Indv[sp].s[year].ave;
						}
					} // End ForEachSpecies
				} // End sppb
				/* ------------ End Accumulate requested output --------------- */	
			} // End for each column
		} // End for each row

		/* ------------------ Average all accumulators ----------------- */
		dist /= grid_Cells;
		wildfire /= grid_Cells;
		ForEachGroup(rg){
			gsize[rg] /= grid_Cells;
			prescribedfire[rg] /= grid_Cells;
		}
		ForEachSpecies(sp){
			spp[sp] /= grid_Cells;
			indv[sp] /= grid_Cells;
		}
		/* --------------- End average all accumulators --------------- */

		/* ----------------- Generate output string ------------------- */
		/* buf will hold the entire string. tbuf will format the output */
		if(BmassFlags.yr){
			sprintf(buf,"%d%c", year+1, sep);
		}
		if(BmassFlags.dist){
			sprintf(tbuf, "%f%c", dist, sep);
			strcat(buf, tbuf);
		}
		if(BmassFlags.ppt){
			sprintf(tbuf, "%f%c%f%c", ppt, sep, pptstd, sep);
			strcat(buf, tbuf);
		}
		if (BmassFlags.pclass) {
			sprintf(tbuf, "\"NA\"%c", sep);
      		strcat(buf, tbuf);
    	}
		if(BmassFlags.tmp){
			sprintf(tbuf, "%f%c%f%c", temp, sep, tempstd, sep);
			strcat(buf, tbuf);
		}
		if(BmassFlags.grpb){
			if(BmassFlags.wildfire){
				sprintf(tbuf, "%f%c", wildfire, sep);
				strcat(buf, tbuf);
			}
			ForEachGroup(rg){
				sprintf(tbuf, "%f%c%f%c", grp[rg], sep, grpstd[rg], sep);
				strcat(buf, tbuf);

				if(BmassFlags.size){
					sprintf(tbuf, "%f%c", gsize[rg], sep);
					strcat(buf, tbuf);
				}
				if(BmassFlags.pr){
					sprintf(tbuf, "%f%c%f%c", gpr[rg], sep, gprstd[rg], sep);
					strcat(buf, tbuf);
				}
				if(BmassFlags.prescribedfire){
					sprintf(tbuf, "%f%c", prescribedfire[rg], sep);
					strcat(buf, tbuf);
				}
			}
		}
		if(BmassFlags.sppb){
			ForEachSpecies(sp){
				sprintf(tbuf, "%f%c", spp[sp], sep);
				strcat(buf, tbuf);
				if(BmassFlags.indv){
					sprintf(tbuf, "%f%c", indv[sp], sep);
					strcat(buf, tbuf);
				}
			}
		}
		/* --------------- End generate output string ---------------- */

		// remove the last (and superfluous) `sep` and replace it with a '\0'
		len_buf = strlen(buf);
		if (len_buf > 1) {
			buf[len_buf - 1] = 0;
		}

		fprintf(file, "%s\n", buf); // Finally, print this line
	} // End for each year

	unload_cell();
	fclose(file); // Close the file
}

/**
 * \brief Output the average mortality across cells. 
 * 
 * \ingroup GRID_PRIVATE
 */
void _Output_AllCellAvgMort(const char* fileName){
    if (!MortFlags.summary) return;

    /* We need a cell loaded for the ForEach loops. */
    load_cell(0,0);

    FILE *file;
    IntS age;
    GrpIndex rg;
    SppIndex sp;
    char sep = MortFlags.sep;

    int row, col, nobs = 0;

    float Gestab[SuperGlobals.max_rgroups],
          Sestab[SuperGlobals.max_spp_per_grp * SuperGlobals.max_rgroups],
          Gmort[SuperGlobals.max_rgroups][Globals->Max_Age],
          Smort[SuperGlobals.max_spp_per_grp * SuperGlobals.max_rgroups][Globals->Max_Age];

    file = OpenFile( fileName, "w");

    /* --------------------- Initialize values ----------------------- */
    ForEachSpecies(sp){
        Sestab[sp] = 0;
        for(age = 0; age < Globals->Max_Age; ++age){
            Smort[sp][age] = 0;
        }
    }
    ForEachGroup(rg){
        Gestab[rg] = 0;
        for(age = 0; age < Globals->Max_Age; ++age){
            Gmort[rg][age] = 0;
        }
    }
    /* ------------------ End initializing values -------------------- */

    /* ----------------------------- Calculate averaged values ------------------------------------- */
    for(row = 0; row < grid_Rows; row++){
        for(col = 0; col < grid_Cols; col++){
            nobs++;
            if (MortFlags.group) {
                ForEachGroup(rg){
                    Gestab[rg] = get_running_mean(nobs, Gestab[rg], gridCells[row][col]._Gestab[rg].s[0].ave);
                }
            }

            if (MortFlags.species) {
                ForEachSpecies(sp){
                    Sestab[sp] = get_running_mean(nobs, Sestab[sp], gridCells[row][col]._Sestab[sp].s[0].ave);
                }
            }

            /* print one line of kill frequencies per age */
            for(age=0; age < Globals->Max_Age; age++) {
                if (MortFlags.group) {
                    ForEachGroup(rg){
                        Gmort[rg][age] = get_running_mean(nobs, Gmort[rg][age], gridCells[row][col]._Gmort[rg].s[age].ave);
                    }
                }
                if (MortFlags.species) {
                    ForEachSpecies(sp) {
                        Smort[sp][age] = get_running_mean(nobs, Smort[sp][age], gridCells[row][col]._Smort[sp].s[age].ave);
                    }
                }
            }
        }
    }
    /* --------------------------- End calculating averaged values -------------------------------- */

    /* --------------------- Print header ------------------ */
    fprintf(file, "Age");
    if (MortFlags.group) {
        ForEachGroup(rg) {
            fprintf(file,"%c%s", sep, RGroup[rg]->name);
        }
    }
    if (MortFlags.species) {
        ForEachSpecies(sp)
            fprintf(file,"%c%s", sep, Species[sp]->name);
    }
    fprintf(file,"\n");
    fprintf(file,"Estabs");
    /* ---------------------- End header ------------------- */

    /* --------------------- Print values ------------------ */
    if(MortFlags.group){
        ForEachGroup(rg){
            fprintf(file,"%c%5.1f", sep, Gestab[rg]);
        }
    }

    if (MortFlags.species) {
        ForEachSpecies(sp){
            fprintf(file,"%c%5.1f", sep, Sestab[sp]);
        }
    }

    fprintf(file,"\n");

    /* print one line of kill frequencies per age */
    for(age=0; age < Globals->Max_Age; age++) {
        fprintf(file,"%d", age+1);
        if (MortFlags.group) {
            ForEachGroup(rg)
                fprintf(file,"%c%5.1f", sep, ( age < GrpMaxAge(rg) )
                                  ? Gmort[rg][age]
                                  : 0.);
        }
        if (MortFlags.species) {
            ForEachSpecies(sp) {
                fprintf(file,"%c%5.1f", sep, ( age < SppMaxAge(sp))
                                ? Smort[sp][age]
                                : 0.);
            }
        }
        fprintf(file,"\n");
    }
    /* ----------------- End printing values --------------- */

    unload_cell();
    CloseFile(&file);
}

/**
 * \brief Separate any SOILWAT2 output into cell-specific files
 * 
 * SOILWAT2 prints output into 3 types of files: yearly, monthly and daily.
 * Each file contains the output for EVERY cell. This function iterates over
 * these files and separates them into cell-specific files.
 * 
 * \author Chandler Haukap
 * \date 18 February 2020
 * \ingroup GRID_PRIVATE
 */
static void _separateSOILWAT2Output(void){
  // If none of the cells requested SOILWAT2 output we can quit now.
  if(!_getNumberSOILWAT2OutputCells()) {
    return;
  }

  int* cells = _getSOILWAT2OutputCells();
  DIR *d;
  struct dirent *dir;
  ChDir("./Output/sw_output/");
  d = opendir(".");
  if (d)
  {
    while ((dir = readdir(d)) != NULL)
    {
      if(strncmp(dir->d_name, "sw2_daily_", 10) == 0) {
        _separateSOILWAT2DailyOutput(dir->d_name, cells);
      } else if(strncmp(dir->d_name, "sw2_monthly_", 12) == 0) {
        _separateSOILWAT2MonthlyOutput(dir->d_name, cells);
      } else if(strncmp(dir->d_name, "sw2_yearly_", 11) == 0) {
        _separateSOILWAT2YearlyOutput(dir->d_name, cells);
      }
    }
    closedir(d);
  }
  Mem_Free(cells);
  ChDir("../..");
}

/**
 * \brief Separate the daily SOILWAT2 output files into cell-specific files.
 * 
 * This function is intended to be called by \ref _separateSOILWAT2Output. It
 * is a function that separates the daily weather files generated by SOILWAT2
 * into cell-specific files.
 * 
 * \param fileName is the name of the file containing the daily data.
 * \param cellNumbers is an array of size \ref _getNumberSOILWAT2OutputCells
 *        containing the cell numbers that requested SOILWAT2 output.
 * 
 * \sa _separateSOILWAT2Output
 * \author Chandler Haukap
 * \date 18 February 2020
 * \ingroup GRID_PRIVATE
 */
static void _separateSOILWAT2DailyOutput(char* fileName, int* cellNumbers) {
  int thisDay, lastDay = 0, junk, outFileIndex;
  size_t bufsize = 9000;
  char junkBuffer[4096];
  int numCells = _getNumberSOILWAT2OutputCells();
  FILE** outFiles = Mem_Calloc(numCells, sizeof(FILE*), 
                               "_separateSOILWAT2DailyOutput");
  char* buffer = Mem_Calloc(bufsize, sizeof(char), 
                            "_separateSOILWAT2DailyOutput");
  FILE* inFile = fopen(fileName, "r");

  if(!inFile) {
    LogError(logfp, LOGFATAL, "Issue while separating SOILWAT2 output.\n\tFile"
             "\"%s\" not found.", fileName);
  }
  
  // Open the output files then copy the header over.
  getline(&buffer, &bufsize, inFile);
  for(outFileIndex = 0; outFileIndex < numCells; ++outFileIndex) {
    sprintf(junkBuffer, "Cell_%d_%s", cellNumbers[outFileIndex], fileName);
    outFiles[outFileIndex] = fopen(junkBuffer, "w");
    fprintf(outFiles[outFileIndex], "%s", buffer);
  }

  // Separate the file
  outFileIndex = 0;
  while(getline(&buffer, &bufsize, inFile) > 0) {
    sscanf(buffer, "%d,%d,%s", &junk, &thisDay, junkBuffer);
    if(lastDay > thisDay) {
      outFileIndex = (outFileIndex + 1) % numCells;
      lastDay = 0;
    } else {
      lastDay = thisDay;
    }
    fprintf(outFiles[outFileIndex], "%s", buffer);
  }

  // Close all of the output files
  for(outFileIndex = 0; outFileIndex < numCells; ++outFileIndex) {
    fclose(outFiles[outFileIndex]);
  }
  fclose(inFile);
  remove(fileName);

  Mem_Free(outFiles);
  Mem_Free(buffer);
}

/**
 * \brief Separate the monthly SOILWAT2 output files into cell-specific files.
 * 
 * This function is intended to be called by \ref _separateSOILWAT2Output. It
 * is a function that separates the monthly weather files generated by SOILWAT2
 * into cell-specific files.
 * 
 * \param fileName is the name of the file containing the monthly data.
 * \param cellNumbers is an array of size \ref _getNumberSOILWAT2OutputCells
 *        containing the cell numbers that requested SOILWAT2 output.
 * 
 * \sa _separateSOILWAT2Output
 * \author Chandler Haukap
 * \date 18 February 2020
 * \ingroup GRID_PRIVATE
 */
static void _separateSOILWAT2MonthlyOutput(char* fileName, int* cellNumbers) {
  int outFileIndex, month = 1;
  size_t bufsize = 9000;
  char junkBuffer[4096];
  int numCells = _getNumberSOILWAT2OutputCells();
  FILE** outFiles = Mem_Calloc(numCells, sizeof(FILE*), 
                               "_separateSOILWAT2MonthlyOutput");
  char* buffer = Mem_Calloc(bufsize, sizeof(char), 
                            "_separateSOILWAT2MonthlyOutput");
  FILE* inFile = fopen(fileName, "r");

  if(!inFile) {
    LogError(logfp, LOGFATAL, "Issue while separating SOILWAT2 output.\n\tFile"
             "\"%s\" not found.", fileName);
  }
  
  // Open the output files then copy the header over.
  getline(&buffer, &bufsize, inFile);
  for(outFileIndex = 0; outFileIndex < numCells; ++outFileIndex) {
    sprintf(junkBuffer, "Cell_%d_%s", cellNumbers[outFileIndex], fileName);
    outFiles[outFileIndex] = fopen(junkBuffer, "w");
    fprintf(outFiles[outFileIndex], "%s", buffer);
  }

  // Separate the file
  outFileIndex = 0;
  while(getline(&buffer, &bufsize, inFile) > 0) {
    fprintf(outFiles[outFileIndex], "%s", buffer);
    month++;
    if(month > 12){
      outFileIndex = (outFileIndex + 1) % numCells;
      month = 1;
    }
  }

  // Close all of the output files
  for(outFileIndex = 0; outFileIndex < numCells; ++outFileIndex) {
    fclose(outFiles[outFileIndex]);
  }
  fclose(inFile);
  remove(fileName);

  Mem_Free(outFiles);
  Mem_Free(buffer);
}

/**
 * \brief Separate the yearly SOILWAT2 output files into cell-specific files.
 * 
 * This function is intended to be called by \ref _separateSOILWAT2Output. It
 * is a function that separates the yearly weather files generated by SOILWAT2
 * into cell-specific files.
 * 
 * \param fileName is the name of the file containing the yearly data.
 * \param cellNumbers is an array of size \ref _getNumberSOILWAT2OutputCells
 *        containing the cell numbers that requested SOILWAT2 output.
 * 
 * \sa _separateSOILWAT2Output
 * \author Chandler Haukap
 * \date 18 February 2020
 * \ingroup GRID_PRIVATE
 */
static void _separateSOILWAT2YearlyOutput(char* fileName, int* cellNumbers) {
  int outFileIndex;
  size_t bufsize = 9000;
  char junkBuffer[4096];
  int numCells = _getNumberSOILWAT2OutputCells();
  FILE** outFiles = Mem_Calloc(numCells, sizeof(FILE*), 
                               "_separateSOILWAT2YearlyOutput");
  char* buffer = Mem_Calloc(bufsize, sizeof(char), 
                            "_separateSOILWAT2YearlyOutput");
  FILE* inFile = fopen(fileName, "r");

  if(!inFile) {
    LogError(logfp, LOGFATAL, "Issue while separating SOILWAT2 output.\n\tFile"
             "\"%s\" not found.", fileName);
  }
  
  // Open the output files then copy the header over.
  getline(&buffer, &bufsize, inFile);
  for(outFileIndex = 0; outFileIndex < numCells; ++outFileIndex) {
    sprintf(junkBuffer, "Cell_%d_%s", cellNumbers[outFileIndex], fileName);
    outFiles[outFileIndex] = fopen(junkBuffer, "w");
    fprintf(outFiles[outFileIndex], "%s", buffer);
  }

  // Separate the file
  outFileIndex = 0;
  while(getline(&buffer, &bufsize, inFile) > 0) {
    fprintf(outFiles[outFileIndex], "%s", buffer);
    outFileIndex = (outFileIndex + 1) % numCells;
  }

  // Close all of the output files
  for(outFileIndex = 0; outFileIndex < numCells; ++outFileIndex) {
    fclose(outFiles[outFileIndex]);
  }
  fclose(inFile);
  remove(fileName);

  Mem_Free(outFiles);
  Mem_Free(buffer);
}

/**
 * \brief Returns the number of [grid cells](\ref gridCells) that requested
 *        SOILWAT2 output.
 * 
 * \return an int greater than or equal to 0.
 * 
 * \author Chandler Haukap
 * \date 18 February 2020
 * \ingroup GRID_PRIVATE
 */
static int _getNumberSOILWAT2OutputCells(void) {
  int i, j, count = 0;
  for(i = 0; i < grid_Rows; ++i) {
    for(j = 0; j < grid_Cols; ++j) {
      if(gridCells[i][j].generateSWOutput) count++;
    }
  }
  return count;
}

/**
 * \brief Get an array of cell numbers that requested SOILWAT2 output.
 * 
 * \return An array of ints of size \ref _getNumberSOILWAT2OutputCells().
 * 
 * \sideeffect Allocated a new array.
 * 
 * \author Chandler Haukap
 * \date 21 February 2020
 * \ingroup GRID_PRIVATE
 */
static int* _getSOILWAT2OutputCells(void) {
  int i, j, count = 0;
  int numCells = _getNumberSOILWAT2OutputCells();
  int* cells = Mem_Calloc(numCells, sizeof(int),
                          "_getSOILWAT2OutputCells");
    
  for(i = 0; i < grid_Rows; ++i) {
    for(j = 0; j < grid_Cols; ++j) {
      if(gridCells[i][j].generateSWOutput) {
        cells[count] = (i * grid_Cols) + (j % grid_Cols);
        count++;
      }
    }
  }

  return cells;
}
