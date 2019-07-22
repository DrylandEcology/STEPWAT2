/********************************************************************************/
//  Source file: ST_grid.c
//  Type: module
//  Application: STEPPE - plant community dynamics simulator
//  Purpose: This module handles the grid.
//  History:
//     (5/24/2013) -- INITIAL CODING - DLM
//     (March - July 2019) -- Overhauled by Chandler Haukap with Fredrick Pierson
/********************************************************************************/
/*
 Summary:
    This module handles the gridded mode of STEPWAT2. To accomplish this we use a grid of cells 
    represented by the CellType struct. The entire grid of cells can be referenced by the 
    gridCells variable which is a 2d array of CellTypes. To allow this module to use the same 
    functions as non-gridded mode the CellType structs must be loaded into the global variables 
    using the load_cell function. As long as a cell is loaded in you can be sure that all 
	functions will work as expected. 

    In addition to all of the functionality of non-gridded mode, gridded mode has two additional 
    features: initialization and seed dispersal. Initialization allows vegetation to establish 
    before the simulation experiments begin. Seed dispersal allows each cell to disperse seeds 
    to nearby cells.
*/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include "ST_grid.h"
#include "ST_steppe.h"
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "ST_globals.h"
#include "ST_stats.h"
#include "rands.h"
#include "sxw_funcs.h"
#include "ST_initialization.h"
#include "ST_progressBar.h"
#include "ST_seedDispersal.h"

char sd_Sep;

int grid_Cells;
int UseDisturbances, UseSoils, sd_DoOutput, sd_MakeHeader; //these are treated like booleans

// these are grids to store the SOILWAT variables... also dynamically allocated/freed
SW_SOILWAT *grid_SW_Soilwat, *spinup_SW_Soilwat;
SW_SITE *grid_SW_Site, *spinup_SW_Site;
SW_VEGPROD *grid_SW_VegProd, *spinup_SW_VegProd;
SW_MODEL *grid_SW_Model, *spinup_SW_Model;

// these two variables are used to store the soil/distubance inputs for each grid cell... also dynamically allocated/freed
SoilType *grid_Soils;

Grid_SD_St **grid_SD; //for seed dispersal

// these variables are used for the soil types in the spinup options
int nSoilTypes, *soilTypes_Array, *grid_SoilTypes;

/***************************** Externed variables **********************************/
/* Note that in an ideal world we wouldn't need to extern any variables because 
   every module would declate them in a header file. Hopefully we can get this
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

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
//from ST_species.c
void proportion_Recovery(void);
void save_annual_species_relsize(void);
void copy_species(const SpeciesType* src, SpeciesType* dest);

//from ST_resgroups.c
void rgroup_Grow(void);
void rgroup_Establish(void);
void rgroup_IncrAges(void);
void rgroup_PartResources(void);
void copy_rgroup(const GroupType* src, GroupType* dest);

//from ST_mortality.c
void mort_Main(Bool *killed);
void mort_EndOfYear(void);
void grazing_EndOfYear(void);
void _kill_annuals(void);
void _kill_maxage(void);
void _kill_extra_growth(void);

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

//functions from ST_stats.c
void stat_Collect(Int year);
void stat_Collect_GMort(void);
void stat_Collect_SMort(void);
void stat_Output_AllMorts(void);
void stat_Output_AllBmass(void);
void Output_AllCellAvgBmass(const char * filename);
void stat_Output_Seed_Dispersal(const char * filename, const char sep,
		Bool makeHeader);
void stat_Copy_Accumulators(StatType* newDist, StatType* newPpt, StatType* newTemp, StatType* newGrp, StatType* newGsize, 
                            StatType* newGpr, StatType* newGmort, StatType* newGestab, StatType* newSpp, StatType* newIndv,
                            StatType* newSmort, StatType* newSestab, StatType* newSrecieved, FireStatsType* newGwf, Bool firstTime);
void _make_header_with_std( char *buf);

/* Functions from sxw.c */
SXW_t* getSXW(void);
SXW_resourceType* getSXWResources(void);
transp_t* getTranspWindow(void);
void copy_sxw_variables(SXW_t* newSXW, SXW_resourceType* newSXWResources, transp_t* newTransp_window);

/*********** Locally Used Function Declarations ************/
/***********************************************************/

static void printGeneralInfo(void);
static void _init_grid_files(void);
static void _init_SXW_inputs(Bool init_SW, char *f_roots);
static void allocate_gridCells(int rows, int cols);
static void allocate_accumulators(void);
static void _load_cell(int row, int col, int year, Bool useAccumulators);
static void _read_disturbances_in(void);
static void _read_soils_in(void);
static void _init_soil_layers(int cell, int isSpinup);
static void _read_init_species(void);
static void _read_maxrgroupspecies(void);
static void _read_grid_setup(void);
static void _read_files(void);
static void _init_stepwat_inputs(void);
static void _init_grid_inputs(void);

/******************** Begin Model Code *********************/
/***********************************************************/

/* Print information about the simulation to stdout. */
static void printGeneralInfo(void){
	/* ------------------- Print some general information to stdout ----------------------- */
    printf("Number of iterations: %d\n", SuperGlobals.runModelIterations);
    printf("Number of years: %d\n", SuperGlobals.runModelYears);
	printf("Number of cells: %d\n\n", grid_Cells);
	if(UseDisturbances) printf("Using grid disturbances file\n");
	if(UseSoils) printf("Using grid soils file\n");
	if(initializationMethod == INIT_WITH_SEEDS){
		printf("Running seed dispersal as initialization\n");
	} else if(initializationMethod == INIT_WITH_SPINUP) { 
		printf("Running Spinup as initialization\n");
	}
	if(initializationMethod != INIT_WITH_NOTHING){
		printf("Number of initialization years: %d\n", SuperGlobals.runInitializationYears);
	}
	if(UseSeedDispersal){
		printf("Dispersing seeds between cells\n");
	}
	printf("\n");
	/* --------------------------- END printing general info -------------------------------- */
}

/* Run gridded mode. */
void runGrid(void)
{
	int i, j;
	Bool killedany;
	IntS year, iter;

	_init_grid_files();				// reads in files.in file
	_read_maxrgroupspecies();       // reads in maxrgroupspecies.in file
	_read_grid_setup();             // reads in grid_setup.in file
    _read_files();                  // reads in Stepwat_Inputs/files.in file
    _init_stepwat_inputs();			// reads the stepwat inputs in
	_init_grid_inputs();			// reads the grid inputs in & initializes the global grid variables
	//SWC hist file prefix needs to be cleared
	Mem_Free(SW_Soilwat.hist.file_prefix);
	SW_Soilwat.hist.file_prefix = NULL;

	printGeneralInfo();

	if(initializationMethod != INIT_WITH_NOTHING){
		runInitialization();
	} else {
		/* If no spinup is requested we still need to reset SXW to set the historical weather
		   file names. */
		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
		SXW_Reset(gridCells[0][0].mySXW->f_watin);
		ChDir("..");
	}

	// SOILWAT resets SW_Weather.name_prefix every iteration. This is not the behavior we want 
	// so the name is stored here.
	char SW_prefix_permanent[2048];
	sprintf(SW_prefix_permanent, "%s/%s", grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS], SW_Weather.name_prefix);

	for (iter = 1; iter <= SuperGlobals.runModelIterations; iter++)
	{ //for each iteration

		/*
		 * 06/15/2016 (akt) Added resetting correct historical weather file path,
		 * as it was resetting to original path value (that was not correct for grid version)from input file after every iteration
		 */
		sprintf(SW_Weather.name_prefix, "%s", SW_prefix_permanent); //updates the directory correctly for the weather files so soilwat can find them

		if (BmassFlags.yearly || MortFlags.yearly){
			parm_Initialize();
			_init_grid_inputs();
		}

		// Initialize the plot for each grid cell
		for (i = 0; i < grid_Rows; i++){
			for (j = 0; j < grid_Cols; j++){
				load_cell(i, j);
				Plot_Initialize();
				Globals->currIter = iter;
			}
		}
		unload_cell(); // Reset the global variables

		// If we used spinup we need to reset to the state of the program right after spinup.
		if (initializationMethod == INIT_WITH_SPINUP){
			loadInitializationConditions();
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
			for (i = 0; i < grid_Rows; i++){
				for (j = 0; j < grid_Cols; j++)
				{ //for each cell

					//fprintf(stderr, "year: %d", year);
					load_cell(i, j);

					Globals->currYear = year;

					/* Seed dispersal needs to take into account last year's precipitation, 
					   so we'll record it before calling Env_Generate(). */
					if (year > 1 && UseSeedDispersal){
						gridCells[i][j].mySeedDispersal->lyppt = gridCells[i][j].myEnvironment.ppt;
					}

					/* The following functions mimic ST_main.c. load_cell(i, j) ensures that all global variables reference the specific cell */

					rgroup_Establish(); 		// Establish individuals. Excludes annuals.

					Env_Generate();				// Generated the SOILWAT environment

					rgroup_PartResources();		// Distribute resources
					rgroup_Grow(); 				// Grow

					mort_Main(&killedany); 		// Mortality that occurs during the growing season

					rgroup_IncrAges(); 			// Increment ages of all plants

					grazing_EndOfYear(); 		// Livestock grazing
					
					save_annual_species_relsize(); // Save annuals before we kill them

					mort_EndOfYear(); 			// End of year mortality.

					stat_Collect(year); 		// Update the accumulators

				    _kill_annuals(); 			// Kill annuals
				    _kill_maxage();             // Kill plants that reach max age
					proportion_Recovery(); 		// Recover from any disturbances
					_kill_extra_growth(); 		// Kill superfluous growth

				} /* end model run for this cell*/
			} /* end model run for this row */
			if (UseSeedDispersal)
				disperseSeeds();
			
			unload_cell(); // Reset the global variables
		}/* end model run for this year*/

		// collects the data appropriately for the mort output... (ie. fills the accumulators in ST_stats.c with the values that they need)
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
				SXW_Reset(gridCells[i][j].mySXW->f_watin);
				unload_cell();
			}
		}
		Mem_Free(SW_Soilwat.hist.file_prefix);
		SW_Soilwat.hist.file_prefix = NULL;
		ChDir("..");

	} /*end iterations */

	if(UseProgressBar){
		logProgress(0, 0, OUTPUT);
	}

	// outputs all of the mort and BMass files for each cell...
	for (i = 0; i < grid_Rows; i++){
		for (j = 0; j < grid_Cols; j++)
		{
			int cell = (j + 1) + (i * grid_Cols) - 1;
			load_cell(i, j);
			char fileMort[1024], fileBMass[1024], fileReceivedProb[1024];

			sprintf(fileReceivedProb, "%s%d.csv", grid_files[GRID_FILE_PREFIX_RECEIVEDPROB], cell);
			sprintf(fileMort, "%s%d.csv", grid_files[GRID_FILE_PREFIX_MORTAVG], cell);
			sprintf(fileBMass, "%s%d.csv", grid_files[GRID_FILE_PREFIX_BMASSAVG], cell);
			parm_SetName(fileMort, F_MortAvg);
			parm_SetName(fileBMass, F_BMassAvg);

			if (MortFlags.summary)
				stat_Output_AllMorts();
			if (BmassFlags.summary)
				stat_Output_AllBmass();
			if (UseSeedDispersal && sd_DoOutput)
				stat_Output_Seed_Dispersal(fileReceivedProb, sd_Sep,
						sd_MakeHeader);
		}
	}
	unload_cell(); // Reset the global variables

	//Here creating grid cells avg values output file
	char fileBMassCellAvg[1024];
	sprintf(fileBMassCellAvg, "%s.csv", grid_files[GRID_FILE_PREFIX_BMASSCELLAVG]);
	if (BmassFlags.summary){
		Output_AllCellAvgBmass(fileBMassCellAvg);
	}

	free_grid_memory();	// Free our allocated memory since we do not need it anymore
	parm_free_memory();		// Free memory allocated to the _files array in ST_params.c
	if(initializationMethod == INIT_WITH_SPINUP) {
		freeInitializationMemory();
	}
	logProgress(0, 0, DONE);
}

/* Read the files.in file which was supplied to the program as an argument.
   This function saves the file names it reads to grid_files and grid_directories. */
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
	if (!strcmp("stdout", grid_files[GRID_FILE_LOGFILE]))
		logfp = stdout;
	else
		logfp = OpenFile(grid_files[GRID_FILE_LOGFILE], "w");

	CloseFile(&f);
}

/* Read all gridded mode files excluding grid_setup.in. This function overrided values specified in non-gridded 
   mode files, so make sure you have read the non-gridded mode files before calling this function.*/
static void _init_grid_inputs(void)
{
    int i, j;

	if (UseDisturbances){
		_read_disturbances_in();
	}
	if (UseSeedDispersal || initializationMethod == INIT_WITH_SEEDS)
	{
		initDispersalParameters();
	}
	if(initializationMethod != INIT_WITH_NOTHING){
		_read_init_species();
	}
	if (UseSoils) {
		_read_soils_in();
	}

    for(i = 0; i < grid_Rows; ++i) {
        for(j = 0; j < grid_Cols; ++j) {
            load_cell(i, j);
            gridCells[i][j].DuringInitialization = FALSE;
        }
    }
    unload_cell();
}

/* Read SXW input files */
static void _init_SXW_inputs(Bool init_SW, char *f_roots)
{
	SXW_Init(init_SW, f_roots);	// initializes soilwat
	if (init_SW == TRUE)
	{
		char aString[2048];
		sprintf(aString, "%s/%s", grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS], SW_Weather.name_prefix);
		sprintf(SW_Weather.name_prefix, "%s", aString); //updates the directory correctly for the weather files so soilwat can find them
	}
}

/* Read in the STEPWAT2 files and populate the grid. This only needs to be called once. 
   DEPENDENCYS: gridCells must be allocated first. */
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
			// Set myTanspWindow to the location of the newly allocated transp window
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
	allocate_accumulators();

	ChDir("..");						// go back to the folder we started in
}

/* Reread input files. Be careful because this function reallocates the grid.
   Make sure you call free_grid_memory before calling this function. */
void rereadInputs(void){
    _read_grid_setup();
    _read_files();
    _init_stepwat_inputs();
    _init_grid_inputs();
}

/* Allocates memory for the grid cells. This only needs to be called once. */
static void allocate_gridCells(int rows, int cols){
	int i, j;
	gridCells = (CellType**) Mem_Calloc(rows, sizeof(CellType*), "allocate_gridCells: rows");
	for(i = 0; i < rows; ++i){
		gridCells[i] = (CellType*) Mem_Calloc(cols, sizeof(CellType), "allocate_gridCells: columns");
	}

	/* Allocate all fields specific to gridded mode. This is not necessary for fields like mySpecies 
	   since they are allocated elsewhere in the code.*/
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			// lyr is a dynamically allocated array
			gridCells[i][j].mySoils.lyr = (Grid_Soil_Lyr*) 
				Mem_Calloc(MAX_LAYERS, sizeof(Grid_Soil_Lyr), "allocate_gridCells: mySoils");
			
			// shouldBeInitialized is a dynamically allocated array
			gridCells[i][j].mySpeciesInit.shouldBeInitialized = (int*)
				Mem_Calloc(MAX_SPECIES, sizeof(int), "allocate_gridCells: mySpeciesInit");

			gridCells[i][j].someKillage = (Bool*) Mem_Calloc(1, sizeof(Bool), "allocate_gridCells: someKillage");
		}
	}
}

/* Initialize each gridCell's accumulators. 
   Must be called after STEPWAT inputs have been read. */
static void allocate_accumulators(void){
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
    			gridCells[i][j]._Dist = (StatType*) Mem_Calloc(1, sizeof(StatType), "allocate_accumulators(Dist)");
		    	gridCells[i][j]._Dist->s = (struct accumulators_st *)
		               		Mem_Calloc( SuperGlobals.runModelYears,
		                           		sizeof(struct accumulators_st),
 		                         		"allocate_accumulators(Dist)");
		  	}
 		 	if (BmassFlags.ppt) {
		    	gridCells[i][j]._Ppt = (StatType*) Mem_Calloc(1, sizeof(StatType), "allocate_accumulators(PPT");
		    	gridCells[i][j]._Ppt->s  = (struct accumulators_st *)
		               		Mem_Calloc( SuperGlobals.runModelYears,
		                           		sizeof(struct accumulators_st),
		                          		"allocate_accumulators(PPT)");
 		 	}
		  	if (BmassFlags.tmp) {
		    	gridCells[i][j]._Temp = (StatType*) Mem_Calloc(1, sizeof(StatType), "allocate_accumulators(Temp)");
		    	gridCells[i][j]._Temp->s = (struct accumulators_st *)
		               		Mem_Calloc( SuperGlobals.runModelYears,
 		                          		sizeof(struct accumulators_st),
 		                         		"allocate_accumulators(Temp)");
		  	}
		  	if (BmassFlags.grpb) {
 			   	gridCells[i][j]._Grp = (struct stat_st *)
           				Mem_Calloc( Globals->grpCount,
                		       sizeof(struct stat_st),
                		      "allocate_accumulators(Grp)");
    			ForEachGroup(rg){
      				gridCells[i][j]._Grp[rg].s = (struct accumulators_st *)
             			Mem_Calloc( SuperGlobals.runModelYears,
                         			sizeof(struct accumulators_st),
                        			"allocate_accumulators(Grp[rg].s)");
				}

    			if (BmassFlags.size) {
      				gridCells[i][j]._Gsize = (struct stat_st *)
             				Mem_Calloc( Globals->grpCount,
                         				sizeof(struct stat_st),
                        				"allocate_accumulators(GSize)");
      				ForEachGroup(rg){
          				gridCells[i][j]._Gsize[rg].s = (struct accumulators_st *)
             							Mem_Calloc( SuperGlobals.runModelYears,
                         							sizeof(struct accumulators_st),
                        							"allocate_accumulators(GSize[rg].s)");
					}
    			}
    			if (BmassFlags.pr) {
      				gridCells[i][j]._Gpr = (struct stat_st *)
             				Mem_Calloc( Globals->grpCount,
                		         		sizeof(struct stat_st),
                		        		"allocate_accumulators(Gpr)");
      				ForEachGroup(rg){
          				gridCells[i][j]._Gpr[rg].s = (struct accumulators_st *)
             							Mem_Calloc( SuperGlobals.runModelYears,
                    		     					sizeof(struct accumulators_st),
                        							"allocate_accumulators(Gpr[rg].s)");
					}
    			}

    			if (BmassFlags.wildfire || BmassFlags.prescribedfire) {
      				gridCells[i][j]._Gwf = (struct fire_st *)
             				Mem_Calloc( 1, sizeof(struct fire_st),
                		        		"allocate_accumulators(Gwf)");

      				gridCells[i][j]._Gwf->wildfire = (int *) Mem_Calloc( 1,
                    		  		 sizeof(int) * SuperGlobals.runModelYears,
                    		  		 "allocate_accumulators(Gwf->wildfire)");
      
      				gridCells[i][j]._Gwf->prescribedFire = (int **) Mem_Calloc( 1,
                      						sizeof(int **) * SuperGlobals.max_rgroups,
                       						"allocate_accumulators(Gwf->prescribedfire");

      				ForEachGroup(rg){
        				gridCells[i][j]._Gwf->prescribedFire[rg] = (int *)
          										   Mem_Calloc( SuperGlobals.runModelYears,
                      										   sizeof(int) * SuperGlobals.runModelYears,
                      										   "allocate_accumulators(Gwf->prescribedFire)");
      				}
    			}
  			}

  			if (MortFlags.group) {
    			gridCells[i][j]._Gestab = (struct stat_st *)
             	  		  Mem_Calloc( Globals->grpCount,
                         			sizeof(struct stat_st),
                         			"allocate_accumulators(Gestab)");

				gridCells[i][j]._Gmort = (struct stat_st *)
           		 		 Mem_Calloc( Globals->grpCount,
                       				sizeof(struct stat_st),
                      				"allocate_accumulators(Gmort)");
    			ForEachGroup(rg){
      				gridCells[i][j]._Gestab[rg].s = (struct accumulators_st *)
                     				Mem_Calloc( 1, sizeof(struct accumulators_st),
                                				"allocate_accumulators(Gestab[rg].s)");
					gridCells[i][j]._Gmort[rg].s = (struct accumulators_st *)
           							Mem_Calloc( GrpMaxAge(rg),
                       							sizeof(struct accumulators_st),
                      							"allocate_accumulators(Gmort[rg].s)");
				}
  			}

  			if (BmassFlags.sppb) {
      			gridCells[i][j]._Spp = (struct stat_st *)
               		   Mem_Calloc( Globals->sppCount,
                           		   sizeof(struct stat_st),
                          		   "allocate_accumulators(Spp)");
      			ForEachSpecies(sp){
        			gridCells[i][j]._Spp[sp].s = (struct accumulators_st *)
               			 		  Mem_Calloc( SuperGlobals.runModelYears,
                           					sizeof(struct accumulators_st),
                          					"allocate_accumulators(Spp[sp].s)");
				}

      			if (BmassFlags.indv) {
        			gridCells[i][j]._Indv = (struct stat_st *)
               		 		Mem_Calloc( Globals->sppCount,
                           		 		sizeof(struct stat_st),
                          		 		"allocate_accumulators(Indv)");
        			ForEachSpecies(sp){
          				gridCells[i][j]._Indv[sp].s = (struct accumulators_st *)
               				  		Mem_Calloc( SuperGlobals.runModelYears,
                           						sizeof(struct accumulators_st),
                          						"allocate_accumulators(Indv[sp].s)");
					}
    			}
  			}
  			if (MortFlags.species) {
    			gridCells[i][j]._Sestab = (struct stat_st *)
           					Mem_Calloc( Globals->sppCount,
                       					sizeof(struct stat_st),
                      					"allocate_accumulators(Sestab)");

				gridCells[i][j]._Smort = (struct stat_st *)
           		 		Mem_Calloc( Globals->sppCount,
                       				sizeof(struct stat_st),
                      				"allocate_accumulators(Smort)");

    			ForEachSpecies(sp){
      				gridCells[i][j]._Sestab[sp].s = (struct accumulators_st *)
                    				Mem_Calloc( 1, sizeof(struct accumulators_st),
                                				"allocate_accumulators(Sestab[sp].s)");
					gridCells[i][j]._Smort[sp].s = (struct accumulators_st *)
            		        		Mem_Calloc( SppMaxAge(sp),
            		                    	   sizeof(struct accumulators_st),
             		                   	   	   "allocate_accumulators(Smort[sp].s)");
				}
  			}

  			if (UseSeedDispersal) {
	  			gridCells[i][j]._Sreceived = Mem_Calloc( Globals->sppCount, sizeof(struct stat_st), "allocate_accumulators(Sreceived)");

	  			ForEachSpecies(sp) {
		  			gridCells[i][j]._Sreceived[sp].s = (struct accumulators_st *)
					  									Mem_Calloc( SuperGlobals.runModelYears,
														  		   sizeof(struct accumulators_st), 
																   "allocate_accumulators(Sreceived[sp].s)");
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

/* Free all memory allocated to the gridded mode during initialization. */
void free_grid_memory(void)
{
	//frees all the memory allocated in this file ST_Grid.c (most of it is dynamically allocated in _init_grid_globals() & _load_grid_globals() functions)
	int i, j, sd_i;
	SppIndex s;

	/* Free memory that we have allocated in ST_grid.c */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			/* Use deallocate_Globals from ST_main to deallocate global variables,
			   and free_all_sxw_memory from sxw to deallocate SXW variables. */
			load_cell(i,j);
			deallocate_Globals(TRUE);
			free_all_sxw_memory();
			// If seed dispersal is on we allocated additional memory
			if(UseSeedDispersal) {
				ForEachSpecies(s) {
					for(sd_i = 0; sd_i < grid_Rows; ++sd_i){
						Mem_Free(gridCells[i][j].mySeedDispersal[s].probabilityOfDispersing[sd_i]);
					}
					Mem_Free(gridCells[i][j].mySeedDispersal[s].probabilityOfDispersing);
				}
			}
			unload_cell();

			Mem_Free(gridCells[i][j].mySpeciesInit.shouldBeInitialized);
			Mem_Free(gridCells[i][j].mySeedDispersal);
			Mem_Free(gridCells[i][j].someKillage);
			Mem_Free(gridCells[i][j].mySoils.lyr);
		}
	}

	for(i = 0; i < grid_Rows; ++i){
		Mem_Free(gridCells[i]);
	}
	Mem_Free(gridCells);
}

/*
static void _load_cell(int row, int col, int year, Bool useAccumulators)
{
	// loads the specified cell into the global variables

	int cell = col + ((row - 1) * grid_Cols) - 1; // converts the row/col into an array index
	int j, k;
	GrpIndex c;
	SppIndex s;
	//fprintf(stderr, " loading cell: %d; ", cell);
	if (useAccumulators)
		stat_Load_Accumulators(cell, year);

	ForEachSpecies(s)
	{
		if (!Species[s]->use_me)
			continue;

		Mem_Free(Species[s]->kills);
		Mem_Free(Species[s]->seedprod);
		_free_head(Species[s]->IndvHead); //free_head() frees the memory allocated by the head and the memory allocated by each part of the linked list

		*Species[s] = grid_Species[s][cell];

		Species[s]->kills = Mem_Calloc(grid_Species[s][cell].max_age,
				sizeof(IntUS), "_load_cell(Species[s]->kills)");
		Species[s]->seedprod = Mem_Calloc(grid_Species[s][cell].viable_yrs,
				sizeof(RealF), "_load_cell(Species[s]->seedprod)");

		memcpy(Species[s]->kills, grid_Species[s][cell].kills,
				grid_Species[s][cell].max_age * sizeof(IntUS));
		memcpy(Species[s]->seedprod, grid_Species[s][cell].seedprod,
				grid_Species[s][cell].viable_yrs * sizeof(RealF));
		Species[s]->IndvHead = _copy_head(grid_Species[s][cell].IndvHead); //copy_head() deep copies the linked list structure (allocating memory when needed)... it will even allocate memory for the head of the list

	}

	ForEachGroup(c)
	{
		if (!RGroup[c]->use_me)
			continue;
		Mem_Free(RGroup[c]->kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)... we need to free it then reallocate it then memcpy it to get the deep copy we want

		*RGroup[c] = grid_RGroup[c][cell]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (ie copy the values in the pointers instead of the addresses) the pointers.  A shallow copy will copy over the values for every non-pointer (C itself does not inherently know how to deep copy, so we must code this behaviour).

		RGroup[c]->kills = Mem_Calloc(grid_RGroup[c][cell].max_age,
				sizeof(IntUS), "_load_cell(RGroup[c]->kills");
		memcpy(RGroup[c]->kills, grid_RGroup[c][cell].kills,
				grid_RGroup[c][cell].max_age * sizeof(IntUS));
	}

	*Succulent = grid_Succulent[cell];
	*Env = grid_Env[cell];
	*Plot = grid_Plot[cell];
	*Globals = grid_Globals[cell];

	Mem_Free(SXW->f_roots);
	Mem_Free(SXW->f_phen);
	Mem_Free(SXW->f_bvt);
	Mem_Free(SXW->f_prod);
	Mem_Free(SXW->f_watin);
	Mem_Free(SXW->transpTotal);
	ForEachVegType(k) {
		Mem_Free(SXW->transpVeg[k]);
	}

	if (SXW->swc != NULL)
		Mem_Free(SXW->swc);
	for (j = 0; j < SW_Site.n_layers + SW_Site.deepdrain; j++)
		Mem_Free(SW_Site.lyr[j]);
	Mem_Free(SW_Site.lyr);

	*SXW = grid_SXW[cell];
	SW_Site = grid_SW_Site[cell];
	SW_Soilwat = grid_SW_Soilwat[cell];
	SW_VegProd = grid_SW_VegProd[cell];

	SXW->transpTotal = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealD),
			"_load_cell(SXW->transp)");
	ForEachVegType(k) {
		SXW->transpVeg[k] = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealD),
				"_load_cell(SXW->transp)");
	}
	SXW->swc = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealF),
			"_load_cell(SXW->swc)");

	SXW->f_roots = Str_Dup(grid_SXW[cell].f_roots);
	SXW->f_phen = Str_Dup(grid_SXW[cell].f_phen);
	SXW->f_bvt = Str_Dup(grid_SXW[cell].f_bvt);
	SXW->f_prod = Str_Dup(grid_SXW[cell].f_prod);
	SXW->f_watin = Str_Dup(grid_SXW[cell].f_watin);
	memcpy(SXW->transpTotal, grid_SXW[cell].transpTotal,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	ForEachVegType(k) {
		memcpy(SXW->transpVeg[k], grid_SXW[cell].transpVeg[k],
				SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	}
	memcpy(SXW->swc, grid_SXW[cell].swc,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealF));

	SW_Site.lyr = Mem_Calloc(
			grid_SW_Site[cell].n_layers + grid_SW_Site[cell].deepdrain,
			sizeof(SW_LAYER_INFO *), "_load_cell(SW_Site.lyr)");
	for (j = 0;
			j < grid_SW_Site[cell].n_layers + grid_SW_Site[cell].deepdrain;
			j++)
	{
		SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
				"_load_cell(SW_Site.lyr[j])");
		memcpy(SW_Site.lyr[j], grid_SW_Site[cell].lyr[j],
				sizeof(SW_LAYER_INFO));
	}
}
*/

/* Load gridCells[row][col] into the globals variables.
   Any call to this function should have an accompanying call to unload_cell(). */
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
	DuringInitialization = gridCells[row][col].DuringInitialization;

	_SomeKillage = gridCells[row][col].someKillage;

	/* Copy this cell's accumulators into the local accumulators in ST_stats.c */
	stat_Copy_Accumulators(gridCells[row][col]._Dist, gridCells[row][col]._Ppt, gridCells[row][col]._Temp,
	                       gridCells[row][col]._Grp, gridCells[row][col]._Gsize, gridCells[row][col]._Gpr,
						   gridCells[row][col]._Gmort, gridCells[row][col]._Gestab, gridCells[row][col]._Spp,
						   gridCells[row][col]._Indv, gridCells[row][col]._Smort, gridCells[row][col]._Sestab,
						   gridCells[row][col]._Sreceived, gridCells[row][col]._Gwf, gridCells[row][col].stats_init);

	/* Copy this cell's SXW variables into the local variables in sxw.c */
	copy_sxw_variables(gridCells[row][col].mySXW, gridCells[row][col].mySXWResources, gridCells[row][col].myTranspWindow);
}

/* Nullify all global variables. This function should appear after every call 
   to load_cell to prevent accidental modification of a grid cell.
   Example usage:
   for i in rows{
       for j in columns{
	       load_cell(i, j)
		   //additional functions
	   }
   }
   unload_cell() */
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

/**************************************************************/
static Bool GetALine2(FILE *f, char buf[], int limit)
{
	//this is similar to the getaline function in filefuncs.c, except this one checks for carriage return characters and doesn't deal with whitespace/... (since excel writes them into .csv files for some aggravating reason)... this one is probably less efficient overall though.
	//this is treating '\r', '\n', and '\r\n' all like they are valid line feed characters... in reality '\r' by itself should never be, but it's not like I can stop excel from outputting .csv files however the heck it feels like...
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

/* Reads the grid disturbance CSV. This function will override disturbance inputs from non-gridded mode. */
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
            num = sscanf(buf, "%d,%d,%d,%d,%hu,%hu,%f,%hu,%hu,%hu,%f", &cell,
                &Globals->pat.use, &Globals->mound.use, &Globals->burrow.use, &RGroup[rg]->killyr, 
				&RGroup[rg]->killfreq_startyr, &RGroup[rg]->killfreq, &RGroup[rg]->extirp, 
				&RGroup[rg]->grazingfrq, &RGroup[rg]->grazingfreq_startyr, &RGroup[rg]->ignition);
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

/***********************************************************/
static int _get_value_index(char* s, char seperator, int nSeperators)
{
	//pretty much this function goes through s until it find nSeperators worth of seperators and then returns the index of the next character
	//this is used to do most of the parsing in the _read_soils_in() function
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

static int _get_value_index2(char* s, int nSeperators)
{
	//pretty much this function goes through s until it find nSeperators worth of seperators and then returns the index of the next character
	//this is used to do most of the parsing in the _read_soils_in() function
	int i = 0, sep = 0;
	while (1)
	{
		i++;
		if (*s++ == '\0') //checks if the current char equals the seperator... and then increments the char pointer
			if (++sep == nSeperators) //needs ++sep, not sep++
				break;
	}
	return i;
}

/* Read the grid soils input */
static void _read_soils_in(void)
{
	FILE *f;
	char buf[4096];
	char rootsin[20];
	char seps[] = ",";
	char *token;
	int i, j, k, cell, num, do_copy, copy_cell, num_layers, depth, depthMin,
			stringIndex, row, col, copy_cell_row, copy_cell_col;
	float d[11];

	nSoilTypes = 0; //initialize our soil types counter

	f = OpenFile(grid_files[GRID_FILE_SOILS], "r");

	GetALine2(f, buf, 4096); // gets rid of the first line (since it just defines the columns)... it's only there for user readability
	for (i = 0; i < grid_Cells; i++)
	{
	    row = i / grid_Cols;
	    col = i % grid_Cols;

	    load_cell(row, col);

		if (!GetALine2(f, buf, 4096))
			break;
		gridCells[row][col].mySoils.rootsFile[0] = '\0';
		rootsin[0] = '\0';

		num = 0;
		token = strtok(buf, seps);
		while (token != NULL && num < 5)
		{
			switch (num)
			{
			case 0:
				cell = atoi(token);
				break;
			case 1:
				do_copy = atoi(token);
				break;
			case 2:
				copy_cell = atoi(token);
				copy_cell_row = copy_cell / grid_Cols;
				copy_cell_col = copy_cell % grid_Cols;
				break;
			case 3:
				num_layers = atoi(token);
				break;
			case 4:
				strcpy(rootsin, token);
				break;
			}
			num++;
			if (num != 5)
				token = strtok(NULL, seps);
		}

		if (num < 5)
			if (!do_copy)
				LogError(logfp, LOGFATAL, "Invalid %s file", grid_files[GRID_FILE_SOILS]);
		if (!do_copy)
			stringIndex = _get_value_index2(buf, 5); //gets us the index of the string that is right after what we just parsed in

		if (num_layers > MAX_LAYERS)
			LogError(logfp, LOGFATAL,
					"Invalid %s file line %d num_layers (%d) exceeds MAX_LAYERS (%d)",
					grid_files[GRID_FILE_SOILS], i + 2, num_layers, MAX_LAYERS);

		if (do_copy == 1 && copy_cell > -1 && copy_cell < grid_Cells
				&& cell != 0 && copy_cell < cell)
		{ //copy this cells values from a previous cell's
			gridCells[row][col].mySoils.lyr = Mem_Calloc(gridCells[copy_cell_row][copy_cell_col].mySoils.num_layers,
					sizeof(Grid_Soil_Lyr), "_read_soils_in()");
			for (j = 0; j < gridCells[copy_cell_row][copy_cell_col].mySoils.num_layers; j++)
				gridCells[row][col].mySoils.lyr[j] = gridCells[copy_cell_row][copy_cell_col].mySoils.lyr[j];
			gridCells[row][col].mySoils.num_layers = gridCells[copy_cell_row][copy_cell_col].mySoils.num_layers;

            // TODO: figure out how to change this code for our new struct
			// grid_SoilTypes[cell] = grid_SoilTypes[copy_cell];

			strcpy(gridCells[row][col].mySoils.rootsFile, gridCells[copy_cell_row][copy_cell_col].mySoils.rootsFile);

			continue;
		}
		else if (do_copy == 1){
			LogError(logfp, LOGFATAL,
					"Invalid %s file line %d invalid copy_cell attempt",
					grid_files[GRID_FILE_SOILS], i + 2);
		}

        // TODO: figure out how to change this code for our new struct
		// grid_SoilTypes[cell] = nSoilTypes;
		// soilTypes_Array[nSoilTypes] = cell;
		nSoilTypes++;

		depthMin = 0;
		gridCells[row][col].mySoils.num_layers = num_layers;
		strcpy(gridCells[row][col].mySoils.rootsFile, rootsin);
		gridCells[row][col].mySoils.lyr = Mem_Calloc(num_layers, sizeof(Grid_Soil_Lyr),
				"_read_soils_in()");
		for (j = 0; j < num_layers; j++)
		{
			//the idea behind using &buf[stringIndex] is that we start scanning at the point in the string that is right after what we just parsed... the & is there because we have to send sscanf the pointer that points to that location
			num = sscanf(&buf[stringIndex],
					"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &depth, &d[0], &d[1],
					&d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9],
					&d[10]);
			if (num != 12)
				LogError(logfp, LOGFATAL,
						"Invalid '%s' file line %d invalid soil layer input",
						grid_files[GRID_FILE_SOILS], i + 2);

			k = stringIndex;
			stringIndex += _get_value_index(&buf[stringIndex], ',', 12); //updates the index of the string that we are at
			if (k == stringIndex)
				LogError(logfp, LOGFATAL,
						"Invalid %s file line %d not enough soil layers",
						grid_files[GRID_FILE_SOILS], i + 2);

			for (k = 0; k < 11; k++)
				gridCells[row][col].mySoils.lyr[j].data[k] = d[k];
			gridCells[row][col].mySoils.lyr[j].width = depth - depthMin;
			depthMin = depth;
		}
	}

	if (i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file, not enough cells",
				grid_files[GRID_FILE_SOILS]);

    unload_cell();
	CloseFile(&f);
}

/*
static void _init_soil_layers(int cell, int isSpinup)
{
	// initializes the soilwat soil layers for the cell correctly based upon the input gathered from our grid_soils input file
	// pretty much takes the data from grid_Soils (read in in _read_soils_in()) and converts it to what SW_Site needs...
	// this function does generally the same things that the _read_layers() function in SW_Site.c does, except that it does it in a way that lets us use it in the grid...
	int i, j;
	i = cell;
	char* errtype;
	Bool evap_ok = TRUE, transp_ok_forb = TRUE, transp_ok_tree = TRUE,
			transp_ok_shrub = TRUE, transp_ok_grass = TRUE; //mitigate gaps in layers
	Bool fail = FALSE;
	RealF fval = 0;

	if (SW_Site.deepdrain)
		SW_Site.n_layers++;
	for (j = 0; j < SW_Site.n_layers; j++)
		Mem_Free(SW_Site.lyr[j]);
	Mem_Free(SW_Site.lyr);

	SW_Site.n_layers = grid_Soils[i].num_layers;
	SW_Site.n_evap_lyrs = SW_Site.n_transp_lyrs[SW_FORBS] =
			SW_Site.n_transp_lyrs[SW_TREES] = SW_Site.n_transp_lyrs[SW_SHRUB] =
					SW_Site.n_transp_lyrs[SW_GRASS] = 0;

	SW_Site.lyr = Mem_Calloc(SW_Site.n_layers + SW_Site.deepdrain,
			sizeof(SW_LAYER_INFO *), "_init_grid_globals()");
	for (j = 0; j < SW_Site.n_layers; j++)
	{
		if (LT(grid_Soils[i].lyr[j].data[0], 0.))
		{
			fail = TRUE;
			fval = grid_Soils[i].lyr[j].data[0];
			errtype = Str_Dup("bulk density");
		}
		else if (LT(grid_Soils[i].lyr[j].data[1],
				0.) || GT(grid_Soils[i].lyr[j].data[1], 1.0))
		{
			fail = TRUE;
			fval = grid_Soils[i].lyr[j].data[1];
			errtype = Str_Dup("gravel content");
		}
		else if (LE(grid_Soils[i].lyr[j].data[7], 0.))
		{
			fail = TRUE;
			fval = grid_Soils[i].lyr[j].data[7];
			errtype = Str_Dup("sand proportion");
		}
		else if (LE(grid_Soils[i].lyr[j].data[8], 0.))
		{
			fail = TRUE;
			fval = grid_Soils[i].lyr[j].data[8];
			errtype = Str_Dup("clay proportion");
		}
		else if (LT(grid_Soils[i].lyr[j].data[9], 0.))
		{
			fail = TRUE;
			fval = grid_Soils[i].lyr[j].data[9];
			errtype = Str_Dup("impermeability");
		}
		if (fail)
		{
			LogError(logfp, LOGFATAL, "Invalid %s (%5.4f) in layer %d.\n",
					errtype, fval, j + 1);
		}

		SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
				"_init_grid_globals()");

		//indexes (for grid_Soils[i].lyr[j].data):
		//0		   1				2		3	  		4			5			6			7	   8		9		10
		//matricd	gravel_content  evco  	trco_grass  trco_shrub  trco_tree  	trco_forb	%sand  %clay imperm soiltemp
		SW_Site.lyr[j]->width = grid_Soils[i].lyr[j].width;
		SW_Site.lyr[j]->soilMatric_density = grid_Soils[i].lyr[j].data[0];
		SW_Site.lyr[j]->fractionVolBulk_gravel = grid_Soils[i].lyr[j].data[1];
		SW_Site.lyr[j]->evap_coeff = grid_Soils[i].lyr[j].data[2];
		SW_Site.lyr[j]->transp_coeff[3] = grid_Soils[i].lyr[j].data[3];
		SW_Site.lyr[j]->transp_coeff[1] = grid_Soils[i].lyr[j].data[4];
		SW_Site.lyr[j]->transp_coeff[0] = grid_Soils[i].lyr[j].data[5];
		SW_Site.lyr[j]->transp_coeff[2] = grid_Soils[i].lyr[j].data[6];
		SW_Site.lyr[j]->fractionWeightMatric_sand =
				grid_Soils[i].lyr[j].data[7];
		SW_Site.lyr[j]->fractionWeightMatric_clay =
				grid_Soils[i].lyr[j].data[8];
		SW_Site.lyr[j]->impermeability = grid_Soils[i].lyr[j].data[9];
		SW_Site.lyr[j]->my_transp_rgn[0] = 0;
		SW_Site.lyr[j]->my_transp_rgn[2] = 0;
		SW_Site.lyr[j]->my_transp_rgn[1] = 0;
		SW_Site.lyr[j]->my_transp_rgn[3] = 0;
		SW_Site.lyr[j]->sTemp = grid_Soils[i].lyr[j].data[10];

		if (evap_ok)
		{
			if (GT(SW_Site.lyr[j]->evap_coeff, 0.0))
				SW_Site.n_evap_lyrs++;
			else
				evap_ok = FALSE;
		}
		if (transp_ok_tree)
		{
			if (GT(SW_Site.lyr[j]->transp_coeff[0], 0.0))
				SW_Site.n_transp_lyrs[SW_TREES]++;
			else
				transp_ok_tree = FALSE;
		}
		if (transp_ok_shrub)
		{
			if (GT(SW_Site.lyr[j]->transp_coeff[1], 0.0))
				SW_Site.n_transp_lyrs[SW_SHRUB]++;
			else
				transp_ok_shrub = FALSE;
		}
		if (transp_ok_grass)
		{
			if (GT(SW_Site.lyr[j]->transp_coeff[3], 0.0))
				SW_Site.n_transp_lyrs[SW_GRASS]++;
			else
				transp_ok_grass = FALSE;
		}
		if (transp_ok_forb)
		{
			if (GT(SW_Site.lyr[j]->transp_coeff[2], 0.0))
				SW_Site.n_transp_lyrs[SW_FORBS]++;
			else
				transp_ok_forb = FALSE;
		}
		water_eqn(SW_Site.lyr[j]->fractionVolBulk_gravel,
				SW_Site.lyr[j]->fractionWeightMatric_sand,
				SW_Site.lyr[j]->fractionWeightMatric_clay, j); //in SW_Site.c, called to initialize some layer data...
		SW_Site.lyr[j]->swcBulk_fieldcap = SW_SWPmatric2VWCBulk(
				SW_Site.lyr[j]->fractionVolBulk_gravel, 0.333, j)
				* SW_Site.lyr[j]->width;
		SW_Site.lyr[j]->swcBulk_wiltpt = SW_SWPmatric2VWCBulk(
				SW_Site.lyr[j]->fractionVolBulk_gravel, 15, j)
				* SW_Site.lyr[j]->width;
		//From calculate_soilBulkDensity in SW_Site.c
		SW_Site.lyr[j]->soilBulk_density = SW_Site.lyr[j]->soilMatric_density
				* (1 - SW_Site.lyr[j]->fractionVolBulk_gravel)
				+ (SW_Site.lyr[j]->fractionVolBulk_gravel * 2.65);
		//already checked for max_layers condition
	}
	if (SW_Site.deepdrain)
	{
		SW_Site.n_layers++;
		SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
				"_init_grid_globals()");
		SW_Site.lyr[j]->width = 1.0;
	}

	init_site_info(); //in SW_Site.c, called to initialize layer data...

	_init_SXW_inputs(FALSE, grid_Soils[i].rootsFile); //we call this so that SXW can set the correct sizes/values up for the memory dynamically allocated in sxw.c
}
*/

/* Read the species initialization CSV. This function only needs to be called if the user requests initialization.*/
static void _read_init_species(void)
{
	FILE *f;
	int i, j, num, cell, do_copy, copy_cell, useInitialization, seeds_Avail,
	    row, col, copy_cell_row, copy_cell_col;
	Bool isAnyCellOnForSpinup = FALSE;
	char buf[4096];

	//open the file/do the reading
	f = OpenFile(grid_files[GRID_FILE_INIT_SPECIES], "r");

	GetALine2(f, buf, 4096); // gets rid of the first line (since it just defines the columns)... it's only there for user readability
	for (i = 0; i < grid_Cells; i++)
	{
		useInitialization = FALSE;
	    row = i / grid_Cols;
	    col = i % grid_Cols;

	    load_cell(row, col);

		if (!GetALine2(f, buf, 4096))
			break;

		num = sscanf(buf, "%d,%d,%d", &cell, &do_copy, &copy_cell);

		copy_cell_row = copy_cell / grid_Cols;
		copy_cell_col = copy_cell % grid_Cols;

		if (num != 3)
			LogError(logfp, LOGFATAL, "Invalid %s file", grid_files[GRID_FILE_INIT_SPECIES]);

		int stringIndex = _get_value_index(buf, ',', 3); //gets us the index of the string that is right after what we just parsed in

		if (do_copy == 1 && copy_cell > -1 && copy_cell < grid_Cells
				&& cell != 0 && copy_cell < cell)
		{ //copy this cells values from a previous cell's
			for (j = 0; j < Globals->sppCount; j++)
				gridCells[row][col].mySpeciesInit.shouldBeInitialized[j] =
				    gridCells[copy_cell_row][copy_cell_col].mySpeciesInit.shouldBeInitialized[j];
			gridCells[row][col].mySpeciesInit.useInitialization =
			    gridCells[copy_cell_row][copy_cell_col].mySpeciesInit.useInitialization;
			continue;
		}
		else if (do_copy == 1)
			LogError(logfp, LOGFATAL,
					"Invalid %s file line %d invalid copy_cell attempt",
					grid_files[GRID_FILE_INIT_SPECIES], i + 2);

		//going through each species
		SppIndex s;
		ForEachSpecies(s)
		{
			num = sscanf(&buf[stringIndex], "%d,", &seeds_Avail);
			if (num != 1){
				LogError(logfp, LOGFATAL, "Invalid %s file line %d invalid species input",
						 grid_files[GRID_FILE_INIT_SPECIES], i + 2);
			}

			if(seeds_Avail){
				useInitialization = TRUE;
				isAnyCellOnForSpinup = TRUE;
			}

			gridCells[row][col].mySpeciesInit.shouldBeInitialized[s] = seeds_Avail;
			stringIndex += _get_value_index(&buf[stringIndex], ',', 1);
		}

		gridCells[row][col].mySpeciesInit.useInitialization = useInitialization;
	}

	if (i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file, not enough cells",
				grid_files[GRID_FILE_INIT_SPECIES]);

	if(!isAnyCellOnForSpinup){
		LogError(logfp, LOGWARN, "Initialization is on, but no species are turned on for initialization inside %s.",
				grid_files[GRID_FILE_INIT_SPECIES]);
	}

    unload_cell();
	CloseFile(&f);
}

/* Read the maxrgroupspecies file. */
static void _read_maxrgroupspecies(void)
{
    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    parm_SetName(grid_files[GRID_FILE_MAXRGROUPSPECIES], F_MaxRGroupSpecies);
    maxrgroupspecies_init();
    ChDir("..");
}

/* Read the non-gridded mode files.in file. */
static void _read_files(void)
{
    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    files_init();
    ChDir("..");
}

/* Reads the grid setup file and allocates gridCells.*/
static void _read_grid_setup(void)
{
    FILE *f;
    char buf[1024], initializationType[1024];
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
    allocate_gridCells(grid_Rows, grid_Cols);

    GetALine(f, buf);
    i = sscanf(buf, "%d", &UseDisturbances);
    if (i != 1)
        LogError(logfp, LOGFATAL,
                 "Invalid grid setup file (disturbances line wrong)");

    GetALine(f, buf);
    i = sscanf(buf, "%d", &UseSoils);
    if (i != 1)
        LogError(logfp, LOGFATAL, "Invalid grid setup file (soils line wrong)");

    GetALine(f, buf);
    i = sscanf(buf, "%d", &j);
    if (i != 1)
        LogError(logfp, LOGFATAL,
                 "Invalid grid setup file (seed dispersal line wrong)");
    UseSeedDispersal = itob(j);

	GetALine(f, buf);
	i = sscanf(buf, "%s", initializationType);
	if(i < 1){
		LogError(logfp, LOGFATAL, "Invalid grid setup file (Initialization line wrong)");
	}
	if(!strncmp(initializationType, "spinup", 6)){
		initializationMethod = INIT_WITH_SPINUP;
	} else if(!strncmp(initializationType, "seeds", 5)){
		initializationMethod = INIT_WITH_SEEDS;
	} else if(!strncmp(initializationType, "none", 4)){
		initializationMethod = INIT_WITH_NOTHING;
	} else {
		LogError(logfp, LOGFATAL, 
		         "Invalid grid setup file (Initialization line wrong. Valid options are \"spinup\", \"seeds\", or \"none\")");
	}

	if(initializationMethod != INIT_WITH_NOTHING){
		GetALine(f, buf);
		i = sscanf(buf, "%hd", &SuperGlobals.runInitializationYears);
		if(i < 1){
			LogError(logfp, LOGFATAL, "Invalid grid setup file (Initialization years line wrong)");
		}
	}

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_DoOutput) != 1)
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal output line\n", grid_files[GRID_FILE_SETUP]);

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_MakeHeader) != 1)
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal make header line\n",
				grid_files[GRID_FILE_SETUP]);

	GetALine(f, buf);
	if (sscanf(buf, "%c", &sd_Sep) != 1)
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal seperator line\n",
				grid_files[GRID_FILE_SETUP]);

	if (sd_Sep == 't') //dealing with tab and space special cases...
		sd_Sep = '\t';
	else if (sd_Sep == 's')
		sd_Sep = ' ';

    CloseFile(&f);
}

/* Output a master .csv file containing the averages across all cells. */
void Output_AllCellAvgBmass(const char * filename){
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

	FILE* file;
	file = fopen(filename, "w");

	buf[0]='\0';

	load_cell(0, 0);

	if(BmassFlags.header){
		_make_header_with_std(buf);
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

		/* ---------------- Average all accumulators ---------------- */
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

		fprintf(file, "%s\n", buf); // Finally, print this line
	} // End for each year

	unload_cell();
	fclose(file); // Close the file
}