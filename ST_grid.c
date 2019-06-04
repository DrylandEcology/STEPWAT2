/********************************************************/
//  Source file: ST_grid.c
//  Type: module
//  Application: STEPPE - plant community dynamics simulator
//  Purpose: This module handles the grid.
//  History:
//     (5/24/2013) -- INITIAL CODING - DLM
//
//	WARNING: This module deals with a LARGE amount of dynamic memory allocation/deallocation (there can be potentially hundreds of thousands of allocations/frees called by this code depending on settings ran).
//			 Be very wary when editing it as even small changes could possibly cause massive memory errors/leaks to occur.  In particular be careful when copy/freeing the linked list of individuals (I would suggest just using the code I wrote for this as it works and it's pretty easy to screw it up).
//			 Always keep in mind that for every time memory is allocated (every time alloc or calloc is called), a corresponding free is required.  I would recommend valgrind or a similar program to debug memory errors as doing it without some kind of tool would be crazy.
/********************************************************/
/*

 ----------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
 (DLM) : 7-16-2012 : General notes about this module and why it does certain things (put in the beginning so that they'd be seen)
 ----------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------

 ----------------------------------------------------------------------------------------------------------------
 the general idea when copying over the dynamically allocated data is (ie. how to deep copy a structure):
 ----------------------------------------------------------------------------------------------------------------

 1.) Free the dynamically allocated memory (use the _free_head() function to free a linked list of individuals)
 2.) Shallow copy the data (ie. *Species[s] = grid_Species[s][cell])... this will correctly copy all of the data of the structure that isn't a pointer.  For the pointers it will simply copy the address (hence why it is a shallow copy).
 3.) Allocate the appropriate amount of memory for the pointers that are being copied to
 4.) Use memcpy to copy the data over to your newly allocated pointer (or use the _copy_head() function to copy the linked list of individuals)
 5.) Be careful at all stages of this process as it is easy to make a simple error that can be very aggravating to try and track down.

 ----------------------------------------------------------------------------------------------------------------
 explanation of why all this trouble is gone through when copying the memory (as it is not obvious by any means):
 ----------------------------------------------------------------------------------------------------------------

 First off, the sizes of all the dynamically allocated memory can be different for every grid cell.  To deal with this we free and reallocate memory every time we copy a grid cell.
 This also gets around the problem that doing a shallow copy via something like "*Species[s] = grid_Species[s][cell]" would overwrite the address of all the pointers contained within and cause a memory leak if they had already been allocated memory.
 The only other way to copy the data would be to individually copy each data member of each structure which would lead to both a much larger code size and code that would break every time a new member is added to any structure in the program while providing likely no performance difference.
 It is also important to know what the difference between a shallow copy and a deep copy is.  When we are copying over all of these variables, what we want is a deep copy (never to be confused with a shallow copy).
 A shallow copy of a structure is when all of the data members are copied over.  However, the pointers are treated specially.  Instead of copying the data contained within the pointers directly it simply copies the address.
 This leaves two pointers that point to the same thing (not what we want) and if the pointer that is copied over was pointing to any memory it would subsequently cause a memory leak.
 A deep copy of a structure will copy over all of the data members and copy the memory that the pointers are pointing to (not the addresses).  It would be great if C inherently knew how to do this, but it does not.
 C's copy constructor cannot know how to do this because it cannot know the length of the pointers.  Also, it would be potentially confusing as it would have to allocate memory which would be unexpected.
 We must code this behavior in so that we copy over the values held by the pointers into the copy's own separate memory.

 ----------------------------------------------------------------------------------------------------------------
 about performance concerns:
 ----------------------------------------------------------------------------------------------------------------

 This module (ST_grid.c) allocates/deallocates a very large amount of memory.  It does so because it must in order to successfully run as desired.
 The performance hit of all this memory management is surprisingly low in CPU execution time.  Luckily, modern day implementations of malloc/free/memcpy are very fast.
 After profiling the code the time spent allocating/deallocating/copying memory is completely negligible compared to the time spent doing calculations.
 Where the approach has it's downsides is that the program requires a TON of memory in order to do large simulations (ie. it took around 2.8 GB for a 10,000 cell grid when I tried it).
 This shouldn't be an issue in most cases.  It is unavoidable though that at some point the number of cells in a simulation will be bounded by the amount of memory available.
 Issues could possibly arise if you're trying to run a simulation that requires more memory then your system has available.  I don't know of a way to easily check for that condition, so just don't do it.

 ----------------------------------------------------------------------------------------------------------------
 If any of the concepts I have been discussing seem confusing (or your knowledge of pointers feels rusty) I would suggest brushing up on your pointers/memory management.
 Some things to go over would be correct free/alloc/memcpy usage (keep in mind that a free is needed for every corresponding alloc call, some people seem not to comprehend that a pointer of pointers (ie. int**) must be freed in multiple steps, otherwise you lose memory), pointer arithmetic, and the difference between arrays & pointers in C.
 ----------------------------------------------------------------------------------------------------------------
 (AKT) : 9-7-2015 : Added extra Grid Cell Avg Output file for biomass values
 */
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include "ST_steppe.h"
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "ST_globals.h"
#include "ST_stats.h"
#include "rands.h"

#include "sxw_funcs.h"
#include "sxw.h"
#include "sxw_vars.h"

#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Model.h"
#include "SW_Weather.h"
#include "sw_src/pcg/pcg_basic.h"

/***************** Structure Declarations ******************/
/***********************************************************/
// represents a single soil layer
struct _grid_soil_lyr_st
{ 
	// Data for this soil layer
	float data[11];
	// Vertical width of this layer
	int width;
}typedef Grid_Soil_Lyr;

//represents the input data for all the soil layers of a cell
struct Soil_st
{
	// Number of soil layers (size of lyr array)
	int num_layers;
	// Name of the roots file belonging to this cell
	char rootsFile[20];
	// Specific layer's information
	Grid_Soil_Lyr* lyr;
}typedef SoilType;


struct _grid_disturb_st
{
	int choices[3]; //used as boolean values (ie flags as to whether or not to use the specified disturbance)
	int kill_yr, /* kill the group in this year; if 0, don't kill, but see killfreq */
	    killfreq_startyr,/* start year for kill frequency*/
	    killfrq, /* kill group at this frequency: <1=prob, >1=# years */
	    extirp, /* year in which group is extirpated (0==ignore) */
	    veg_prod_type, /* type of VegProd.  1 for tree, 2 for shrub, 3 for grass, 4 for forb */
		grazingfreq_startyr,/* start year for grazing frequency*/
	    grazing_frq; /* grazing effect on group at this frequency: <1=prob, >1=# years */

	RealF xgrow; /* ephemeral growth = mm extra ppt * xgrow */

}typedef Grid_Disturb_St;

struct _grid_sd_struct
{ //for seed dispersal
	//the idea is that this structure contains all of the cells that the cell represented can possibly disperse seeds to and the probability of the happening for each cell
	int size, seeds_present, seeds_received, *cells; //seeds_present && seeds_received are treated as boolean values... cells & prob are to be the length of size
	float *prob, lyppt;
}typedef Grid_SD_St;

struct _grid_sxw_st
{ //holds pointers dynamically allocated by SXW->c
	RealD *roots_max, *rootsXphen, *roots_active, *roots_active_rel,
			*roots_active_sum, *phen, *prod_bmass, *prod_pctlive;
}typedef Grid_SXW_St;

struct _grid_init_species_st
{
	/* TRUE is this cell should use spinup */
	int use_SpinUp;
	/* Array of Boolean values. TRUE if given species
	   should be included in spinup */
	int *species_seed_avail;
}typedef Grid_Init_Species_St;

/* Struct to hold all plot-specific parameters */
struct grid_cell_st
{
	/* RGroup coresponding to this cell */
	GroupType **myGroup;
	/* Species corresponding to this cell */
	SpeciesType **mySpecies;
	/* Succulents corresponding to this cell */
	SucculentType mySucculent;
	/* This cell's environment. We expect each cell to
	 * have slightly different weather each year */
	EnvType myEnvironment;
	/* Cell's plot data */
	PlotType myPlot;
	/* Global variables corresponding to this cell */ 
	ModelType myGlobals;
	/* If TRUE this cell should use seed dispersal */
	Bool useSeedDispersal;
	/* TRUE if this cell is in spinup mode */
	Bool duringSpinup;
	/* species spinup information */
	Grid_Init_Species_St mySpeciesInit;
	/* seed dispersal information corresponding to this cell */
	Grid_SD_St *mySeedDispersal;

	Bool* someKillage;
	
	/* ---------------- accumulators -------------------- */
	StatType *_Dist, *_Ppt, *_Temp,
  		*_Grp, *_Gsize, *_Gpr, *_Gmort, *_Gestab,
  		*_Spp, *_Indv, *_Smort, *_Sestab, *_Sreceived;
	FireStatsType *_Gwf;
	Bool stats_init;
	/* -------------- end accumulators ------------------ */

	/* -------------------- SXW ------------------------- */
	transp_t* myTranspWindow;
	SXW_t* mySXW;
	SXW_resourceType* mySXWResources;
	/* ------------------ End SXW ----------------------- */

	/* ------------------- Soils ------------------------ */
	// Soil layer information for this cell.
	SoilType mySoils;
	/* ------------------ End Soils --------------------- */
} typedef CellType;

/************ Module Variable Declarations ***************/
/***********************************************************/

/* Indices for grid_directories go here */
enum
{
    GRID_DIRECTORY_STEPWAT_INPUTS,

    /* Automatically generate number of directories since enums start at 0 */
    N_GRID_DIRECTORIES
};

/* Indices for grid_files go here */
enum
{
    GRID_FILE_LOGFILE,
    GRID_FILE_SETUP,
    GRID_FILE_DISTURBANCES,
    GRID_FILE_SOILS,
    GRID_FILE_SEED_DISPERSAL,
    GRID_FILE_INIT_SPECIES,
    GRID_FILE_FILES,
    GRID_FILE_MAXRGROUPSPECIES,

    GRID_FILE_PREFIX_BMASSAVG,
    GRID_FILE_PREFIX_MORTAVG,
    GRID_FILE_PREFIX_RECEIVEDPROB,
    GRID_FILE_PREFIX_BMASSCELLAVG,

    /* Automatically generate number of files since enums start at 0 */
    N_GRID_FILES
};

char *grid_files[N_GRID_FILES], *grid_directories[N_GRID_DIRECTORIES], sd_Sep;

int grid_Cols, grid_Rows, grid_Cells, sd_NYearsSeedsAvailable;
int UseDisturbances, UseSoils, sd_DoOutput, sd_MakeHeader, sd_Option1a,
		sd_Option1b, sd_Option2a, sd_Option2b; //these are treated like booleans

// these variables are for storing the globals in STEPPE... they are dynamically allocated/freed
SpeciesType **grid_Species, **spinup_Species;
GroupType **grid_RGroup, **spinup_RGroup;
SucculentType *grid_Succulent, *spinup_Succulent;
EnvType *grid_Env, *spinup_Env;
PlotType *grid_Plot, *spinup_Plot;
ModelType *grid_Globals, *spinup_Globals;

CellType** gridCells;

// these two variables are for storing SXW variables... also dynamically allocated/freed
SXW_t *grid_SXW, *spinup_SXW;
Grid_SXW_St *grid_SXW_ptrs, *spinup_SXW_ptrs;

// these are SOILWAT variables that we need...
extern SW_SOILWAT SW_Soilwat;
extern SW_SITE SW_Site;
extern SW_VEGPROD SW_VegProd;
extern SW_WEATHER SW_Weather;
extern pcg32_random_t grid_rng; //this file's unique random number generator

/* We need to seed these RNGs when using the gridded mode but do not use them in this file. */
extern pcg32_random_t environs_rng;
extern pcg32_random_t mortality_rng;
extern pcg32_random_t resgroups_rng;
extern pcg32_random_t species_rng;
extern pcg32_random_t markov_rng;

//This is Rgroup structure pointer that will read rgroup disturbance value, will be used in
// grid disturbance
//extern GroupType     *RGroup [MAX_RGROUPS]; // don't need to extern here, because this is done in ST_Globals->h

// these are grids to store the SOILWAT variables... also dynamically allocated/freed
SW_SOILWAT *grid_SW_Soilwat, *spinup_SW_Soilwat;
SW_SITE *grid_SW_Site, *spinup_SW_Site;
SW_VEGPROD *grid_SW_VegProd, *spinup_SW_VegProd;
SW_MODEL *grid_SW_Model, *spinup_SW_Model;

// these two variables are used to store the soil/distubance inputs for each grid cell... also dynamically allocated/freed
SoilType *grid_Soils;
Grid_Disturb_St *grid_Disturb;
Grid_Init_Species_St *grid_initSpecies;

Grid_SD_St **grid_SD; //for seed dispersal

// these variables are used for the soil types in the spinup options
int nSoilTypes, *soilTypes_Array, *grid_SoilTypes;

extern Bool UseProgressBar;
extern Bool* _SomeKillage;

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/

// Declare functions defined elsewhere:
void runGrid(void); //to be called from ST_main.c
void _kill_annuals(void);
void _kill_extra_growth(void);
void rgroup_Extirpate(GrpIndex rg);

/************* External Function Declarations **************/
/***********************************************************/

//from ST_species.c
void proportion_Recovery(void);
void save_annual_species_relsize(void);


//from ST_resgroups.c
void rgroup_Grow(void);
void rgroup_Establish(void);
void rgroup_IncrAges(void);
void rgroup_PartResources(void);
//void rgroup_ResPartIndiv(void);

//from ST_mortality.c
void mort_Main(Bool *killed);
void mort_EndOfYear(void);
void grazing_EndOfYear(void);

//functions from ST_params.c
void parm_Initialize(void);
void parm_SetFirstName(char *s);
void parm_SetName(char *s, int which);
void parm_free_memory(void);

//from ST_main.c
void Plot_Initialize(void);

//functions from ST_stats.c
void stat_Collect(Int year);
void stat_Collect_GMort(void);
void stat_Collect_SMort(void);
void stat_Output_AllMorts(void);
void stat_Output_AllBmass(void);
//Adding functions for creating grid cells avg values output file
void stat_Output_AllBmassAvg(void);
void stat_Output_AllCellAvgBmass(const char * filename);
void stat_Output_Seed_Dispersal(const char * filename, const char sep,
		Bool makeHeader);
void stat_Load_Accumulators(int cell, int year);
void stat_Save_Accumulators(int cell, int year);
void stat_Free_Accumulators(void);
void stat_Init_Accumulators(void);
void stat_Copy_Accumulators(StatType* newDist, StatType* newPpt, StatType* newTemp, StatType* newGrp, StatType* newGsize, 
                            StatType* newGpr, StatType* newGmort, StatType* newGestab, StatType* newSpp, StatType* newIndv,
                            StatType* newSmort, StatType* newSestab, StatType* newSrecieved, FireStatsType* newGwf, Bool firstTime);

//functions from sxw.c
//void free_sxw_memory( void );
//void free_all_sxw_memory(void);
void load_sxw_memory(RealD * grid_roots_max, RealD* grid_rootsXphen,
		RealD* grid_roots_active, RealD* grid_roots_active_rel,
		RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass,
		RealD* grid_prod_pctlive);
void save_sxw_memory(RealD * grid_roots_max, RealD* grid_rootsXphen,
		RealD* grid_roots_active, RealD* grid_roots_active_rel,
		RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass,
		RealD* grid_prod_pctlive);
void _deallocate_memory(void);
//void SXW_init( Bool init_SW );

void copy_sxw_variables(SXW_t* newSXW, SXW_resourceType* newSXWResources, transp_t* newTransp_window);

void maxrgroupspecies_init(void);
void files_init(void);

SXW_t* getSXW(void);
SXW_resourceType* getSXWResources(void);
transp_t* getTranspWindow(void);

/*********** Locally Used Function Declarations ************/
/***********************************************************/

static int _load_bar(char* prefix, clock_t start, int x, int n, int r, int w);
static double _time_remaining(clock_t start, char* timeChar, double percentDone);
static void _run_spinup(void);
static void _init_grid_files(void);
static void _init_grid_inputs(void);
static void _init_SXW_inputs(Bool init_SW, char *f_roots);
static void _init_stepwat_inputs(void);
static void _init_grid_globals(void);
static void _init_spinup_globals(void);
static void allocate_gridCells(int rows, int cols);
static void allocate_accumulators(void);
static void _load_grid_globals(void);
static void _load_spinup_globals(void);
static void _free_grid_memory(void);
static void _free_spinup_memory(void);
static void _free_grid_globals(void);
static void _free_spinup_globals(void);
static void _load_cell(int row, int col, int year, Bool useAccumulators);
static void load_cell(int row, int col);
static void unload_cell(void);
static void _load_spinup_cell(int cell);
static void _save_cell(int row, int col, int year, Bool useAccumulators);
static void _save_spinup_cell(int cell);
static void _read_disturbances_in(void);
static void _read_soils_in(void);
static void _init_soil_layers(int cell, int isSpinup);
static float _read_a_float(FILE *f, char *buf, const char *filename,
		const char *descriptor);
static float _cell_dist(int row1, int row2, int col1, int col2, float cellLen);
static void _read_seed_dispersal_in(void);
static void _do_seed_dispersal(void);
static void _set_sd_lyppt(int row, int col);
static void _kill_groups_and_species(void);
static int _do_grid_disturbances(int row, int col);
static void _read_init_species(void);
static void _do_groups_and_species_extirpate(void);
static void _do_grid_proportion_Recovery(int row, int col);
static void _do_grid_grazing_EndOfYear(int row, int col);

static IndivType* _create_empty_indv(void); //these 5 functions are used for copying/freeing the linked list of individuals correctly...
static void _free_individuals(IndivType *head);
static void _free_head(IndivType *head);
static IndivType* _copy_individuals(IndivType *head);
static IndivType* _copy_head(IndivType *head);

static void _read_maxrgroupspecies(void);
static void _read_grid_setup(void);
static void _read_files(void);

/******************** Begin Model Code *********************/
/***********************************************************/

/***********************************************************/
double _time_remaining(clock_t start, char* timeChar, double percentDone)
{
	// for use in _load_bar() function... pretty much it returns the amount of time left and a character in timeChar representing the units
	// percent done must be from 0 to 1, with 0 representing no work done (0%) and 1 representing all of the work done (ie 100%)
	clock_t timeElapsed = clock() - start; //gets the time since we started (this is in CPU cycles).  We convert to seconds in the next part
	double result = (((double) timeElapsed) / CLOCKS_PER_SEC) / percentDone
			* (1.0 - percentDone); //gets estimated time left until the work is complete based upon timeElapsed and percentDone

	*timeChar = 's'; //get the biggest time units that we can use (ie no use displaying 120 seconds if we can display 2 minutes)
	if (result > 60)
	{
		result /= 60.0;
		*timeChar = 'm';
		if (result > 60)
		{
			result /= 60.0;
			*timeChar = 'h';
			if (result > 24)
			{
				result /= 24.0;
				*timeChar = 'd';
				if (result > 7)
				{
					result /= 7.0;
					*timeChar = 'w';
				}
			}
		}
	}

	return result;
}

/***********************************************************/
static int _load_bar(char* prefix, clock_t start, int x, int n, int r, int w)
{
	//loads a progress bar, x is how much progress has been made, n is how much progress to be done, r is how many times to update, and w is width of the progress bar
	//wouldn't suggest using this if your stdout is being redirected somewhere else (it doesn't work right in that case and gives unwanted output)
	//gotten from: http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/

	//modified to include an estimate of the time remaining... start must be a clock started (using clock() function) right before beginning the work

	// Only update r times.
	if (x % (n / r) != 0)
		return (int) ((x / (float) n) * 100);

	int i;

	printf("\r"); //we output a carriage-return character to put us back at the beginning of the line
	if (prefix != NULL)
		printf("%s", prefix);

	// Calculuate the ratio of complete-to-incomplete.
	float ratio = x / (float) n;
	int c = ratio * w;
	int result = (int) (ratio * 100);

	if (result > 1)
	{ //we don't give an estimate if less than 2% of the work is complete because we don't have any data to go off of...
		char timeChar;
		double timeLeft = _time_remaining(start, &timeChar, (double) ratio);
		if (timeLeft < 10)
			printf("(est 0%.2f%c) ", timeLeft, timeChar);
		else
			printf("(est %.2f%c) ", timeLeft, timeChar);
	}

	// Show the percentage complete.
	printf("%3d%% [", (int) (ratio * 100));

	// Show the load bar.
	for (i = 0; i < c; i++)
		printf("=");

	for (i = c; i < w; i++)
		printf(" ");

	printf("]");
	fflush(stdout); //we do this to flush (ie. print) the output, since stdout typically waits for a newline character before flushing
	return result;
}

/* Runs the gridded version of the code */
void runGrid(void)
{
	int i, j;
	Bool killedany;
	IntS year, iter;
	double prog_Percent = 0.0, prog_Incr, prog_Acc = 0.0;
	char prog_Prefix[32];
	clock_t prog_Time;

	_init_grid_files();				// reads in files.in file
	_read_maxrgroupspecies();       // reads in maxrgroupspecies.in file
	_read_grid_setup();             // reads in grid_setup.in file
    _read_files();                  // reads in Stepwat_Inputs/files.in file
    _init_stepwat_inputs();			// reads the stepwat inputs in
	_init_grid_inputs();			// reads the grid inputs in & initializes the global grid variables
	//SWC hist file prefix needs to be cleared
	Mem_Free(SW_Soilwat.hist.file_prefix);
	SW_Soilwat.hist.file_prefix = NULL;

	/* Print some general information to stdout */
	load_cell(0,0);
	printf("Number of layers: %d\n", SW_Site.n_layers);
    printf("Number of iterations: %d\n", Globals->runModelIterations);
    printf("Number of years: %d\n", Globals->runModelYears);
	printf("Number of cells: %d\n\n", grid_Cells);
	if(UseDisturbances) printf("Using grid Disturbances file.\n");
	if(UseSoils) printf("Using grid soils file\n");
	if(UseSeedDispersal){
		 printf("Seeds availible for %d years at the start of the simulation.\n",sd_NYearsSeedsAvailable);
		 if(sd_Option1a) printf("All species and cells are eligible for seed dispersal.\n");
	}
	if(sd_Option2a) printf("Running Spinup for %s cells.\n", sd_Option2b ? "all" : "specified");
	unload_cell();
	/* END printing general info */

	if (sd_Option2a || sd_Option2b)
		_run_spinup();				// does the initial spinup

	// SOILWAT resets SW_Weather.name_prefix every iteration. This is not the behavior we want 
	// so the name is stored here.
	char SW_prefix_permanent[2048];
	sprintf(SW_prefix_permanent, "%s/%s", grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS], SW_Weather.name_prefix);

	for (iter = 1; iter <= gridCells[0][0].myGlobals.runModelIterations; iter++)
	{ //for each iteration

		/*
		 * 06/15/2016 (akt) Added resetting correct historical weather file path,
		 * as it was resetting to original path value (that was not correct for grid version)from input file after every iteration
		 */
		sprintf(SW_Weather.name_prefix, "%s", SW_prefix_permanent); //updates the directory correctly for the weather files so soilwat can find them

		if (BmassFlags.yearly || MortFlags.yearly)
			parm_Initialize();

		// Initialize the plot for each grid cell
		for (i = 0; i < grid_Rows; i++){
			for (j = 0; j < grid_Cols; j++){
				load_cell(i, j);
				Plot_Initialize();
				Globals->currIter = iter;
			}
		}
		unload_cell(); // Reset the global variables

		RandSeed(gridCells[0][0].myGlobals.randseed, &environs_rng);
		RandSeed(gridCells[0][0].myGlobals.randseed, &mortality_rng);
		RandSeed(gridCells[0][0].myGlobals.randseed, &resgroups_rng);
		RandSeed(gridCells[0][0].myGlobals.randseed, &species_rng);
		RandSeed(gridCells[0][0].myGlobals.randseed, &grid_rng);
		RandSeed(gridCells[0][0].myGlobals.randseed, &markov_rng);

		for (year = 1; year <= gridCells[0][0].myGlobals.runModelYears; year++)
		{ //for each year
			for (i = 0; i < grid_Rows; i++){
				for (j = 0; j < grid_Cols; j++)
				{ //for each cell

					//fprintf(stderr, "year: %d", year);
					load_cell(i, j);

					Globals->currYear = year;

					if (year > 1 && UseSeedDispersal)
						_set_sd_lyppt(i, j);

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
					proportion_Recovery(); 		// Recover from any disturbances
					_kill_extra_growth(); 		// Kill superfluous growth

				} /* end model run for this cell*/
			} /* end model run for this row */
			if (UseSeedDispersal)
				_do_seed_dispersal();
			
			unload_cell(); // Reset the global variables
		}/* end model run for this year*/

		// collects the data appropriately for the mort output... (ie. fills the accumulators in ST_stats.c with the values that they need)
		if (MortFlags.summary)
			for (i = 1; i <= grid_Rows; i++)
				for (j = 1; j <= grid_Cols; j++)
				{
					_load_cell(i, j, Globals->runModelYears, TRUE);
					stat_Collect_GMort();
					stat_Collect_SMort();
					_save_cell(i, j, Globals->runModelYears, TRUE);
				}
		//reset soilwat to initial condition
		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
		SXW_Reset();
		Mem_Free(SW_Soilwat.hist.file_prefix);
		SW_Soilwat.hist.file_prefix = NULL;
		ChDir("..");

	} /*end iterations */

	// outputs all of the mort and BMass files for each cell...
	for (i = 0; i < grid_Rows; i++){
		for (j = 0; j < grid_Cols; j++)
		{
			int cell = j + ((i - 1) * grid_Cols) - 1;
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
				stat_Output_AllBmassAvg();
			if (UseSeedDispersal && sd_DoOutput)
				stat_Output_Seed_Dispersal(fileReceivedProb, sd_Sep,
						sd_MakeHeader);
		}
	}
	unload_cell(); // Reset the global variables

	//Here creating grid cells avg values output file
	char fileBMassCellAvg[1024];
	sprintf(fileBMassCellAvg, "%s.csv", grid_files[GRID_FILE_PREFIX_BMASSCELLAVG]);
	if (BmassFlags.summary)
		stat_Output_AllCellAvgBmass(fileBMassCellAvg);

	_free_grid_memory(); // free our allocated memory since we do not need it anymore
	_deallocate_memory(); // sxw memory.
}

/* "Spinup" the model by running for Globals->runModelYears without seed dispersal or statistics outputs. */
static void _run_spinup(void)
{
	/* Dummy accumulators to ensure we do not collect statistics */
	StatType *dummy_Dist, *dummy_Ppt, *dummy_Temp,
  		*dummy_Grp, *dummy_Gsize, *dummy_Gpr, *dummy_Gmort, *dummy_Gestab,
  		*dummy_Spp, *dummy_Indv, *dummy_Smort, *dummy_Sestab, *dummy_Sreceived;
	FireStatsType *dummy_Gwf;

	/* ints used for iterating over gridCells */
	int i, j;

	/* ints used for iterating over years and iterations */
	IntS year, iter;

	/* For iterating over RGroup and Species */
	GrpIndex rg;
	SppIndex sp;
	IntS k;

	/* 'Globals' is not allocated in gridded mode. Therefore we need to implicitely remember years
	   and iterations. gridCells[0][0] is guaranteed to exist because it is populated first. */
	int iterations = gridCells[0][0].myGlobals.runModelIterations;
	int total_years = gridCells[0][0].myGlobals.runModelYears;

	/* killedany for mortality functions, temporary_storage for swapping variables */
	Bool killedany, temporary_storage;

	DuringSpinup = TRUE;

	if (!UseSoils)
	{ // if we're not using inputting soils then there is simply one soil type as all the soils are the same
		nSoilTypes = 1;
		soilTypes_Array[0] = 0;
		for (i = 0; i < grid_Cells; i++)
			grid_SoilTypes[i] = 0;
	}

	for (iter = 1; iter <= 1; iter++)
	{ //for each iteration... only 1 iteration allowed for now

		/* Since this is technically an iteration so we need to seed the RNGs. */
		RandSeed(Globals->randseed, &environs_rng);
		RandSeed(Globals->randseed, &mortality_rng);
		RandSeed(Globals->randseed, &resgroups_rng);
		RandSeed(Globals->randseed, &species_rng);
		RandSeed(Globals->randseed, &grid_rng);
		RandSeed(Globals->randseed, &markov_rng);

		if (BmassFlags.yearly || MortFlags.yearly)
			parm_Initialize();

		Plot_Initialize();

		/* Before we start iterating we need to swap Species[sp]->use_me and mySpeciesInit.species_seed_avail[sp].
		   species_seed_avail is an array of booleans that represent whether the given species should be used 
		   in spinup. use_me is a boolean that represents whether the given species should be used in production.
		   By swaping them we save space, but we have to remember to swap them back before the production run. */
		for(i = 0; i < grid_Rows; ++i){
			for(j = 0; j < grid_Cols; ++j){
				load_cell(i, j);	/* We could do this without loading the cell, but there would be no guarantee
									   that ForEachGroup would iterate correctly */

				Globals->currIter = iter; // We need to do this before the "years" loop.

				/* Begin swaping variables */
				ForEachGroup(rg){
					ForEachGroupSpp(sp, rg, k){
						// Temporarily store use_me
						temporary_storage = Species[sp]->use_me;
						// Swap use_me
						Species[sp]->use_me = gridCells[i][j].mySpeciesInit.species_seed_avail[sp]; 
						// Swap species_seed_avail[sp]
						gridCells[i][j].mySpeciesInit.species_seed_avail[sp] = temporary_storage;
					} /* End for each species */
				} /* End for each group */
			} /* End for each column */
		} /* End for each row */
		unload_cell(); // Reset the global variables

		for (year = 1; year <= total_years; year++)
		{ //for each year
			for (i = 0; i < grid_Rows; ++i)
			{ // for each row
				for(j = 0; j < grid_Cols; ++j)
				{ // for each column
					// If we should run spinup on this cell
					if(gridCells[i][j].mySpeciesInit.use_SpinUp){
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

					rgroup_Establish(); 		// Establish individuals. Excludes annuals.

					Env_Generate();				// Generated the SOILWAT environment

					rgroup_PartResources();		// Distribute resources
					rgroup_Grow(); 				// Grow

					mort_Main(&killedany); 		// Mortality that occurs during the growing season

					rgroup_IncrAges(); 			// Increment ages of all plants

					grazing_EndOfYear(); 		// Livestock grazing
					
					//save_annual_species_relsize(); // See stat_Collect below for why this is commented out.

					mort_EndOfYear(); 			// End of year mortality.

					/* stat_Collect was called in spinup, but I have no idea why. We definitely
					   do not want to use the accumulator values, so why collect them? -Chandler */
					//stat_Collect(year);

				    _kill_annuals(); 			// Kill annuals
					proportion_Recovery(); 		// Recover from any disturbances
					_kill_extra_growth(); 		// Kill superfluous growth			
				} /* end column */
			} /* end row */

			unload_cell(); // Reset the global variables
		} /* end model run for this year*/

		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
		SXW_Reset();
		//TODO: This is a shortcut. swc history is not used and shouldn't be until this is fixed.
		Mem_Free(SW_Soilwat.hist.file_prefix);
		SW_Soilwat.hist.file_prefix = NULL;
		ChDir("..");
		//_free_grid_globals(); //free's the grid variables that change every iter
	} /* End iterations */

	/* Swap back Species[sp]->use_me and species_seed_avail[sp]. */
	for(i = 0; i < grid_Rows; ++i){
		for(j = 0; j < grid_Cols; ++j){
			load_cell(i, j);	/* We could do this without loading the cell, but there would be no guarantee
								   that ForEachGroup would iterate correctly */

			/* Begin swaping variables */
			ForEachGroup(rg){
				ForEachGroupSpp(sp, rg, k){
					// Temporarily store use_me
					temporary_storage = Species[sp]->use_me;
					// Swap use_me
					Species[sp]->use_me = gridCells[i][j].mySpeciesInit.species_seed_avail[sp]; 
					// Swap species_seed_avail[sp]
					gridCells[i][j].mySpeciesInit.species_seed_avail[sp] = temporary_storage;
				} /* End for each species */
			} /* End for each group */
		} /* End for each column */
	} /* End for each row */
	unload_cell(); // Reset the global variables

	DuringSpinup = FALSE;
}

/***********************************************************/
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

	/*printf("stepwat dir: %s\n", grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
	 for(i = 0; i < N_GRID_FILES; i++)
	 printf("%d : %s\n", i, grid_files[i]);*/

	CloseFile(&f);
}

/***********************************************************/
static void _init_grid_inputs(void)
{
    int i, j;

    _init_grid_globals();

	if (UseDisturbances)
		_read_disturbances_in();
	if (UseSeedDispersal)
	{
		_read_seed_dispersal_in();
		_read_init_species();
	}
	if (UseSoils)
		_read_soils_in();

    for(i = 0; i < grid_Rows; ++i) {
        for(j = 0; j < grid_Cols; ++j) {
            load_cell(i, j);
            gridCells[i][j].duringSpinup = FALSE;
        }
    }

    unload_cell();
}
/***********************************************************/
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
   DEPENDENCYS: allocate_gridCells() must be called first. */
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
			
			// species_seed_avail is a dynamically allocated array
			gridCells[i][j].mySpeciesInit.species_seed_avail = (int*)
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
	for(int i = 0; i < grid_Rows; ++i){
		for(int j = 0; j < grid_Cols; ++j){
			/* load_cell is not necessary for the actual accumulators, but it is necessary for
			   the ForEach loops. We still have to refer to the accumulators as
			   gridCells[i][j].<accumulator> because the ST_stats accumulators are local. */
			load_cell(i,j);

  			if (BmassFlags.dist) {
    			gridCells[i][j]._Dist = (StatType*) Mem_Calloc(1, sizeof(StatType), "allocate_accumulators(Dist)");
		    	gridCells[i][j]._Dist->s = (struct accumulators_st *)
		               		Mem_Calloc( Globals->runModelYears,
		                           		sizeof(struct accumulators_st),
 		                         		"allocate_accumulators(Dist)");
		  	}
 		 	if (BmassFlags.ppt) {
		    	gridCells[i][j]._Ppt = (StatType*) Mem_Calloc(1, sizeof(StatType), "allocate_accumulators(PPT");
		    	gridCells[i][j]._Ppt->s  = (struct accumulators_st *)
		               		Mem_Calloc( Globals->runModelYears,
		                           		sizeof(struct accumulators_st),
		                          		"allocate_accumulators(PPT)");
 		 	}
		  	if (BmassFlags.tmp) {
		    	gridCells[i][j]._Temp = (StatType*) Mem_Calloc(1, sizeof(StatType), "allocate_accumulators(Temp)");
		    	gridCells[i][j]._Temp->s = (struct accumulators_st *)
		               		Mem_Calloc( Globals->runModelYears,
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
             			Mem_Calloc( Globals->runModelYears,
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
             							Mem_Calloc( Globals->runModelYears,
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
             							Mem_Calloc( Globals->runModelYears,
                    		     					sizeof(struct accumulators_st),
                        							"allocate_accumulators(Gpr[rg].s)");
					}
    			}

    			if (BmassFlags.wildfire || BmassFlags.prescribedfire) {
      				gridCells[i][j]._Gwf = (struct fire_st *)
             				Mem_Calloc( 1, sizeof(struct fire_st),
                		        		"allocate_accumulators(Gwf)");

      				gridCells[i][j]._Gwf->wildfire = (int *) Mem_Calloc( 1,
                    		  		 sizeof(int) * Globals->runModelYears,
                    		  		 "allocate_accumulators(Gwf->wildfire)");
      
      				gridCells[i][j]._Gwf->prescribedFire = (int **) Mem_Calloc( 1,
                      						sizeof(int **) * SuperGlobals.max_rgroups,
                       						"allocate_accumulators(Gwf->prescribedfire");

      				ForEachGroup(rg){
        				gridCells[i][j]._Gwf->prescribedFire[rg] = (int *)
          										   Mem_Calloc( Globals->runModelYears,
                      										   sizeof(int) * Globals->runModelYears,
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
               			 		  Mem_Calloc( Globals->runModelYears,
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
               				  		Mem_Calloc( Globals->runModelYears,
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
					  									Mem_Calloc( Globals->runModelYears,
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

/***********************************************************/
static IndivType* _create_empty_indv(void)
{
	//simply allocates memory for an individual and returns the pointer to it

	return Mem_Calloc(1, sizeof(IndivType), "_create_empty_indv()");

}

/***********************************************************/
static void _free_individuals(IndivType *head)
{
	//frees the memory allocated in the linked list of individuals pointed to by head
	if (head->Next != NULL)
		_free_individuals(head->Next); //recursively call itself in order to free the next obj in the list (Mem_Free() will end up getting called on the objects in the list from last to first)
	Mem_Free(head);
}

/***********************************************************/
static void _free_head(IndivType * head)
{
	if (head == NULL)
		return; //if head == NULL then there is no memory to free, so just return
	_free_individuals(head);
}

/***********************************************************/
static IndivType* _copy_individuals(IndivType *head)
{
	//performs a deep copy (not to be confused with a shallow copy) since that's what we want
	IndivType *n = _create_empty_indv(); //allocate space for the individual to be copied to

	*n = *head;

	if (head->Next != NULL)
	{
		n->Next = _copy_individuals(head->Next); //recursively call itself in order to copy the next obj in the list
		n->Next->Prev = n;
	}
	else
		n->Next = NULL;

	return n;
}

/***********************************************************/
static IndivType* _copy_head(IndivType *head)
{
	if (head == NULL)
		return NULL; //if head == NULL then there is nothing to copy...
	return _copy_individuals(head);
}

/***********************************************************/
static void _init_grid_globals(void)
{
	//initializes grid variables, allocating the memory necessary for them (this step is only needed to be done once)

    // Keep these until I know if the CellType struct accounts for these variables.
	soilTypes_Array = Mem_Calloc(grid_Cells, sizeof(int*),
			"_init_grid_globals()");
	grid_SoilTypes = Mem_Calloc(grid_Cells, sizeof(int*),
			"_init_grid_globals()");
}

/***********************************************************/
static void _init_spinup_globals(void)
{
	//initializes spinup variables, allocating the memory necessary for them (this step is only needed to be done once)

	GrpIndex c;
	SppIndex s;

	spinup_Succulent = Mem_Calloc(nSoilTypes, sizeof(SucculentType),
			"_init_spinup_globals()");
	spinup_Env = Mem_Calloc(nSoilTypes, sizeof(EnvType),
			"_init_spinup_globals()");
	spinup_Plot = Mem_Calloc(nSoilTypes, sizeof(PlotType),
			"_init_spinup_globals()");
	spinup_Globals = Mem_Calloc(nSoilTypes, sizeof(ModelType),
			"_init_spinup_globals()");

	ForEachSpecies(s)
		if (Species[s]->use_me)
			spinup_Species[s] = Mem_Calloc(nSoilTypes, sizeof(SpeciesType),
					"_init_spinup_globals()");
	ForEachGroup(c)
		if (RGroup[c]->use_me)
			spinup_RGroup[c] = Mem_Calloc(nSoilTypes, sizeof(GroupType),
					"_init_spinup_globals()");

	spinup_SXW = Mem_Calloc(nSoilTypes, sizeof(SXW_t),
			"_init_spinup_globals()");
	spinup_SW_Soilwat = Mem_Calloc(nSoilTypes, sizeof(SW_SOILWAT),
			"_init_spinup_globals()");
	spinup_SW_Site = Mem_Calloc(nSoilTypes, sizeof(SW_SITE),
			"_init_spinup_globals()");
	spinup_SW_VegProd = Mem_Calloc(nSoilTypes, sizeof(SW_VEGPROD),
			"_init_spinup_globals()");
	if (UseSoils)
	{
		spinup_SXW_ptrs = Mem_Calloc(nSoilTypes, sizeof(Grid_SXW_St),
				"_init_spinup_globals()");
	}
}

/***********************************************************/
static void _load_grid_globals(void)
{
	//this initializes/allocates memory needed... this step is needed to be done for every iteration

	int i, j, k;
	GrpIndex c;
	SppIndex s;

	if (UseSoils)
		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]); //change the directory for _init_soil_layers()
	for (i = 0; i < grid_Cells; i++)
	{

		ForEachSpecies(s)
		{ //macros defined in ST_defines.h
			if (!Species[s]->use_me)
				continue;
			grid_Species[s][i] = *Species[s];

			grid_Species[s][i].kills = Mem_Calloc(Species[s]->max_age,
					sizeof(IntUS), "_init_grid_globals()");
			grid_Species[s][i].seedprod = Mem_Calloc(Species[s]->viable_yrs,
					sizeof(RealF), "_init_grid_globals()");

			memcpy(grid_Species[s][i].kills, Species[s]->kills,
					Species[s]->max_age * sizeof(IntUS));
			memcpy(grid_Species[s][i].seedprod, Species[s]->seedprod,
					Species[s]->viable_yrs * sizeof(RealF));

			grid_Species[s][i].IndvHead = _copy_head(Species[s]->IndvHead); //copy_head() deep copies the structure (allocating memory when needed)... it will even allocate memory for the head of the list
		}

		ForEachGroup(c)
		{
			if (!RGroup[c]->use_me)
				continue;
			grid_RGroup[c][i] = *RGroup[c];
			grid_RGroup[c][i].kills = Mem_Calloc(RGroup[c]->max_age,
					sizeof(IntUS), "_init_grid_globals()");

			memcpy(grid_RGroup[c][i].kills, RGroup[c]->kills,
					RGroup[c]->max_age * sizeof(IntUS));
			if (UseDisturbances)
			{
				grid_RGroup[c][i].killyr = grid_Disturb[i].kill_yr;
				grid_RGroup[c][i].killfreq = grid_Disturb[i].killfrq;
				grid_RGroup[c][i].extirp = grid_Disturb[i].extirp;
			}
		}

		grid_Succulent[i] = *Succulent;
		grid_Env[i] = *Env;
		grid_Plot[i] = *Plot;
		grid_Globals[i] = *Globals;

		if (UseDisturbances)
		{
			grid_Globals[i].pat.use = grid_Disturb[i].choices[0];
			grid_Globals[i].mound.use = grid_Disturb[i].choices[1];
			grid_Globals[i].burrow.use = grid_Disturb[i].choices[2];
		}
		if (UseSoils)
			_init_soil_layers(i, 0);

		grid_SXW[i] = *SXW;

		grid_SXW[i].f_roots = Str_Dup(SXW->f_roots);
		grid_SXW[i].f_phen = Str_Dup(SXW->f_phen);
		grid_SXW[i].f_bvt = Str_Dup(SXW->f_bvt);
		grid_SXW[i].f_prod = Str_Dup(SXW->f_prod);
		grid_SXW[i].f_watin = Str_Dup(SXW->f_watin);

		grid_SXW[i].transpTotal = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
				sizeof(RealD), "_init_grid_globals()");
		ForEachVegType(k) {
			grid_SXW[i].transpVeg[k] = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
					sizeof(RealD), "_init_grid_globals()");
		}
		grid_SXW[i].swc = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealF),
				"_init_grid_globals()");

		memcpy(grid_SXW[i].transpTotal, SXW->transpTotal,
				SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
		ForEachVegType(k) {
			memcpy(grid_SXW[i].transpVeg[k], SXW->transpVeg[k],
					SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
		}
		memcpy(grid_SXW[i].swc, SXW->swc,
				SXW->NPds * SXW->NSoLyrs * sizeof(RealF));

		grid_SW_Soilwat[i] = SW_Soilwat;
		grid_SW_Site[i] = SW_Site;
		grid_SW_VegProd[i] = SW_VegProd;

		grid_SW_Site[i].lyr = Mem_Calloc(
				SW_Site.n_layers + SW_Site.deepdrain,
				sizeof(SW_LAYER_INFO *), "_init_grid_globals()");
		for (j = 0; j < SW_Site.n_layers + SW_Site.deepdrain; j++)
		{
			grid_SW_Site[i].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
					"_init_grid_globals()");
			memcpy(grid_SW_Site[i].lyr[j], SW_Site.lyr[j],
					sizeof(SW_LAYER_INFO));
		}
	}
	if (UseSoils)
		ChDir(".."); //get back to our previous directory

}

/***********************************************************/
static void _load_spinup_globals(void)
{
	//this initializes/allocates memory needed... this step is needed to be done for every iteration

	int i, j, k;
	GrpIndex c;
	SppIndex s;

	if (UseSoils)
		ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]); //change the directory for _init_soil_layers()
	for (i = 0; i < nSoilTypes; i++)
	{

		int cell = soilTypes_Array[i]; //this is the cell number of the first cell representing this soil type

		ForEachSpecies(s)
		{ //macros defined in ST_defines.h
			if (!Species[s]->use_me)
				continue;
			spinup_Species[s][i] = *Species[s];

			spinup_Species[s][i].kills = Mem_Calloc(Species[s]->max_age,
					sizeof(IntUS), "_init_spinup_globals()");
			spinup_Species[s][i].seedprod = Mem_Calloc(Species[s]->viable_yrs,
					sizeof(RealF), "_init_spinup_globals()");

			memcpy(spinup_Species[s][i].kills, Species[s]->kills,
					Species[s]->max_age * sizeof(IntUS));
			memcpy(spinup_Species[s][i].seedprod, Species[s]->seedprod,
					Species[s]->viable_yrs * sizeof(RealF));

			spinup_Species[s][i].IndvHead = _copy_head(Species[s]->IndvHead); //copy_head() deep copies the structure (allocating memory when needed)... it will even allocate memory for the head of the list
		}

		ForEachGroup(c)
		{
			if (!RGroup[c]->use_me)
				continue;
			spinup_RGroup[c][i] = *RGroup[c];
			spinup_RGroup[c][i].kills = Mem_Calloc(RGroup[c]->max_age,
					sizeof(IntUS), "_init_spinup_globals()");

			memcpy(spinup_RGroup[c][i].kills, RGroup[c]->kills,
					RGroup[c]->max_age * sizeof(IntUS));
			if (UseDisturbances)
			{
				spinup_RGroup[c][i].killyr = grid_Disturb[cell].kill_yr;
				spinup_RGroup[c][i].killfreq = grid_Disturb[cell].killfrq;
				spinup_RGroup[c][i].extirp = grid_Disturb[cell].extirp;
			}
		}

		spinup_Succulent[i] = *Succulent;
		spinup_Env[i] = *Env;
		spinup_Plot[i] = *Plot;
		spinup_Globals[i] = *Globals;

		if (UseDisturbances)
		{
			spinup_Globals[i].pat.use = grid_Disturb[cell].choices[0];
			spinup_Globals[i].mound.use = grid_Disturb[cell].choices[1];
			spinup_Globals[i].burrow.use = grid_Disturb[cell].choices[2];
		}
		if (UseSoils) //TODO: DM should i be cell value after lookup from soilTypes_Array
			_init_soil_layers(i, 1);

		spinup_SXW[i] = *SXW;

		spinup_SXW[i].f_roots = Str_Dup(SXW->f_roots);
		spinup_SXW[i].f_phen = Str_Dup(SXW->f_phen);
		spinup_SXW[i].f_bvt = Str_Dup(SXW->f_bvt);
		spinup_SXW[i].f_prod = Str_Dup(SXW->f_prod);
		spinup_SXW[i].f_watin = Str_Dup(SXW->f_watin);

		spinup_SXW[i].transpTotal = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
				sizeof(RealD), "_init_spinup_globals()");
		ForEachVegType(k) {
			spinup_SXW[i].transpVeg[k] = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
					sizeof(RealD), "_init_grid_globals()");
		}
		spinup_SXW[i].swc = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
				sizeof(RealF), "_init_grid_globals()");

		memcpy(spinup_SXW[i].transpTotal, SXW->transpTotal,
				SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
		ForEachVegType(k) {
			memcpy(spinup_SXW[i].transpVeg[k], SXW->transpVeg[k],
					SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
		}
		memcpy(spinup_SXW[i].swc, SXW->swc,
				SXW->NPds * SXW->NSoLyrs * sizeof(RealF));

		spinup_SW_Soilwat[i] = SW_Soilwat;
		spinup_SW_Site[i] = SW_Site;
		spinup_SW_VegProd[i] = SW_VegProd;

		spinup_SW_Site[i].lyr = Mem_Calloc(
				SW_Site.n_layers + SW_Site.deepdrain,
				sizeof(SW_LAYER_INFO *), "_init_grid_globals()");
		for (j = 0; j < SW_Site.n_layers + SW_Site.deepdrain; j++)
		{
			spinup_SW_Site[i].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
					"_init_grid_globals()");
			memcpy(spinup_SW_Site[i].lyr[j], SW_Site.lyr[j],
					sizeof(SW_LAYER_INFO));
		}
	}
	if (UseSoils)
		ChDir(".."); //get back to our previous directory

}

/***********************************************************/
static void _free_grid_globals(void)
{
	//frees memory allocated in _load_grid_globals() function.
	int i, j, k;
	GrpIndex c;
	SppIndex s;

	for (i = 0; i < grid_Cells; i++)
	{

		ForEachSpecies(s)
		{
			if (!Species[s]->use_me)
				continue;
			Mem_Free(grid_Species[s][i].kills);
			Mem_Free(grid_Species[s][i].seedprod);
			_free_head(grid_Species[s][i].IndvHead);
		}

		ForEachGroup(c)
			if (RGroup[c]->use_me)
				Mem_Free(grid_RGroup[c][i].kills);

		Mem_Free(grid_SXW[i].f_roots);
		Mem_Free(grid_SXW[i].f_phen);
		Mem_Free(grid_SXW[i].f_bvt);
		Mem_Free(grid_SXW[i].f_prod);
		Mem_Free(grid_SXW[i].f_watin);
		Mem_Free(grid_SXW[i].transpTotal);
		ForEachVegType(k) {
			Mem_Free(grid_SXW[i].transpVeg[k]);
		}
		Mem_Free(grid_SXW[i].swc);
		if (UseSoils)
		{
			Mem_Free(grid_SXW_ptrs[i].roots_max);
			Mem_Free(grid_SXW_ptrs[i].rootsXphen);
			Mem_Free(grid_SXW_ptrs[i].roots_active);
			Mem_Free(grid_SXW_ptrs[i].roots_active_rel);
			Mem_Free(grid_SXW_ptrs[i].roots_active_sum);
			Mem_Free(grid_SXW_ptrs[i].phen);
			Mem_Free(grid_SXW_ptrs[i].prod_bmass);
			Mem_Free(grid_SXW_ptrs[i].prod_pctlive);
		}

		for (j = 0;
				j < grid_SW_Site[i].n_layers + grid_SW_Site[i].deepdrain;
				j++)
			Mem_Free(grid_SW_Site[i].lyr[j]);
		Mem_Free(grid_SW_Site[i].lyr);

	}

}

/***********************************************************/
static void _free_spinup_globals(void)
{
	//frees memory allocated in _load_spinup_globals() function.
	int i, j, k;
	GrpIndex c;
	SppIndex s;

	for (i = 0; i < nSoilTypes; i++)
	{

		ForEachSpecies(s)
		{
			if (!Species[s]->use_me)
				continue;
			Mem_Free(spinup_Species[s][i].kills);
			Mem_Free(spinup_Species[s][i].seedprod);
			_free_head(spinup_Species[s][i].IndvHead);
		}

		ForEachGroup(c)
			if (RGroup[c]->use_me)
				Mem_Free(spinup_RGroup[c][i].kills);

		Mem_Free(spinup_SXW[i].f_roots);
		Mem_Free(spinup_SXW[i].f_phen);
		Mem_Free(spinup_SXW[i].f_bvt);
		Mem_Free(spinup_SXW[i].f_prod);
		Mem_Free(spinup_SXW[i].f_watin);
		Mem_Free(spinup_SXW[i].transpTotal);
		ForEachVegType(k) {
			Mem_Free(spinup_SXW[i].transpVeg[k]);
		}
		Mem_Free(spinup_SXW[i].swc);
		if (UseSoils)
		{
			Mem_Free(spinup_SXW_ptrs[i].roots_max);
			Mem_Free(spinup_SXW_ptrs[i].rootsXphen);
			Mem_Free(spinup_SXW_ptrs[i].roots_active);
			Mem_Free(spinup_SXW_ptrs[i].roots_active_rel);
			Mem_Free(spinup_SXW_ptrs[i].roots_active_sum);
			Mem_Free(spinup_SXW_ptrs[i].phen);
			Mem_Free(spinup_SXW_ptrs[i].prod_bmass);
			Mem_Free(spinup_SXW_ptrs[i].prod_pctlive);
		}

		for (j = 0;
				j < spinup_SW_Site[i].n_layers + spinup_SW_Site[i].deepdrain;
				j++)
			Mem_Free(spinup_SW_Site[i].lyr[j]);
		Mem_Free(spinup_SW_Site[i].lyr);
	}

}

/***********************************************************/
static void _free_grid_memory(void)
{
	//frees all the memory allocated in this file ST_Grid.c (most of it is dynamically allocated in _init_grid_globals() & _load_grid_globals() functions)

	int i;
	GrpIndex c;
	SppIndex s;

	_free_grid_globals();
	if (sd_Option2a || sd_Option2b)
		_free_spinup_memory();

	ForEachSpecies(s)
		if (Species[s]->use_me)
			Mem_Free(grid_Species[s]);
	ForEachGroup(c)
		if (RGroup[c]->use_me)
			Mem_Free(grid_RGroup[c]);

	Mem_Free(grid_Succulent);
	Mem_Free(grid_Env);
	Mem_Free(grid_Plot);
	Mem_Free(grid_Globals);

	Mem_Free(grid_SXW);
	Mem_Free(grid_SW_Soilwat);
	Mem_Free(grid_SW_Site);
	Mem_Free(grid_SW_VegProd);

	if (UseSoils)
	{
		free_all_sxw_memory();
		Mem_Free(grid_SXW_ptrs);
		for (i = 0; i < grid_Cells; i++)
			Mem_Free(grid_Soils[i].lyr);
		Mem_Free(grid_Soils);
	}
	if (UseDisturbances)
		Mem_Free(grid_Disturb);
	if (UseSeedDispersal)
	{
		ForEachSpecies(s)
			if (Species[s]->use_me && Species[s]->use_dispersal)
			{
				for (i = 0; i < grid_Cells; i++)
				{
					Mem_Free(grid_SD[s][i].cells);
					Mem_Free(grid_SD[s][i].prob);
					grid_SD[s][i].size = 0;
				}
				Mem_Free(grid_SD[s]);
			}

		for (i = 0; i < grid_Cells; i++)
			Mem_Free(grid_initSpecies[i].species_seed_avail);
		Mem_Free(grid_initSpecies);
	}

	stat_Free_Accumulators(); //free our memory we allocated for all the accumulators now that they're unnecessary to have

	for (i = 0; i < N_GRID_DIRECTORIES; i++) //frees the strings allocated in _init_grid_files()
		Mem_Free(grid_directories[i]);
	for (i = 0; i < N_GRID_FILES; i++)
		Mem_Free(grid_files[i]);

	// freeing random memory that other parts of steppe/soilwat allocate... this isn't quite everything but it's a good start
	parm_free_memory(); //frees memory allocated in ST_params.c

	ForEachSpecies(s)
	{
		if (!Species[s]->use_me)
			continue;
		_free_head(Species[s]->IndvHead);
		Mem_Free(Species[s]->kills);
		Mem_Free(Species[s]->seedprod);
	}

	ForEachGroup(c)
		if (RGroup[c]->use_me)
			Mem_Free(RGroup[c]->kills);

	for (i = 0; i < SW_Site.n_layers + SW_Site.deepdrain; i++)
		Mem_Free(SW_Site.lyr[i]);
	Mem_Free(SW_Site.lyr);

}

/***********************************************************/
static void _free_spinup_memory(void)
{
	// frees spinup memory

	GrpIndex c;
	SppIndex s;

	_free_spinup_globals();

	ForEachSpecies(s)
		if (Species[s]->use_me)
			Mem_Free(spinup_Species[s]);
	ForEachGroup(c)
		if (RGroup[c]->use_me)
			Mem_Free(spinup_RGroup[c]);

	Mem_Free(spinup_Succulent);
	Mem_Free(spinup_Env);
	Mem_Free(spinup_Plot);
	Mem_Free(spinup_Globals);

	Mem_Free(spinup_SXW);
	Mem_Free(spinup_SW_Soilwat);
	Mem_Free(spinup_SW_Site);
	Mem_Free(spinup_SW_VegProd);

	if (UseSoils)
	{
		Mem_Free(spinup_SXW_ptrs);
	}

	//if(sd_Option2a || sd_Option2b) {
	Mem_Free(soilTypes_Array);
	Mem_Free(grid_SoilTypes);
	//}
}

/***********************************************************/
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

	if (UseSoils)
		load_sxw_memory(grid_SXW_ptrs[cell].roots_max,
				grid_SXW_ptrs[cell].rootsXphen,
				grid_SXW_ptrs[cell].roots_active,
				grid_SXW_ptrs[cell].roots_active_rel,
				grid_SXW_ptrs[cell].roots_active_sum,
				grid_SXW_ptrs[cell].phen, grid_SXW_ptrs[cell].prod_bmass,
				grid_SXW_ptrs[cell].prod_pctlive);
}

/* Load gridCells[row][col] into the globals variables.
   Any call to this function should have an accompanying call to unload_cell(). */
static void load_cell(int row, int col){
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
	DuringSpinup = gridCells[row][col].duringSpinup;

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

/* Nullify all global variables. This function should appear after every call to load_cell to prevent
   accidental modification of a grid cell.
   Example usage:
   for i in rows{
       for j in columns{
	       load_cell(i, j)
		   //additional functions
	   }
   }
   unload_cell() */
static void unload_cell(){
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

/***********************************************************/
static void _load_spinup_cell(int cell)
{
	// loads the specified cell into the global variables (from the spinup)

	int j, k;
	GrpIndex c;
	SppIndex s;

	ForEachSpecies(s)
	{
		if (!Species[s]->use_me)
			continue;

		Mem_Free(Species[s]->kills);
		Mem_Free(Species[s]->seedprod);
		_free_head(Species[s]->IndvHead); //free_head() frees the memory allocated by the head and the memory allocated by each part of the linked list

		*Species[s] = spinup_Species[s][cell];

		Species[s]->kills = Mem_Calloc(spinup_Species[s][cell].max_age,
				sizeof(IntUS), "_load_spinup_cell(Species[s]->kills)");
		Species[s]->seedprod = Mem_Calloc(spinup_Species[s][cell].viable_yrs,
				sizeof(RealF), "_load_spinup_cell(Species[s]->seedprod)");

		memcpy(Species[s]->kills, spinup_Species[s][cell].kills,
				spinup_Species[s][cell].max_age * sizeof(IntUS));
		memcpy(Species[s]->seedprod, spinup_Species[s][cell].seedprod,
				spinup_Species[s][cell].viable_yrs * sizeof(RealF));
		Species[s]->IndvHead = _copy_head(spinup_Species[s][cell].IndvHead); //copy_head() deep copies the linked list structure (allocating memory when needed)... it will even allocate memory for the head of the list

	}

	ForEachGroup(c)
	{
		if (!RGroup[c]->use_me)
			continue;
		Mem_Free(RGroup[c]->kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)... we need to free it then reallocate it then memcpy it to get the deep copy we want

		*RGroup[c] = spinup_RGroup[c][cell]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (ie copy the values in the pointers instead of the addresses) the pointers.  A shallow copy will copy over the values for every non-pointer (C itself does not inherently know how to deep copy, so we must code this behaviour).

		RGroup[c]->kills = Mem_Calloc(spinup_RGroup[c][cell].max_age,
				sizeof(IntUS), "_load_spinup_cell(RGroup[c]->kills");
		memcpy(RGroup[c]->kills, spinup_RGroup[c][cell].kills,
				spinup_RGroup[c][cell].max_age * sizeof(IntUS));
	}

	*Succulent = spinup_Succulent[cell];
	*Env = spinup_Env[cell];
	*Plot = spinup_Plot[cell];
	*Globals = spinup_Globals[cell];

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

	*SXW = spinup_SXW[cell];
	SW_Site = spinup_SW_Site[cell];
	SW_Soilwat = spinup_SW_Soilwat[cell];
	SW_VegProd = spinup_SW_VegProd[cell];

	SXW->transpTotal = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealD),
			"_load_spinup_cell(SXW->transpTotal)");
	ForEachVegType(k) {
		SXW->transpVeg[k] = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealD),
				"_load_spinup_cell(SXW->transpVeg)");
	}
	SXW->swc = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealF),
			"_load_spinup_cell(SXW->swc)");

	SXW->f_roots = Str_Dup(spinup_SXW[cell].f_roots);
	SXW->f_phen = Str_Dup(spinup_SXW[cell].f_phen);
	SXW->f_bvt = Str_Dup(spinup_SXW[cell].f_bvt);
	SXW->f_prod = Str_Dup(spinup_SXW[cell].f_prod);
	SXW->f_watin = Str_Dup(spinup_SXW[cell].f_watin);
	memcpy(SXW->transpTotal, spinup_SXW[cell].transpTotal,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	ForEachVegType(k) {
		memcpy(SXW->transpVeg[k], spinup_SXW[cell].transpVeg[k],
				SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	}
	memcpy(SXW->swc, spinup_SXW[cell].swc,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealF));

	SW_Site.lyr = Mem_Calloc(
			spinup_SW_Site[cell].n_layers + spinup_SW_Site[cell].deepdrain,
			sizeof(SW_LAYER_INFO *), "_load_spinup_cell(SW_Site.lyr)");
	for (j = 0;
			j < spinup_SW_Site[cell].n_layers
							+ spinup_SW_Site[cell].deepdrain; j++)
	{
		SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
				"_load_spinup_cell(SW_Site.lyr[j])");
		memcpy(SW_Site.lyr[j], spinup_SW_Site[cell].lyr[j],
				sizeof(SW_LAYER_INFO));
	}

	if (UseSoils)
		load_sxw_memory(spinup_SXW_ptrs[cell].roots_max,
				spinup_SXW_ptrs[cell].rootsXphen,
				spinup_SXW_ptrs[cell].roots_active,
				spinup_SXW_ptrs[cell].roots_active_rel,
				spinup_SXW_ptrs[cell].roots_active_sum,
				spinup_SXW_ptrs[cell].phen,
				spinup_SXW_ptrs[cell].prod_bmass,
				spinup_SXW_ptrs[cell].prod_pctlive);
}

/***********************************************************/
static void _save_cell(int row, int col, int year, Bool useAccumulators)
{
	// saves the specified cell into the grid variables

	int cell = col + ((row - 1) * grid_Cols) - 1; // converts the row/col into an array index
	int j, k;
	GrpIndex c;
	SppIndex s;
	//fprintf(stderr, "saving cell: %d\n", cell);

	if (useAccumulators)
		stat_Save_Accumulators(cell, year);

	ForEachSpecies(s)
	{
		if (!Species[s]->use_me)
			continue;

		Mem_Free(grid_Species[s][cell].kills);
		Mem_Free(grid_Species[s][cell].seedprod);
		_free_head(grid_Species[s][cell].IndvHead);

		grid_Species[s][cell] = *Species[s];

		grid_Species[s][cell].kills = Mem_Calloc(Species[s]->max_age,
				sizeof(IntUS), "_save_cell(grid_Species[cell][s].kills)");
		grid_Species[s][cell].seedprod = Mem_Calloc(Species[s]->viable_yrs,
				sizeof(RealF), "_save_cell(grid_Species[cell][s].seedprod)");

		memcpy(grid_Species[s][cell].kills, Species[s]->kills,
				Species[s]->max_age * sizeof(IntUS));
		memcpy(grid_Species[s][cell].seedprod, Species[s]->seedprod,
				Species[s]->viable_yrs * sizeof(RealF));
		grid_Species[s][cell].IndvHead = _copy_head(Species[s]->IndvHead);
	}

	ForEachGroup(c)
	{
		if (!RGroup[c]->use_me)
			continue;
		Mem_Free(grid_RGroup[c][cell].kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)

		grid_RGroup[c][cell] = *RGroup[c]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (i.e. copy the values in the pointers instead of the addresses) the pointers

		grid_RGroup[c][cell].kills = Mem_Calloc(RGroup[c]->max_age,
				sizeof(IntUS), "_save_cell(grid_RGroup[cell][c].kills)");
		memcpy(grid_RGroup[c][cell].kills, RGroup[c]->kills,
				RGroup[c]->max_age * sizeof(IntUS));
	}

	grid_Succulent[cell] = *Succulent;
	grid_Env[cell] = *Env;
	grid_Plot[cell] = *Plot;
	grid_Globals[cell] = *Globals;

	Mem_Free(grid_SXW[cell].f_roots);
	Mem_Free(grid_SXW[cell].f_phen);
	Mem_Free(grid_SXW[cell].f_bvt);
	Mem_Free(grid_SXW[cell].f_prod);
	Mem_Free(grid_SXW[cell].f_watin);
	Mem_Free(grid_SXW[cell].transpTotal);
	ForEachVegType(k) {
		Mem_Free(grid_SXW[cell].transpVeg[k]);
	}
	Mem_Free(grid_SXW[cell].swc);
	for (j = 0;
			j < grid_SW_Site[cell].n_layers + grid_SW_Site[cell].deepdrain;
			j++)
		Mem_Free(grid_SW_Site[cell].lyr[j]);
	Mem_Free(grid_SW_Site[cell].lyr);

	grid_SXW[cell] = *SXW;
	grid_SW_Site[cell] = SW_Site;
	grid_SW_Soilwat[cell] = SW_Soilwat;
	grid_SW_VegProd[cell] = SW_VegProd;

	grid_SXW[cell].transpTotal = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
			sizeof(RealD), "_save_cell(grid_SXW[cell].transp)");
	ForEachVegType(k) {
		grid_SXW[cell].transpVeg[k] = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
				sizeof(RealD), "_save_cell(grid_SXW[cell].transp)");
	}
	grid_SXW[cell].swc = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealF),
			"_save_cell(grid_SXW[cell].swc)");

	grid_SXW[cell].f_roots = Str_Dup(SXW->f_roots);
	grid_SXW[cell].f_phen = Str_Dup(SXW->f_phen);
	grid_SXW[cell].f_bvt = Str_Dup(SXW->f_bvt);
	grid_SXW[cell].f_prod = Str_Dup(SXW->f_prod);
	grid_SXW[cell].f_watin = Str_Dup(SXW->f_watin);
	memcpy(grid_SXW[cell].transpTotal, SXW->transpTotal,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	ForEachVegType(k) {
		memcpy(grid_SXW[cell].transpVeg[k], SXW->transpVeg[k],
				SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	}
	memcpy(grid_SXW[cell].swc, SXW->swc,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealF));

	grid_SW_Site[cell].lyr = Mem_Calloc(
			SW_Site.n_layers + SW_Site.deepdrain, sizeof(SW_LAYER_INFO *),
			"_save_cell(grid_SW_Site[cell].lyr[j])");
	for (j = 0; j < SW_Site.n_layers + SW_Site.deepdrain; j++)
	{
		grid_SW_Site[cell].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
				"_save_cell(grid_SW_Site[cell].lyr[j])");
		memcpy(grid_SW_Site[cell].lyr[j], SW_Site.lyr[j],
				sizeof(SW_LAYER_INFO));
	}

	if (UseSoils)
		save_sxw_memory(grid_SXW_ptrs[cell].roots_max,
				grid_SXW_ptrs[cell].rootsXphen,
				grid_SXW_ptrs[cell].roots_active,
				grid_SXW_ptrs[cell].roots_active_rel,
				grid_SXW_ptrs[cell].roots_active_sum,
				grid_SXW_ptrs[cell].phen, grid_SXW_ptrs[cell].prod_bmass,
				grid_SXW_ptrs[cell].prod_pctlive);
}

/***********************************************************/
static void _save_spinup_cell(int cell)
{
	// saves the specified cell into the grid variables (from the spinup)

	int j, k;
	GrpIndex c;
	SppIndex s;

	ForEachSpecies(s)
	{
		if (!Species[s]->use_me)
			continue;

		Mem_Free(spinup_Species[s][cell].kills);
		Mem_Free(spinup_Species[s][cell].seedprod);
		_free_head(spinup_Species[s][cell].IndvHead);

		spinup_Species[s][cell] = *Species[s];

		spinup_Species[s][cell].kills = Mem_Calloc(Species[s]->max_age,
				sizeof(IntUS), "_save_cell(grid_Species[cell][s].kills)");
		spinup_Species[s][cell].seedprod = Mem_Calloc(Species[s]->viable_yrs,
				sizeof(RealF), "_save_cell(grid_Species[cell][s].seedprod)");

		memcpy(spinup_Species[s][cell].kills, Species[s]->kills,
				Species[s]->max_age * sizeof(IntUS));
		memcpy(spinup_Species[s][cell].seedprod, Species[s]->seedprod,
				Species[s]->viable_yrs * sizeof(RealF));
		spinup_Species[s][cell].IndvHead = _copy_head(Species[s]->IndvHead);
	}

	ForEachGroup(c)
	{
		if (!RGroup[c]->use_me)
			continue;
		Mem_Free(spinup_RGroup[c][cell].kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)

		spinup_RGroup[c][cell] = *RGroup[c]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (i.e. copy the values in the pointers instead of the addresses) the pointers

		spinup_RGroup[c][cell].kills = Mem_Calloc(RGroup[c]->max_age,
				sizeof(IntUS), "_save_cell(grid_RGroup[cell][c].kills)");
		memcpy(spinup_RGroup[c][cell].kills, RGroup[c]->kills,
				RGroup[c]->max_age * sizeof(IntUS));
	}

	spinup_Succulent[cell] = *Succulent;
	spinup_Env[cell] = *Env;
	spinup_Plot[cell] = *Plot;
	spinup_Globals[cell] = *Globals;

	Mem_Free(spinup_SXW[cell].f_roots);
	Mem_Free(spinup_SXW[cell].f_phen);
	Mem_Free(spinup_SXW[cell].f_bvt);
	Mem_Free(spinup_SXW[cell].f_prod);
	Mem_Free(spinup_SXW[cell].f_watin);
	Mem_Free(spinup_SXW[cell].transpTotal);
	ForEachVegType(k) {
		Mem_Free(spinup_SXW[cell].transpVeg[k]);
	}
	Mem_Free(spinup_SXW[cell].swc);
	for (j = 0;
			j
					< spinup_SW_Site[cell].n_layers
							+ spinup_SW_Site[cell].deepdrain; j++)
		Mem_Free(spinup_SW_Site[cell].lyr[j]);
	Mem_Free(spinup_SW_Site[cell].lyr);

	spinup_SXW[cell] = *SXW;
	spinup_SW_Site[cell] = SW_Site;
	spinup_SW_Soilwat[cell] = SW_Soilwat;
	spinup_SW_VegProd[cell] = SW_VegProd;

	spinup_SXW[cell].transpTotal = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
			sizeof(RealD), "_save_cell(grid_SXW[cell].transp)");
	ForEachVegType(k) {
		spinup_SXW[cell].transpVeg[k] = Mem_Calloc(SXW->NPds * SXW->NSoLyrs,
				sizeof(RealD), "_save_cell(grid_SXW[cell].transp)");
	}
	spinup_SXW[cell].swc = Mem_Calloc(SXW->NPds * SXW->NSoLyrs, sizeof(RealF),
			"_save_cell(grid_SXW[cell].swc)");

	spinup_SXW[cell].f_roots = Str_Dup(SXW->f_roots);
	spinup_SXW[cell].f_phen = Str_Dup(SXW->f_phen);
	spinup_SXW[cell].f_bvt = Str_Dup(SXW->f_bvt);
	spinup_SXW[cell].f_prod = Str_Dup(SXW->f_prod);
	spinup_SXW[cell].f_watin = Str_Dup(SXW->f_watin);
	memcpy(spinup_SXW[cell].transpTotal, SXW->transpTotal,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	ForEachVegType(k) {
		memcpy(spinup_SXW[cell].transpVeg[k], SXW->transpVeg[k],
				SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	}
	memcpy(spinup_SXW[cell].swc, SXW->swc,
			SXW->NPds * SXW->NSoLyrs * sizeof(RealF));

	spinup_SW_Site[cell].lyr = Mem_Calloc(
			SW_Site.n_layers + SW_Site.deepdrain, sizeof(SW_LAYER_INFO *),
			"_save_cell(grid_SW_Site[cell].lyr[j])");
	for (j = 0; j < SW_Site.n_layers + SW_Site.deepdrain; j++)
	{
		spinup_SW_Site[cell].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO),
				"_save_cell(grid_SW_Site[cell].lyr[j])");
		memcpy(spinup_SW_Site[cell].lyr[j], SW_Site.lyr[j],
				sizeof(SW_LAYER_INFO));
	}

	if (UseSoils)
		save_sxw_memory(spinup_SXW_ptrs[cell].roots_max,
				spinup_SXW_ptrs[cell].rootsXphen,
				spinup_SXW_ptrs[cell].roots_active,
				spinup_SXW_ptrs[cell].roots_active_rel,
				spinup_SXW_ptrs[cell].roots_active_sum,
				spinup_SXW_ptrs[cell].phen,
				spinup_SXW_ptrs[cell].prod_bmass,
				spinup_SXW_ptrs[cell].prod_pctlive);
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

/***********************************************************/
static void _read_disturbances_in(void)
{
	// reads the grid disturbances input file
	// the file should be something like: "cell,use_fecal_pats,use_ant_mounds,use_animal_burrows,kill_yr"
	// there should be no spaces in between, just commas separating the values
	// kill_yr will overwrite the kill year for each RGroup in the cell (0 means don't use, a # > 0 means kill everything at this year)

	FILE *f;
	char buf[1024];
	int i, row, col, cell, num = 0, choices[3];
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

        // xgrow and veg_prod_type not found in grid_disturbances.csv
        ForEachGroup(rg)
            num = sscanf(buf, "%d,%d,%d,%d,%hd,%f,%hd,%hd,%hu,%hd", &cell,
                &choices[0], &choices[1],
                &choices[2], &RGroup[rg]->killyr,
                &RGroup[rg]->killfreq, &RGroup[rg]->extirp,
                &RGroup[rg]->killfreq_startyr,
                &RGroup[rg]->grazingfrq, &RGroup[rg]->grazingfreq_startyr);

        // Comes from _load_grid_globals
        if (UseDisturbances) {
            Globals->pat.use = choices[0];
            Globals->mound.use = choices[1];
            Globals->burrow.use = choices[2];
        }

        // xgrow and veg_prod_type not found in grid_disturbances.csv
		if (num != 10)
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

/***********************************************************/
static void _read_soils_in(void)
{
	// reads the grid soils input
	// the file should be something like: "cell,copy_cell,copy_which,num_layers,..."
	// there should be no spaces in between, just commas separating the values
	// this function reads in pretty much a .csv file, but it will not account for all of the possibilities that a .csv file could be written as (as accounting for all these possibilities would take a while to code and be unproductive) so keep that in mind

	FILE *f;
	char buf[4096];
	char rootsin[20];
	char seps[] = ",";
	char *token;
	int i, j, k, cell, num, do_copy, copy_cell, num_layers, depth, depthMin,
			stringIndex, row, col, copy_cell_row, copy_cell_col;
	float d[11];

	if (sd_Option2a || sd_Option2b)
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
			if (sd_Option2a || sd_Option2b)
			    grid_SoilTypes[cell] = grid_SoilTypes[copy_cell];

			strcpy(gridCells[row][col].mySoils.rootsFile, gridCells[copy_cell_row][copy_cell_col].mySoils.rootsFile);

			continue;
		}
		else if (do_copy == 1)
			LogError(logfp, LOGFATAL,
					"Invalid %s file line %d invalid copy_cell attempt",
					grid_files[GRID_FILE_SOILS], i + 2);

        // TODO: figure out how to change this code for our new struct
		if (sd_Option2a || sd_Option2b)
		{
			grid_SoilTypes[cell] = nSoilTypes;
			soilTypes_Array[nSoilTypes] = cell;
			nSoilTypes++;
		}

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

/***********************************************************/
static void _init_soil_layers(int cell, int isSpinup)
{
	// initializes the soilwat soil layers for the cell correctly based upon the input gathered from our grid_soils input file
	// pretty much takes the data from grid_Soils (read in in _read_soils_in()) and converts it to what SW_Site needs...
	// this function does generally the same things that the _read_layers() function in SW_Site.c does, except that it does it in a way that lets us use it in the grid...
	int i, j;
	i = cell;
	char* errtype;
	Bool evap_ok = TRUE, transp_ok_forb = TRUE, transp_ok_tree = TRUE,
			transp_ok_shrub = TRUE, transp_ok_grass = TRUE; /* mitigate gaps in layers */
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

	if (!isSpinup)
	{

		grid_SXW_ptrs[i].roots_max = Mem_Calloc(SXW->NGrps * SXW->NTrLyrs,
				sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].rootsXphen = Mem_Calloc(
				SXW->NGrps * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		grid_SXW_ptrs[i].roots_active = Mem_Calloc(
				SXW->NGrps * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		grid_SXW_ptrs[i].roots_active_rel = Mem_Calloc(
				SXW->NGrps * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		grid_SXW_ptrs[i].roots_active_sum = Mem_Calloc(
				4 * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		grid_SXW_ptrs[i].phen = Mem_Calloc(SXW->NGrps * MAX_MONTHS,
				sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].prod_bmass = Mem_Calloc(SXW->NGrps * MAX_MONTHS,
				sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].prod_pctlive = Mem_Calloc(SXW->NGrps * MAX_MONTHS,
				sizeof(RealD), "_init_soil_layers()");

		save_sxw_memory(grid_SXW_ptrs[i].roots_max, grid_SXW_ptrs[i].rootsXphen,
				grid_SXW_ptrs[i].roots_active,
				grid_SXW_ptrs[i].roots_active_rel,
				grid_SXW_ptrs[i].roots_active_sum, grid_SXW_ptrs[i].phen,
				grid_SXW_ptrs[i].prod_bmass, grid_SXW_ptrs[i].prod_pctlive);
	}
	else
	{

		spinup_SXW_ptrs[i].roots_max = Mem_Calloc(SXW->NGrps * SXW->NTrLyrs,
				sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].rootsXphen = Mem_Calloc(
				SXW->NGrps * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		spinup_SXW_ptrs[i].roots_active = Mem_Calloc(
				SXW->NGrps * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		spinup_SXW_ptrs[i].roots_active_rel = Mem_Calloc(
				SXW->NGrps * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		spinup_SXW_ptrs[i].roots_active_sum = Mem_Calloc(
				4 * SXW->NPds * SXW->NTrLyrs, sizeof(RealD),
				"_init_soil_layers()");
		spinup_SXW_ptrs[i].phen = Mem_Calloc(SXW->NGrps * MAX_MONTHS,
				sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].prod_bmass = Mem_Calloc(SXW->NGrps * MAX_MONTHS,
				sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].prod_pctlive = Mem_Calloc(SXW->NGrps * MAX_MONTHS,
				sizeof(RealD), "_init_soil_layers()");

		save_sxw_memory(spinup_SXW_ptrs[i].roots_max,
				spinup_SXW_ptrs[i].rootsXphen, spinup_SXW_ptrs[i].roots_active,
				spinup_SXW_ptrs[i].roots_active_rel,
				spinup_SXW_ptrs[i].roots_active_sum, spinup_SXW_ptrs[i].phen,
				spinup_SXW_ptrs[i].prod_bmass, spinup_SXW_ptrs[i].prod_pctlive);
	}
}

/***********************************************************/
static float _read_a_float(FILE *f, char *buf, const char *filename,
		const char *descriptor)
{
	//small function to reduce code duplication in the _read_seed_dispersal_in() function...
	//f should be already open, and all of the character arrays should be pre-allocated before calling the function...
	float result;

	if (!GetALine(f, buf))
		LogError(logfp, LOGFATAL, "Invalid %s file: %s", filename, descriptor);
	if (sscanf(buf, "%f", &result) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: %s", filename, descriptor);

	return result;
}

/***********************************************************/
static float _cell_dist(int row1, int row2, int col1, int col2, float cellLen)
{
	//returns the distance between the two grid cells
	if (row1 == row2)
	{
		return (abs(col1-col2) * cellLen);
	}
	else if (col1 == col2)
	{
		return (abs(row1-row2) * cellLen);
	}
	else
	{ // row1 != row2 && col1 != col2
		//the problem can be thought of in terms of a right triangle...
		//using the pythagorean theorem: c = sqrt(a^2 + b^2)... c (the hypotenuse) represents the distance that we need.  a is the distance between columns and b is the distance between rows.
		return sqrt(
				pow(abs(col1-col2)*cellLen, 2.0) + pow(abs(row1-row2)*cellLen, 2.0));
	}
}

/***********************************************************/
static void _read_seed_dispersal_in(void)
{
	// reads the grid seed dispersal input file and sets up grid_SD with the correct values and probabilities

	FILE *f;
	char buf[1024];
	float sd_Rate, H, VW, VT, MAXD, plotLength, d, pd;
	int maxCells, i, j, k, MAXDP, row, col, cell;
	SppIndex s;

	// read in the seed dispersal input file to get the constants that we need
	f = OpenFile(grid_files[GRID_FILE_SEED_DISPERSAL], "r");

	VW = _read_a_float(f, buf, grid_files[GRID_FILE_SEED_DISPERSAL], "VW line");

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_DoOutput) != 1)
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal output line\n", grid_files[GRID_FILE_SEED_DISPERSAL]);

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_MakeHeader) != 1)
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal make header line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	GetALine(f, buf);
	if (sscanf(buf, "%c", &sd_Sep) != 1)
		LogError(logfp, LOGFATAL,
				"Invalid %s file: seed dispersal seperator line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	if (sd_Sep == 't') //dealing with tab and space special cases...
		sd_Sep = '\t';
	else if (sd_Sep == 's')
		sd_Sep = ' ';

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_NYearsSeedsAvailable) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 1 line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_Option1a) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 1a line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_Option1b) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 1b line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_Option2a) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 2a line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	GetALine(f, buf);
	if (sscanf(buf, "%d", &sd_Option2b) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 2b line\n",
				grid_files[GRID_FILE_SEED_DISPERSAL]);

	if ((sd_Option1a && (sd_Option1b || sd_Option2a || sd_Option2b))
			|| (sd_Option2a && (sd_Option1a || sd_Option1b || sd_Option2b)))
		LogError(logfp, LOGFATAL,
				"Invalid %s file: conflicting options chosen\n", grid_files[GRID_FILE_SEED_DISPERSAL]);

	CloseFile(&f);

    for (i = 0; i < grid_Cells; i++) {
        row = i / grid_Cols;
        col = i % grid_Cols;

        gridCells[row][col].mySeedDispersal = Mem_Calloc(MAX_SPECIES, sizeof(Grid_SD_St), "_read_seed_dispersal_in");
    }

	/*
	 * The following variables should be global (not cell-specific):
	 *     Globals->sppCount
	 *     Species[s]->sd_H
	 *     Species[s]->sd_VT
	 *     Globals->plotsize
	 *     Species[s]->use_me
	 *     Species[s]->use_dispersal
	 *
	 * Since the values of these variables do not change between cells, we can call load_cell(0, 0)
	 * to load their values so that ForEachSpecies and other lines of code do not segfault.
	 */
	load_cell(0, 0);

	// Begin reading cell-specific fields
	ForEachSpecies(s)
	{
		// set up grid_SD with the seed dispersal probabilities needed later on...
		H = Species[s]->sd_H;
		VT = Species[s]->sd_VT;
		MAXD = ((H * VW) / VT) / 100.0; // divide by 100.0 because we want it in meters, not centimeters
		sd_Rate = -(log(0.05) / MAXD); //sd_Rate is the seed dispersal rate... 0.05 = exp(-RATE*MAXD) => RATE = -(ln(0.05)/MAXD) See Coffin et al. 1993

		plotLength = sqrt(Globals->plotsize);
		MAXDP = (int) ceil(MAXD / plotLength); //MAXD in terms of plots... rounds up to the nearest integer
		maxCells = (int) pow((MAXDP * 2) + 1.0, 2.0); //gets the maximum number of cells that a grid cell can possibly disperse seeds to... it ends up being more then the maximum actually...
		if (grid_Cells < maxCells)
			maxCells = grid_Cells;
		if (!(Species[s]->use_me && Species[s]->use_dispersal))
			continue;

		for (i = 0; i < grid_Cells; i++)
		{
		    row = i / grid_Cols;
		    col = i % grid_Cols;

			gridCells[row][col].mySeedDispersal[s].cells = Mem_Calloc(maxCells, sizeof(int),
					"_read_seed_dispersal_in()"); //the cell number
			gridCells[row][col].mySeedDispersal[s].prob = Mem_Calloc(maxCells, sizeof(float),
					"_read_seed_dispersal_in()"); //the probability that the cell will disperse seeds to this distance
			gridCells[row][col].mySeedDispersal[s].size = 0; //refers to the number of cells reachable...
		}

		for (row = 1; row <= grid_Rows; row++)
			for (col = 1; col <= grid_Cols; col++)
			{

				cell = col + ((row - 1) * grid_Cols) - 1;
				k = 0;

				for (i = 1; i <= grid_Rows; i++)
					for (j = 1; j <= grid_Cols; j++)
					{
						if (i == row && j == col)
							continue;

						d = _cell_dist(i, row, j, col, plotLength); //distance
						pd = (d > MAXD) ? (0.0) : (exp(-sd_Rate * d)); //dispersal probability

						if (!ZRO(pd))
						{
							gridCells[row - 1][col - 1].mySeedDispersal[s].cells[k] = i
									+ ((j - 1) * grid_Cols) - 1;
							gridCells[row - 1][col - 1].mySeedDispersal[s].prob[k] = pd;
							gridCells[row - 1][col - 1].mySeedDispersal[s].size++;
							k++;
						}
					}
			}

		for (i = 0; i < grid_Cells; i++) {
		    row = i / grid_Cols;
		    col = i % grid_Cols;

			if (gridCells[row][col].mySeedDispersal[s].size > 0)
			{
                gridCells[row][col].mySeedDispersal[s].cells = Mem_ReAlloc(gridCells[row][col].mySeedDispersal[s].cells,
                        gridCells[row][col].mySeedDispersal[s].size * sizeof(int));
                gridCells[row][col].mySeedDispersal[s].prob = Mem_ReAlloc(gridCells[row][col].mySeedDispersal[s].prob,
                        gridCells[row][col].mySeedDispersal[s].size * sizeof(float));
			}
		}
	}

	unload_cell();
}

/***********************************************************/
static void _do_seed_dispersal(void)
{
	float biomass, randomN, LYPPT, presentProb, receivedProb;
	int i, j, germ, sgerm, year;
	SppIndex s;

	if (Globals->currYear == 1 && !sd_Option1a && !sd_Option1b)
	{ //since we have no previous data to go off of, use the current years...
		for (i = 0; i < grid_Cells; i++)
			ForEachSpecies(s)
			{
				if (!(Species[s]->use_me && Species[s]->use_dispersal))
					continue;
				grid_Species[s][i].allow_growth = grid_Species[s][i].sd_sgerm =
						1;// since it's the first year, we have to allow growth...
				if (UseDisturbances)
					if (1 == grid_Disturb[i].kill_yr)
						grid_Species[s][i].allow_growth = 0;
				grid_SD[s][i].lyppt = grid_Env[i].ppt;
			}
	}
	else
	{
		// figure out whether or not to allow growth for the current year... based upon whether the species already has plants or germination allowed this year and seeds received last year...
		ForEachSpecies(s)
		{

			if (!(Species[s]->use_me && Species[s]->use_dispersal))
				continue;

			// germination probability
			randomN = RandUni(&grid_rng);
			germ = LE(randomN, Species[s]->seedling_estab_prob);

			year = Globals->currYear - 1;

			for (i = 0; i < grid_Cells; i++)
			{

				if (sd_Option1a && Globals->currYear <= sd_NYearsSeedsAvailable)
				{
					grid_SD[s][i].seeds_present = 1;
				}
				else if (sd_Option1b
						&& Globals->currYear <= sd_NYearsSeedsAvailable
						&& grid_initSpecies[i].species_seed_avail[s])
				{
					grid_SD[s][i].seeds_present = 1;
				}

				sgerm = (grid_SD[s][i].seeds_present
						|| grid_SD[s][i].seeds_received) && germ; //refers to whether the species has seeds available from the previous year and conditions are correct for germination this year
				grid_Species[s][i].allow_growth = FALSE;
				biomass = grid_Species[s][i].relsize
						* grid_Species[s][i].mature_biomass;

				if (UseDisturbances)
				{
					if ((sgerm || year < grid_Disturb[i].kill_yr
							|| grid_Disturb[i].kill_yr <= 0 || GT(biomass, 0.0))
					/*&& (year != grid_Disturb[i].kill_yr)*/)
					{
						//commented above one condition as it was causing a bug there, next year of killing year will make
						//allow_growth flag to false as 	year = Globals->currYear - 1 , so for example if killing year= 6 and Globals->currYear=7 then here
						// year variable will be 7-1 =6 that is equal to killing year 6, so this condition (year != grid_Disturb[i].kill_yr)
						//will fail and allow_growth will not become TRUE, then when Globals->currYear=8 this allow_growth= FALSE will carry forward and there will no call
						// to other functions like Species_Update_Newsize() so new size will not be updated and last year size will carry forward so in final output year 7 and year 8 will
						// have same output that is not correct.
						grid_Species[s][i].allow_growth = TRUE;
					}

				}
				else if (sgerm || GT(biomass, 0.0))
					grid_Species[s][i].allow_growth = TRUE;
				grid_Species[s][i].sd_sgerm = sgerm; //based upon whether we have received/produced seeds that germinated
				//if(grid_Species[s][i].allow_growth == TRUE &&  i == 52 && s == 0 && Globals->currIter == 1)
				//	printf("%s allow_growth:%d year:%d sgerm:%d iter:%d\n", grid_Species[s][i].name, grid_Species[s][i].allow_growth, year, sgerm, Globals->currIter);
			}
		}

	}

	// calculate whether or not seeds were received/produced this year, this data is used the next time the function is called
	ForEachSpecies(s)
	{
		if (!(Species[s]->use_me && Species[s]->use_dispersal))
			continue;

		IndivType* indiv;

		// figure out which species in each cell produced seeds...
		for (i = 0; i < grid_Cells; i++)
		{
			grid_SD[s][i].seeds_present = grid_SD[s][i].seeds_received =
					grid_Species[s][i].received_prob = 0;

			biomass = 0;	//getting the biggest individual in the species...
			ForEachIndiv(indiv, &grid_Species[s][i])
				if (indiv->relsize * grid_Species[s][i].mature_biomass
						> biomass)
					biomass = indiv->relsize
							* grid_Species[s][i].mature_biomass;

			if (GE(biomass,
					grid_Species[s][i].mature_biomass
							* grid_Species[s][i].sd_Param1))
			{
				randomN = RandUni(&grid_rng);

				LYPPT = grid_SD[s][i].lyppt;
				float PPTdry = grid_Species[s][i].sd_PPTdry, PPTwet =
						grid_Species[s][i].sd_PPTwet;
				float Pmin = grid_Species[s][i].sd_Pmin, Pmax =
						grid_Species[s][i].sd_Pmax;

				//p3 = Pmin, if LYPPT < PPTdry
				//p3 = 1 - (1-Pmin) * exp(-d * (LYPPT - PPTdry)) with d = - ln((1 - Pmax)/(1 - Pmin)) / (PPTwet - PPTdry), if PPTdry <= LYPPT <= PPTwet
				//p3 = Pmax, if LYPPT > PPTwet

				presentProb = 0.0;
				if (PPTdry <= LYPPT && LYPPT <= PPTwet)
				{
					float d = -log(((1 - Pmax) / (1 - Pmin)))
							/ (PPTwet - PPTdry); //log is the natural log in STD c's math.h
					presentProb = 1 - (1 - Pmin) * exp((-d * (LYPPT - PPTdry)));
				}
				else if (LYPPT < PPTdry)
					presentProb = Pmin;
				else if (LYPPT > PPTwet)
					presentProb = Pmax;

				if (LE(randomN, presentProb))
					grid_SD[s][i].seeds_present = 1;
			}
			//if(i == 0) printf("cell: %d lyppt: %f\n", i, grid_SD[i].lyppt);
		}

		// figure out which species in each cell received seeds...
		for (i = 0; i < grid_Cells; i++)
		{
			if (grid_SD[s][i].seeds_present)
				continue;
			receivedProb = 0;

			for (j = 0; j < grid_SD[s][i].size; j++)
				if (grid_SD[s][grid_SD[s][i].cells[j]].seeds_present)
					receivedProb += grid_SD[s][i].prob[j];

			randomN = RandUni(&grid_rng);
			if (LE(randomN, receivedProb) && !ZRO(receivedProb))
				grid_SD[s][i].seeds_received = 1;
			else
				grid_SD[s][i].seeds_received = 0;

			grid_Species[s][i].received_prob = receivedProb;
		}
	}
}

/***********************************************************/
static void _set_sd_lyppt(int row, int col)
{
	int cell = col + ((row - 1) * grid_Cols) - 1;
	SppIndex s;

	ForEachSpecies(s)
		if (Species[s]->use_me && Species[s]->use_dispersal)
			grid_SD[s][cell].lyppt = grid_Env[cell].ppt;
}


/*
 *
 */

static void _do_groups_and_species_extirpate(void)
{
	//printf("inside _do_groups_and_species_extirpate()\n");
	GrpIndex rg;
	ForEachGroup(rg)
	{
		rgroup_Extirpate(rg);
	}
}

/***********************************************************/
static void _kill_groups_and_species(void)
{
	/* (AKT /Kyle) 11/17/2016
	 *we may need to revisit this setting PR to 0 when we kill everything
	 *ForEachGroup(c)
	 *RGroup[c]->pr = 0.0; // reset the pr, so our output doesn't look weird
	 */

	//printf("inside _kill_groups_and_species()\n");
	GrpIndex rg;
	ForEachGroup(rg)
	{
		if (Globals->currYear < RGroup[rg]->startyr)
		{
			/* don't start trying to kill until RGroup[rg]->startyr year as nothing grow till now */
			continue;
		}
		Int i;

		ForEachEstSpp2( rg, i)
		{
			if (!Species[RGroup[rg]->est_spp[i]]->use_me)
			{
				continue;
			}
			else
			{
			//	printf("calling Species_Proportion_Kill() with rgroup name= %s , RGroup[%d]->proportion_killed =%f for Species[%d]->name= %s \n",RGroup[rg]->name,rg, RGroup[rg]->proportion_killed, i, Species[RGroup[rg]->est_spp[i]]->name);
				Species_Proportion_Kill(RGroup[rg]->est_spp[i], 6, RGroup[rg]->proportion_killed);
			}

		}

	}


}

/***********************************************************/
static int _do_grid_disturbances(int row, int col)
{
	// return 1 if a disturbance occurs, else return 0

	if(UseDisturbances)
	{
		int cell = col + ((row - 1) * grid_Cols) - 1;
//		printf( "inside _do_grid_disturbances Globals->currYear =%d, cell=%d, grid_Disturb[cell].kill_yr =%d \n",
//							Globals->currYear, cell, grid_Disturb[cell].kill_yr);
		if ((Globals->currYear >=grid_Disturb[cell].killfreq_startyr) && GT((float)grid_Disturb[cell].killfrq, 0.))
		{
			if (LT((float)grid_Disturb[cell].killfrq, 1.0))
			{
				if (RandUni(&grid_rng) <= grid_Disturb[cell].killfrq)
				{
					grid_Disturb[cell].kill_yr = Globals->currYear;
				}

			}
			else if (((Globals->currYear - grid_Disturb[cell].killfreq_startyr) % (IntU) grid_Disturb[cell].killfrq) == 0)
			{
				grid_Disturb[cell].kill_yr = Globals->currYear;
			}

		}

		if (Globals->currYear == grid_Disturb[cell].extirp)
		{
			_do_groups_and_species_extirpate();
			return 1;
		}
		else if (Globals->currYear == grid_Disturb[cell].kill_yr)
		{
//			printf( "current year matched with cell kill_year so calling  _kill_groups_and_species() Globals->currYear =%d, cell=%d, grid_Disturb[cell].kill_yr =%d \n",
//							Globals->currYear, cell, grid_Disturb[cell].kill_yr);
			_kill_groups_and_species();
			return 1;
		}


	}
	return 0;
}


static void _do_grid_proportion_Recovery(int row, int col)
{
	/*======================================================*/
	/* PURPOSE
       Perform the sorts of proportion_Recovery one might expect at next year after the killing year
	   HISTORY
	   Nov 22 2016 -AKT -Added Species Proportion Recovery for Grid Version */
	/*======================================================*/
	if (UseDisturbances)
	{
		int cell = col + ((row - 1) * grid_Cols) - 1;
//		printf( "inside _do_grid_proportion_Recovery Globals->currYear =%d, cell=%d, grid_Disturb[cell].kill_yr =%d \n",
//							Globals->currYear, cell, grid_Disturb[cell].kill_yr);
		if ((Globals->currYear >= grid_Disturb[cell].killfreq_startyr) && GT((float)grid_Disturb[cell].killfrq, 0.))
		{
			if (LT((float)grid_Disturb[cell].killfrq, 1.0))
			{
				if (RandUni(&grid_rng) <= grid_Disturb[cell].killfrq)
				{
					grid_Disturb[cell].kill_yr = Globals->currYear;
				}

			}
			else if (((Globals->currYear - grid_Disturb[cell].killfreq_startyr) % (IntU) grid_Disturb[cell].killfrq) == 0)
			{
				grid_Disturb[cell].kill_yr = Globals->currYear;
			}

		}

		//rgroup proportion recovery
		if (Globals->currYear == grid_Disturb[cell].kill_yr)
		{
			GrpIndex rg;
			ForEachGroup(rg)
			{
				if (Globals->currYear < RGroup[rg]->startyr)
				{
					/* don't start trying to grow until RGroup[rg]->startyr year */
					continue;
				}

				Int i;

				ForEachEstSpp2( rg, i)
				{
					if (!Species[RGroup[rg]->est_spp[i]]->use_me)
					{
						continue;
					}
					else
					{
					//	printf( "calling Species_Proportion_Recovery() with rgroup name= %s , RGroup[%d]->proportion_recovered =%f for Species[%d]->name= %s \n", RGroup[rg]->name, rg, RGroup[rg]->proportion_recovered, i, Species[RGroup[rg]->est_spp[i]]->name);
						Species_Proportion_Recovery(RGroup[rg]->est_spp[i], 6,
								                    RGroup[rg]->proportion_recovered,
								                    RGroup[rg]->proportion_killed);
					}

				}

			}

		}

	}

}

static void _do_grid_grazing_EndOfYear(int row, int col)
{

	/*======================================================*/
	/* PURPOSE
	 * Perform the sorts of grazing one might expect at end of year, it is based on grazing frequency
	 * HISTORY
	 * Nov 22 2016 -AKT  -Added Species grazing EndOfYear for grid model  
	 */
	/*======================================================*/

	if (UseDisturbances)
	{
		IntU grazingyr = 0;
		int cell = col + ((row - 1) * grid_Cols) - 1;

		if ((Globals->currYear >=grid_Disturb[cell].grazingfreq_startyr) && grid_Disturb[cell].grazing_frq > 0.)
		{
			if (grid_Disturb[cell].grazing_frq < 1.0)
			{
				if (RandUni(&grid_rng) <= grid_Disturb[cell].grazing_frq)
				{
					grazingyr = Globals->currYear;
				}

			}
			else if (((Globals->currYear - grid_Disturb[cell].grazingfreq_startyr) % (IntU) grid_Disturb[cell].grazing_frq) == 0)
			{
				grazingyr = Globals->currYear;
			}

		}

		//rgroup proportion grazing
		if (Globals->currYear == grazingyr)
		{
			GrpIndex rg;
			ForEachGroup(rg)
			{
				if (Globals->currYear < RGroup[rg]->startyr)
				{
					/* don't start trying to grow or do grazing until RGroup[rg]->startyr year */
					continue;
				}

				Int i;

				ForEachEstSpp2( rg, i)
				{
					if (!Species[RGroup[rg]->est_spp[i]]->use_me)
					{
						continue;
					}
					else
					{
						// printf( "year = %d, calling Species_Proportion_Grazing() with rgroup name= %s , RGroup[%d]->proportion_grazing =%f for Species[%d]->name= %s \n", Globals->currYear,RGroup[rg]->name, rg,RGroup[rg]->proportion_grazing, i, Species[RGroup[rg]->est_spp[i]]->name);
			     		Species_Proportion_Grazing(RGroup[rg]->est_spp[i], RGroup[rg]->proportion_grazing);
					}

				}

			}

		}

	}

}

/***********************************************************/
static void _read_init_species(void)
{
	// reads the grid init species input
	// the file should be something like: "cell,copy_cell,copy_which,use_SpinUp,(all the species names seperated by a comma)"
	// there should be no spaces in between, just commas separating the values (ie it should be a .csv file, but does not account for all of the possibilities that a .csv file could be)

	FILE *f;
	int i, j, num, cell, do_copy, copy_cell, use_SpinUp, seeds_Avail,
	    row, col, copy_cell_row, copy_cell_col;
	char buf[4096];

	//open the file/do the reading
	f = OpenFile(grid_files[GRID_FILE_INIT_SPECIES], "r");

	GetALine2(f, buf, 4096); // gets rid of the first line (since it just defines the columns)... it's only there for user readability
	for (i = 0; i < grid_Cells; i++)
	{
	    row = i / grid_Cols;
	    col = i % grid_Cols;

	    load_cell(row, col);

		if (!GetALine2(f, buf, 4096))
			break;

		num = sscanf(buf, "%d,%d,%d,%d", &cell, &do_copy, &copy_cell,
				&use_SpinUp);

		copy_cell_row = copy_cell / grid_Cols;
		copy_cell_col = copy_cell % grid_Cols;

		if (num != 4)
			LogError(logfp, LOGFATAL, "Invalid %s file", grid_files[GRID_FILE_INIT_SPECIES]);

		gridCells[row][col].mySpeciesInit.use_SpinUp = use_SpinUp;

		int stringIndex = _get_value_index(buf, ',', 4); //gets us the index of the string that is right after what we just parsed in

		if (do_copy == 1 && copy_cell > -1 && copy_cell < grid_Cells
				&& cell != 0 && copy_cell < cell)
		{ //copy this cells values from a previous cell's
			for (j = 0; j < Globals->sppCount; j++)
				gridCells[row][col].mySpeciesInit.species_seed_avail[j] =
				    gridCells[copy_cell_row][copy_cell_col].mySpeciesInit.species_seed_avail[j];
			gridCells[row][col].mySpeciesInit.use_SpinUp =
			    gridCells[copy_cell_row][copy_cell_col].mySpeciesInit.use_SpinUp;
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
			if (num != 1)
				LogError(logfp, LOGFATAL,
						"Invalid %s file line %d invalid species input",
						grid_files[GRID_FILE_INIT_SPECIES], i + 2);

			gridCells[row][col].mySpeciesInit.species_seed_avail[s] = seeds_Avail;
			stringIndex += _get_value_index(&buf[stringIndex], ',', 1);
		}
	}

	if (i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file, not enough cells",
				grid_files[GRID_FILE_INIT_SPECIES]);

    unload_cell();
	CloseFile(&f);
}

/***********************************************************/
static void _read_maxrgroupspecies(void)
{
    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    parm_SetName(grid_files[GRID_FILE_MAXRGROUPSPECIES], F_MaxRGroupSpecies);
    maxrgroupspecies_init();
    ChDir("..");
}

/***********************************************************/
static void _read_files(void)
{
    ChDir(grid_directories[GRID_DIRECTORY_STEPWAT_INPUTS]);
    files_init();
    ChDir("..");
}

/***********************************************************/
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

    CloseFile(&f);
}