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
#include "rands.h"

#include "sxw_funcs.h"
#include "sxw.h"
#include "sxw_vars.h"

#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Model.h"
#include "SW_Weather.h"

/***************** Structure Declarations ******************/
/***********************************************************/

struct _grid_soil_lyr_st { // represents a single soil layer
	float data[11];
	int width;
} typedef Grid_Soil_Lyr;

struct _grid_soil_st { //represents the input data for all the soil layers of a cell
	int num_layers;
	Grid_Soil_Lyr* lyr;
} typedef Grid_Soil_St;

struct _grid_disturb_st {
	int choices[3]; //used as boolean values (ie flags as to whether or not to use the specified disturbance)
	int kill_yr, killfrq, extirp;
} typedef Grid_Disturb_St;

struct _grid_sd_struct { //for seed dispersal
	//the idea is that this structure contains all of the cells that the cell represented can possibly disperse seeds to and the probability of the happening for each cell
	int size, seeds_present, seeds_received, *cells; //seeds_present && seeds_received are treated as boolean values... cells & prob are to be the length of size
	float *prob, lyppt;
} typedef Grid_SD_St;

struct _grid_sxw_st { //holds pointers dynamically allocated by SXW.c
	RealD *roots_max, *rootsXphen, *roots_active, *roots_active_rel, *roots_active_sum, *phen, *prod_bmass, *prod_pctlive;
} typedef Grid_SXW_St;

struct _grid_init_species_st {
	int use_SpinUp;
	int *species_seed_avail;
} typedef Grid_Init_Species_St;

/************ Module Variable Declarations ***************/
/***********************************************************/

#define N_GRID_FILES 10
#define N_GRID_DIRECTORIES 1

char *grid_files[N_GRID_FILES], *grid_directories[N_GRID_DIRECTORIES], sd_Sep;

int grid_Cols, grid_Rows, grid_Cells, sd_NYearsSeedsAvailable;
int UseDisturbances, UseSoils, sd_DoOutput, sd_MakeHeader, sd_Option1a, sd_Option1b, sd_Option2a, sd_Option2b; //these are treated like booleans

// these variables are for storing the globals in STEPPE... they are dynamically allocated/freed
SpeciesType	*grid_Species[MAX_SPECIES], *spinup_Species[MAX_SPECIES];
GroupType	*grid_RGroup [MAX_RGROUPS], *spinup_RGroup[MAX_RGROUPS];
SucculentType	*grid_Succulent, *spinup_Succulent;
EnvType		*grid_Env, *spinup_Env;
PlotType	*grid_Plot, *spinup_Plot;
ModelType	*grid_Globals, *spinup_Globals;

// these two variables are for storing SXW variables... also dynamically allocated/freed
SXW_t *grid_SXW, *spinup_SXW;
Grid_SXW_St *grid_SXW_ptrs, *spinup_SXW_ptrs;

// these are SOILWAT variables that we need...
extern SW_SOILWAT SW_Soilwat;
extern SW_SITE SW_Site;
extern SW_VEGPROD SW_VegProd;
extern SW_WEATHER SW_Weather;

// these are grids to store the SOILWAT variables... also dynamically allocated/freed
SW_SOILWAT *grid_SW_Soilwat, *spinup_SW_Soilwat;
SW_SITE *grid_SW_Site, *spinup_SW_Site;
SW_VEGPROD *grid_SW_VegProd, *spinup_SW_VegProd;
SW_MODEL *grid_SW_Model, *spinup_SW_Model;

// these two variables are used to store the soil/distubance inputs for each grid cell... also dynamically allocated/freed
Grid_Soil_St *grid_Soils;
Grid_Disturb_St *grid_Disturb;
Grid_Init_Species_St *grid_initSpecies;

Grid_SD_St *grid_SD[MAX_SPECIES]; //for seed dispersal

// these variables are used for the soil types in the spinup options
int nSoilTypes, *soilTypes_Array, *grid_SoilTypes;

// these are both declared and set in the ST_main.c module
extern Bool UseSoilwat;
extern Bool UseProgressBar;

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/

void runGrid( void ); //to be called from ST_main.c

/************* External Function Declarations **************/
/***********************************************************/

//from ST_resgroups.c
void rgroup_Grow( void);
void rgroup_Establish( void) ;
void rgroup_IncrAges( void);
void rgroup_PartResources( void);
void rgroup_ResPartIndiv(void);

//from ST_mortality.c
void mort_Main( Bool *killed);
void mort_EndOfYear( void);

//functions from ST_params.c
void parm_Initialize( Int);
void parm_SetFirstName( char *s);
void parm_SetName( char *s, int which);
void parm_free_memory( void );

//from ST_main.c
void Plot_Initialize( void );

//functions from ST_stats.c
void stat_Collect(Int year);
void stat_Collect_GMort( void );
void stat_Collect_SMort( void );
void stat_Output_AllMorts( void) ;
void stat_Output_AllBmass(void) ;
void stat_Output_Seed_Dispersal(const char * filename, const char sep, Bool makeHeader); 
void stat_Load_Accumulators(int cell, int year);
void stat_Save_Accumulators(int cell, int year);
void stat_Free_Accumulators( void );
void stat_Init_Accumulators( void );

//functions from sxw.c
void free_sxw_memory( void ); 
void free_all_sxw_memory( void );
void load_sxw_memory( RealD * grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive );
void save_sxw_memory( RealD * grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive );
void SXW_init( Bool init_SW );

//from SW_Site.c (both needed to initialize the soil layers properly)
void water_eqn(RealD fractionGravel, RealD sand, RealD clay, LyrIndex n);
void init_site_info(void); 

/*********** Locally Used Function Declarations ************/
/***********************************************************/

static int _load_bar(char* prefix, clock_t start, int x, int n, int r, int w);
static double _time_remaining(clock_t start, char* timeChar, double percentDone); 
static void _run_spinup( void );
static void _init_grid_files( void );
static void _init_grid_inputs( void );
static void _init_SXW_inputs( Bool init_SW );
static void _init_stepwat_inputs( void );
static void _init_grid_globals( void );
static void _init_spinup_globals( void );
static void _load_grid_globals( void );
static void _load_spinup_globals( void );
static void _free_grid_memory( void );
static void _free_spinup_memory( void );
static void _free_grid_globals( void );
static void _free_spinup_globals( void );
static void _load_cell( int row, int col, int year, Bool useAccumulators );
static void _load_spinup_cell( int cell );
static void _save_cell( int row, int col, int year, Bool useAccumulators );
static void _save_spinup_cell(int cell );
static void _read_disturbances_in( void );
static void _read_soils_in( void );
static void _init_soil_layers(int cell, int isSpinup);
static float _read_a_float(FILE *f, char *buf, const char *filename, const char *descriptor);
static float _cell_dist(int row1, int row2, int col1, int col2, float cellLen);
static void _read_seed_dispersal_in( void );
static void _do_seed_dispersal( void );
static void _set_sd_lyppt(int row, int col);
static void _kill_groups_and_species( void );
static int  _do_grid_disturbances(int row, int col);
static void _read_init_species( void );

//static void _copy_species(SpeciesType* to, SpeciesType* from, Bool deep);

static IndivType* _create_empty_indv( void ); //these 5 functions are used for copying/freeing the linked list of individuals correctly...
static void _free_individuals( IndivType *head );
static void _free_head( IndivType *head );
static IndivType* _copy_individuals( IndivType *head );
static IndivType* _copy_head( IndivType *head );

/******************** Begin Model Code *********************/
/***********************************************************/

/***********************************************************/
double _time_remaining(clock_t start, char* timeChar, double percentDone) {
	// for use in _load_bar() function... pretty much it returns the amount of time left and a character in timeChar representing the units
	// percent done must be from 0 to 1, with 0 representing no work done (0%) and 1 representing all of the work done (ie 100%)
	clock_t timeElapsed = clock() - start; //gets the time since we started (this is in CPU cycles).  We convert to seconds in the next part
	double result = (((double)timeElapsed)/CLOCKS_PER_SEC) / percentDone * (1.0 - percentDone); //gets estimated time left until the work is complete based upon timeElapsed and percentDone 

	*timeChar = 's'; //get the biggest time units that we can use (ie no use displaying 120 seconds if we can display 2 minutes)
	if(result > 60) {
		result /= 60.0;
		*timeChar = 'm';
		if(result > 60) {
			result /= 60.0;
			*timeChar = 'h';
			if(result > 24) {
				result /= 24.0;
				*timeChar = 'd';
				if(result > 7) {
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
	if ( x % (n/r) != 0 ) return (int) ((x/(float)n)*100);	
	
    	int i;    

	printf("\r"); //we output a carriage-return character to put us back at the beginning of the line	
    	if(prefix != NULL) printf("%s", prefix);
 
    	// Calculuate the ratio of complete-to-incomplete.
    	float ratio = x/(float)n;
    	int   c     = ratio * w;
	int result = (int)(ratio*100);
	
	if(result > 1) { //we don't give an estimate if less than 2% of the work is complete because we don't have any data to go off of...
		char timeChar;
		double timeLeft = _time_remaining(start, &timeChar, (double) ratio);
		if(timeLeft < 10)
			printf("(est 0%.2f%c) ", timeLeft, timeChar);
		else
			printf("(est %.2f%c) ", timeLeft, timeChar);
	}

	// Show the percentage complete.
    	printf("%3d%% [", (int)(ratio*100) );
 
    	// Show the load bar.
    	for (i=0; i<c; i++)
       		printf("=");
 
    	for (i=c; i<w; i++)
       		printf(" ");
 
    	printf("]");
	fflush(stdout); //we do this to flush (ie. print) the output, since stdout typically waits for a newline character before flushing 
	return result;
}

/***********************************************************/
void runGrid( void ) {
	// this function sets up & runs the grid
	
	_init_grid_files();				// reads in files.in file
	_init_stepwat_inputs();				// reads the stepwat inputs in
	_init_grid_inputs();				// reads the grid inputs in & initializes the global grid variables
	
	if(sd_Option2a || sd_Option2b)
		_run_spinup();				// does the initial spinup

	double prog_Percent = 0.0, prog_Incr, prog_Acc = 0.0;
	char prog_Prefix[32];
	clock_t prog_Time;
	int i, j;
	Bool killedany;
	IntS year, iter;
	if(UseProgressBar) {
		prog_Incr = (((double)1)/ ((double)((Globals.runModelYears*grid_Cells)*Globals.runModelIterations)));  //gets how much progress we'll make in one year towards our goal of iter*years*cells	
		prog_Time = clock();  //used for timing
		sprintf(prog_Prefix, "simulations: ");
	}

	for(iter = 1; iter <= Globals.runModelIterations; iter++) { //for each iteration
	
		if (BmassFlags.yearly || MortFlags.yearly)
        		parm_Initialize( iter);
        	
		Plot_Initialize();
		if(iter > 1) _free_grid_globals(); //frees the memory from when we called _load_grid_globals() last time... (doesn't need to be called on the first iteration because the memory hasn't been allocated yet)
		
		Globals.currIter = iter;
		_load_grid_globals(); //allocates/initializes grid variables (specifically the ones that are going to change every iter)
		
		for( year=1; year <= Globals.runModelYears; year++) {//for each year
			for(i = 1; i <= grid_Rows; i++)
				for(j = 1; j <= grid_Cols; j++) { //for each cell

					//fprintf(stderr, "year: %d", year);
					_load_cell(i, j, year, TRUE);
					
					if(year == 1 && (sd_Option2a || sd_Option2b)) {
						// finish this...
							
						int cell = j + ( (i-1) * grid_Cols) - 1;  // converts the row/col into an array index
						int spinup_cell = grid_SoilTypes[cell]; // gets us the actual cell we want to load up from the spinup

						if( (sd_Option2b && grid_initSpecies[cell].use_SpinUp) || sd_Option2a )
							_load_spinup_cell(spinup_cell); // loads the spinup cell into the global variables
					}

	          			Globals.currYear = year;
				
					if(year > 1 && UseSeedDispersal)
						_set_sd_lyppt(i, j);	

					_do_grid_disturbances(i, j);
					
					rgroup_Establish();  /* excludes annuals */

          				Env_Generate(); //if UseSoilwat it calls : SXW_Run_SOILWAT() which calls : _sxw_sw_run() which calls : SW_CTL_run_current_year()
          				
					rgroup_PartResources();
					rgroup_Grow();
					
					mort_Main( &killedany);
					
					rgroup_IncrAges();
						
					
         				stat_Collect(year);
					mort_EndOfYear();
         				
         				_save_cell(i, j, year, TRUE);
         		
         				if(UseProgressBar) {
         					prog_Percent += prog_Incr; //updating our percent done
         				if(prog_Percent > prog_Acc) { //only update if 1% progress or more has been made since the last time we updated (this check is so it doesn't waste processing time that could be spent running the simulations by updating all of the time)
         					prog_Acc += 0.01;
         					_load_bar(prog_Prefix, prog_Time, (int) (100 * prog_Percent), 100, 100, 10); //display our bar to the console
         				}
         			}
    			} /* end model run for this cell*/
    			
			if(UseSeedDispersal)
				_do_seed_dispersal();
		}/* end model run for this year*/	
    	
		// collects the data appropriately for the mort output... (ie. fills the accumulators in ST_stats.c with the values that they need)
		if(MortFlags.summary)
			for( i = 1; i <= grid_Rows; i++)
				for( j = 1; j <= grid_Cols; j++) {
					_load_cell(i, j, Globals.runModelYears, TRUE);
        				stat_Collect_GMort();
        				stat_Collect_SMort();
   					_save_cell(i, j, Globals.runModelYears, TRUE);
   				}
   					
	} /*end iterations */

    	if(UseProgressBar) printf("\rsimulations took approximately: %.2f seconds\n", ((double)(clock() - prog_Time) / CLOCKS_PER_SEC));
    
	if(UseProgressBar) {
 		prog_Percent = prog_Acc = 0.0;
 		prog_Incr = ((double)1)/ ((double)grid_Cells);
 		prog_Time = clock();
 		sprintf(prog_Prefix, "outputting: ");
	}	
 
 	// outputs all of the mort and BMass files for each cell...
	for(i = 1; i <= grid_Rows; i++)
  		for(j = 1; j <= grid_Cols; j++) {	
  	
  			int cell = j + ( (i-1) * grid_Cols) - 1;
  			_load_cell(i, j, 1, TRUE);
  			for( year=2; year <= Globals.runModelYears; year++) // _load_cell gets the first years accumulators loaded, so we start at 2...
  				stat_Load_Accumulators(cell, year);
  				
  			char fileMort[1024], fileBMass[1024], fileReceivedProb[1024];
  		
			sprintf(fileReceivedProb, "%s%d.out", grid_files[9], cell);
  			sprintf(fileMort, "%s%d.out", grid_files[8], cell);
  			sprintf(fileBMass, "%s%d.out", grid_files[7], cell);
  			parm_SetName(fileMort, F_MortAvg);
  			parm_SetName(fileBMass, F_BMassAvg);
  		
  			if (MortFlags.summary)
    				stat_Output_AllMorts();
  			if (BmassFlags.summary)
    				stat_Output_AllBmass();
			if (UseSeedDispersal && sd_DoOutput)
				stat_Output_Seed_Dispersal(fileReceivedProb, sd_Sep, sd_MakeHeader); 

        		if(UseProgressBar) {
        			prog_Percent += prog_Incr;
        			if(prog_Percent > prog_Acc) {
        				prog_Acc += 0.01;
         				_load_bar(prog_Prefix, prog_Time, 100 * prog_Percent, 100, 100, 10);
         			}
         		}
    		
    	}
	if(UseProgressBar) printf("\routputting files took approximately %.2f seconds\n", ((double)(clock() - prog_Time) / CLOCKS_PER_SEC));
	_free_grid_memory(); // free our allocated memory since we do not need it anymore
	/*if(UseProgressBar)*/ printf("!\n");
}


/***********************************************************/
static void _run_spinup( void ) {

	//does the spinup, it's pretty much like running the grid, except some differences like no seed dispersal and no need to deal with output accumulators

	int i, spinup_Cell;
	Bool killedany;
	IntS year, iter;

	DuringSpinup = TRUE;

	if(!UseSoils || !UseSoilwat) { // if we're not using inputting soils then there is simply one soil type as all the soils are the same
		nSoilTypes = 1;
		soilTypes_Array[0] = 0;
		for(i = 0; i < grid_Cells; i++)
			grid_SoilTypes[i] = 0;
	}

	_init_spinup_globals();
	_load_spinup_globals();

	for(iter = 1; iter <= 1; iter++) { //for each iteration... only 1 iteration allowed for now
	
		if (BmassFlags.yearly || MortFlags.yearly)
        		parm_Initialize( iter);
        	
		Plot_Initialize();
		
		Globals.currIter = iter;
		//_load_grid_globals(); //allocates/initializes grid variables (specifically the ones that are going to change every iter)

		for( year=1; year <= Globals.runModelYears; year++) { //for each year
			for(spinup_Cell = 0; spinup_Cell < nSoilTypes; spinup_Cell++) { // for each different soil type

					int cell = soilTypes_Array[spinup_Cell]; // this is the first cell of this soiltypes actual cell number
					int j = cell % grid_Cols + 1; // this is the column of the first cell of this soiltype
					int i = ((cell + 1 - j) / grid_Cols) + 1; // this is the row of the first cell of this soiltype 

					_load_spinup_cell(spinup_Cell);
	          			Globals.currYear = year;

					_do_grid_disturbances(i, j);
					
					rgroup_Establish();  /* excludes annuals */

          				Env_Generate(); //if UseSoilwat it calls : SXW_Run_SOILWAT() which calls : _sxw_sw_run() which calls : SW_CTL_run_current_year()
          				
					rgroup_PartResources();
					rgroup_Grow();
					
					mort_Main( &killedany);
					
					rgroup_IncrAges();
						
					
         				stat_Collect(year);
					mort_EndOfYear();

					_save_spinup_cell(spinup_Cell);
         				
    			} /* end model run for this cell*/
    			
		}/* end model run for this year*/	
    		
		//_free_grid_globals(); //free's the grid variables that change every iter
	} /*end iterations */

	DuringSpinup = FALSE;
}

/***********************************************************/
static void _init_grid_files( void ) {
	// reads the files.in file

	FILE *f;
	char buf[1024];
	int i;
	
	f = OpenFile(Parm_name(F_First), "r");
    
    	for(i = 0; i < N_GRID_DIRECTORIES; i++) { //0 is stepwat directory
    		if(!GetALine(f, buf)) break;
    		grid_directories[i] = Str_Dup(Str_TrimLeftQ(buf));
    	}
    	if(i != N_GRID_DIRECTORIES) LogError(stderr, LOGFATAL, "Invalid files.in");
    
    	for(i = 0; i < N_GRID_FILES; i++) {
    		if(!GetALine(f, buf)) break;
    		grid_files[i] = Str_Dup(Str_TrimLeftQ(buf));
    	}
    	if(i !=  N_GRID_FILES) LogError(stderr, LOGFATAL, "Invalid files.in");
    
    	// opens the log file...
    	if ( !strcmp("stdout", grid_files[0]) )
    		logfp = stdout;
  	 else
    		logfp = OpenFile(grid_files[0], "w"); //grid_files[0] is the logfile to use
    	
	 /*printf("stepwat dir: %s\n", grid_directories[0]);
    	for(i = 0; i < N_GRID_FILES; i++)
    		printf("%d : %s\n", i, grid_files[i]);*/
    
    	CloseFile(&f);
}

/***********************************************************/
static void _init_grid_inputs( void ) {
	// reads in the grid input file
	
	FILE *f;
	char buf[1024];
	int i, j;

	f = OpenFile(grid_files[1], "r"); //grid_files[1] is the grid inputs file
	
	GetALine(f, buf);
	i=sscanf( buf, "%d %d", &grid_Rows, &grid_Cols);
	if(i != 2)
		LogError(logfp, LOGFATAL, "Invalid grid setup file (rows/cols line wrong)");
	
	grid_Cells = grid_Cols * grid_Rows;
	
	if(grid_Cells > MAX_CELLS)
		LogError(logfp, LOGFATAL, "Number of cells in grid exceeds MAX_CELLS defined in ST_defines.h");
		
	Globals.nCells = (grid_Cols * grid_Rows);
	
	GetALine(f, buf);
	i=sscanf( buf, "%d", &UseDisturbances );
	if(i != 1)
		LogError(logfp, LOGFATAL, "Invalid grid setup file (disturbances line wrong)");
		
	GetALine(f, buf);
	i=sscanf( buf, "%d", &UseSoils );
	if(i != 1)
		LogError(logfp, LOGFATAL, "Invalid grid setup file (soils line wrong)");
		
	GetALine(f, buf);
	i=sscanf( buf, "%d", &j );
	if(i != 1)
		LogError(logfp, LOGFATAL, "Invalid grid setup file (seed dispersal line wrong)");
	UseSeedDispersal = itob(j);	

	CloseFile(&f);
	
	_init_grid_globals(); // initializes the global grid variables
	if(UseDisturbances)	
		_read_disturbances_in();
	if(UseSeedDispersal) {
		_read_seed_dispersal_in();
		_read_init_species();
	}
	if(UseSoils && UseSoilwat)
		_read_soils_in();

	DuringSpinup = FALSE;
}

/***********************************************************/
static void _init_SXW_inputs( Bool init_SW ) {
    	SXW_Init(init_SW);	// initializes soilwat
    	if(init_SW == TRUE) {
    		char aString[2048];
    		sprintf(aString, "%s/%s", grid_directories[0], SW_Weather.name_prefix);
    		sprintf(SW_Weather.name_prefix, "%s", aString);  //updates the directory correctly for the weather files so soilwat can find them
    	}
}

/***********************************************************/
static void _init_stepwat_inputs( void ) {
	// reads in the stepwat inputs
	ChDir(grid_directories[0]);			// changes to the folder that the stepwat input is in
	
	parm_SetFirstName(grid_files[6]);		// correctly sets the name of the stepwat files.in file
	parm_Initialize( 0);				// loads stepwat input files
	
	if(UseSoilwat)
		_init_SXW_inputs(TRUE);
	
	ChDir("..");					// goes back to the folder that we were in
}

/***********************************************************/
static IndivType* _create_empty_indv( void ) { 
  	//simply allocates memory for an individual and returns the pointer to it

  	return Mem_Calloc(1, sizeof(IndivType), "_create_empty_indv()");

}

/***********************************************************/
static void _free_individuals( IndivType *head ) {
	//frees the memory allocated in the linked list of individuals pointed to by head
	if(head->Next != NULL)
		_free_individuals(head->Next); //recursively call itself in order to free the next obj in the list (Mem_Free() will end up getting called on the objects in the list from last to first)
	Mem_Free(head);
}

/***********************************************************/
static void _free_head( IndivType * head ) {
	if(head == NULL) return; //if head == NULL then there is no memory to free, so just return
	_free_individuals(head);
}

/***********************************************************/
static IndivType* _copy_individuals( IndivType *head ) { 
	//performs a deep copy (not to be confused with a shallow copy) since that's what we want
	IndivType *n = _create_empty_indv(); //allocate space for the individual to be copied to
	
	*n = *head;

	if(head->Next != NULL) {
		n->Next = _copy_individuals(head->Next); //recursively call itself in order to copy the next obj in the list
		n->Next->Prev = n;
	} else
		n->Next = NULL;
	
	return n;
}

/***********************************************************/
static IndivType* _copy_head( IndivType *head ) {
	if(head == NULL) return NULL; //if head == NULL then there is nothing to copy...
	return _copy_individuals(head);
}

/**********************************************************/
/*
static void _copy_species(SpeciesType* to, SpeciesType* from, Bool deep) {
	if(deep) {
		Mem_Free(to->kills);
		Mem_Free(to->seedprod);
		_free_head(to->IndvHead); //free_head() frees the memory allocated by the head and the memory allocated by each part of the linked list
	}

	//*to = *from;
	to->est_count = from->est_count;
	to->estabs = from->estabs;
	to->relsize = from->relsize;
	to->extragrowth = from->extragrowth;
	to->allow_growth = from->allow_growth;

	if(deep) {
		to->kills = Mem_Calloc(from->max_age, sizeof(IntUS), "_load_cell(Species[s]->kills)");
		to->seedprod = Mem_Calloc(from->viable_yrs, sizeof(RealF), "_load_cell(Species[s]->seedprod)");

		memcpy(to->kills, from->kills, from->max_age * sizeof(IntUS));
		memcpy(to->seedprod, from->seedprod, from->viable_yrs * sizeof(RealF));
		to->IndvHead = _copy_head(from->IndvHead); //copy_head() deep copies the linked list structure (allocating memory when needed)... it will even allocate memory for the head of the list
	}
}
*/
/***********************************************************/
static void _init_grid_globals( void ) {
	//initializes grid variables, allocating the memory necessary for them (this step is only needed to be done once)
	
	GrpIndex c;
	SppIndex s;
	int i;
	
	grid_Succulent = Mem_Calloc(grid_Cells, sizeof(SucculentType), "_init_grid_globals()");
	grid_Env = Mem_Calloc(grid_Cells, sizeof(EnvType), "_init_grid_globals()");
	grid_Plot = Mem_Calloc(grid_Cells, sizeof(PlotType), "_init_grid_globals()");
	grid_Globals = Mem_Calloc(grid_Cells, sizeof(ModelType), "_init_grid_globals()");
	
	ForEachSpecies(s)
		if(Species[s]->use_me) grid_Species[s] = Mem_Calloc(grid_Cells, sizeof(SpeciesType), "_init_grid_globals()");
	ForEachGroup(c)
		if(RGroup[c]->use_me) grid_RGroup[c] = Mem_Calloc(grid_Cells, sizeof(GroupType), "_init_grid_globals()");
	
	if(UseSoilwat) {
		grid_SXW = Mem_Calloc(grid_Cells, sizeof(SXW_t), "_init_grid_globals()");
		grid_SW_Soilwat = Mem_Calloc(grid_Cells, sizeof(SW_SOILWAT), "_init_grid_globals()");
		grid_SW_Site = Mem_Calloc(grid_Cells, sizeof(SW_SITE), "_init_grid_globals()");
		grid_SW_VegProd = Mem_Calloc(grid_Cells, sizeof(SW_VEGPROD), "_init_grid_globals()");
		if(UseSoils) {
			grid_Soils = Mem_Calloc(grid_Cells, sizeof(Grid_Soil_St), "_init_grid_globals()");
			for(i = 0; i < grid_Cells; i++)
				grid_Soils[i].num_layers = 0;
			grid_SXW_ptrs = Mem_Calloc(grid_Cells, sizeof(Grid_SXW_St), "_init_grid_globals()");
		}
	}
	
	if(UseDisturbances)
		grid_Disturb = Mem_Calloc(grid_Cells, sizeof(Grid_Disturb_St), "_init_grid_globals()");
	if(UseSeedDispersal) {
		ForEachSpecies(s)
			if(Species[s]->use_me && Species[s]->use_dispersal) grid_SD[s] = Mem_Calloc(grid_Cells, sizeof(Grid_SD_St), "_init_grid_globals()");

		grid_initSpecies = Mem_Calloc(grid_Cells, sizeof(Grid_Init_Species_St), "_init_grid_globals()");
		for(i = 0; i < grid_Cells; i++)
			grid_initSpecies[i].species_seed_avail = Mem_Calloc(Globals.sppCount, sizeof(int), "_init_grid_globals()");
	}


	//if(sd_Option2a || sd_Option2b) {
		soilTypes_Array = Mem_Calloc(grid_Cells, sizeof(int*), "_init_grid_globals()");
		grid_SoilTypes = Mem_Calloc(grid_Cells, sizeof(int*), "_init_grid_globals()");
	//}
	
	stat_Init_Accumulators();
}

/***********************************************************/
static void _init_spinup_globals( void ) {
	//initializes spinup variables, allocating the memory necessary for them (this step is only needed to be done once)
	
	GrpIndex c;
	SppIndex s;
	
	spinup_Succulent = Mem_Calloc(nSoilTypes, sizeof(SucculentType), "_init_spinup_globals()");
	spinup_Env = Mem_Calloc(nSoilTypes, sizeof(EnvType), "_init_spinup_globals()");
	spinup_Plot = Mem_Calloc(nSoilTypes, sizeof(PlotType), "_init_spinup_globals()");
	spinup_Globals = Mem_Calloc(nSoilTypes, sizeof(ModelType), "_init_spinup_globals()");
	
	ForEachSpecies(s)
		if(Species[s]->use_me) spinup_Species[s] = Mem_Calloc(nSoilTypes, sizeof(SpeciesType), "_init_spinup_globals()");
	ForEachGroup(c)
		if(RGroup[c]->use_me) spinup_RGroup[c] = Mem_Calloc(nSoilTypes, sizeof(GroupType), "_init_spinup_globals()");
	
	if(UseSoilwat) {
		spinup_SXW = Mem_Calloc(nSoilTypes, sizeof(SXW_t), "_init_spinup_globals()");
		spinup_SW_Soilwat = Mem_Calloc(nSoilTypes, sizeof(SW_SOILWAT), "_init_spinup_globals()");
		spinup_SW_Site = Mem_Calloc(nSoilTypes, sizeof(SW_SITE), "_init_spinup_globals()");
		spinup_SW_VegProd = Mem_Calloc(nSoilTypes, sizeof(SW_VEGPROD), "_init_spinup_globals()");
		if(UseSoils) {
			spinup_SXW_ptrs = Mem_Calloc(nSoilTypes, sizeof(Grid_SXW_St), "_init_spinup_globals()");
		}
	}
}

/***********************************************************/
static void _load_grid_globals( void ) {
	//this initializes/allocates memory needed... this step is needed to be done for every iteration

	int i, j;
	GrpIndex c;
	SppIndex s;
	
	if(UseSoils && UseSoilwat) ChDir(grid_directories[0]); //change the directory for _init_soil_layers()
	for(i = 0; i < grid_Cells; i++) {
		
		ForEachSpecies(s) { //macros defined in ST_defines.h
			if(!Species[s]->use_me) continue;
			grid_Species[s][i] = *Species[s];

			grid_Species[s][i].kills = Mem_Calloc(Species[s]->max_age, sizeof(IntUS), "_init_grid_globals()");
			grid_Species[s][i].seedprod = Mem_Calloc(Species[s]->viable_yrs, sizeof(RealF), "_init_grid_globals()");
			
			memcpy(grid_Species[s][i].kills, Species[s]->kills, Species[s]->max_age * sizeof(IntUS));
			memcpy(grid_Species[s][i].seedprod, Species[s]->seedprod, Species[s]->viable_yrs * sizeof(RealF));
			
			grid_Species[s][i].IndvHead = _copy_head(Species[s]->IndvHead); //copy_head() deep copies the structure (allocating memory when needed)... it will even allocate memory for the head of the list
		}
		
		ForEachGroup(c) {
			if(!RGroup[c]->use_me) continue;
			grid_RGroup [c][i] = *RGroup[c];
			grid_RGroup [c][i].kills = Mem_Calloc(RGroup[c]->max_age, sizeof(IntUS), "_init_grid_globals()");
			
			memcpy(grid_RGroup[c][i].kills, RGroup[c]->kills, RGroup[c]->max_age * sizeof(IntUS));
			if(UseDisturbances) { 
				grid_RGroup[c][i].killyr = grid_Disturb[i].kill_yr;
				grid_RGroup[c][i].killfreq = grid_Disturb[i].killfrq;
				grid_RGroup[c][i].extirp = grid_Disturb[i].extirp;
			}
		}
			
		grid_Succulent[i] = Succulent;
		grid_Env[i] = Env;
		grid_Plot[i] = Plot;
		grid_Globals[i] = Globals;
		
		if(UseDisturbances) {
			grid_Globals[i].pat.use = grid_Disturb[i].choices[0];
			grid_Globals[i].mound.use = grid_Disturb[i].choices[1];
			grid_Globals[i].burrow.use = grid_Disturb[i].choices[2];
		}
		if(UseSoils && UseSoilwat)
			_init_soil_layers(i, 0);
				
		if(UseSoilwat) {
			grid_SXW[i] = SXW; 
		
			grid_SXW[i].transpTotal = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealD), "_init_grid_globals()");
			grid_SXW[i].swc = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealF), "_init_grid_globals()");
			
			memcpy(grid_SXW[i].transpTotal, SXW.transpTotal, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
			memcpy(grid_SXW[i].swc, SXW.swc, SXW.NPds * SXW.NSoLyrs * sizeof(RealF));
			
			grid_SW_Soilwat[i] = SW_Soilwat;
			grid_SW_Site[i] = SW_Site;
			grid_SW_VegProd[i] = SW_VegProd;
			
        		grid_SW_Site[i].lyr = Mem_Calloc(SW_Site.n_layers, sizeof(SW_LAYER_INFO *), "_init_grid_globals()");
        		for(j = 0; j < SW_Site.n_layers; j++) {
        			grid_SW_Site[i].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_init_grid_globals()");
        			memcpy(grid_SW_Site[i].lyr[j], SW_Site.lyr[j], sizeof(SW_LAYER_INFO));
        		}
		}
	}
	if(UseSoils && UseSoilwat) ChDir(".."); //get back to our previous directory
	
}

/***********************************************************/
static void _load_spinup_globals( void ) {
	//this initializes/allocates memory needed... this step is needed to be done for every iteration

	int i, j;
	GrpIndex c;
	SppIndex s;
	
	if(UseSoils && UseSoilwat) ChDir(grid_directories[0]); //change the directory for _init_soil_layers()
	for(i = 0; i < nSoilTypes; i++) {
	
		int cell = soilTypes_Array[i]; //this is the cell number of the first cell representing this soil type

		ForEachSpecies(s) { //macros defined in ST_defines.h
			if(!Species[s]->use_me) continue;
			spinup_Species[s][i] = *Species[s];

			spinup_Species[s][i].kills = Mem_Calloc(Species[s]->max_age, sizeof(IntUS), "_init_spinup_globals()");
			spinup_Species[s][i].seedprod = Mem_Calloc(Species[s]->viable_yrs, sizeof(RealF), "_init_spinup_globals()");
			
			memcpy(spinup_Species[s][i].kills, Species[s]->kills, Species[s]->max_age * sizeof(IntUS));
			memcpy(spinup_Species[s][i].seedprod, Species[s]->seedprod, Species[s]->viable_yrs * sizeof(RealF));
			
			spinup_Species[s][i].IndvHead = _copy_head(Species[s]->IndvHead); //copy_head() deep copies the structure (allocating memory when needed)... it will even allocate memory for the head of the list
		}
		
		ForEachGroup(c) {
			if(!RGroup[c]->use_me) continue;
			spinup_RGroup [c][i] = *RGroup[c];
			spinup_RGroup [c][i].kills = Mem_Calloc(RGroup[c]->max_age, sizeof(IntUS), "_init_spinup_globals()");
			
			memcpy(spinup_RGroup[c][i].kills, RGroup[c]->kills, RGroup[c]->max_age * sizeof(IntUS));
			if(UseDisturbances) { 
				spinup_RGroup[c][i].killyr = grid_Disturb[cell].kill_yr;
				spinup_RGroup[c][i].killfreq = grid_Disturb[cell].killfrq;
				spinup_RGroup[c][i].extirp = grid_Disturb[cell].extirp;
			}
		}
			
		spinup_Succulent[i] = Succulent;
		spinup_Env[i] = Env;
		spinup_Plot[i] = Plot;
		spinup_Globals[i] = Globals;
		
		if(UseDisturbances) {
			spinup_Globals[i].pat.use = grid_Disturb[cell].choices[0];
			spinup_Globals[i].mound.use = grid_Disturb[cell].choices[1];
			spinup_Globals[i].burrow.use = grid_Disturb[cell].choices[2];
		}
		if(UseSoils && UseSoilwat)//TODO: should i be cell value after lookup from soilTypes_Array
			_init_soil_layers(i, 1);
				
		if(UseSoilwat) {
			spinup_SXW[i] = SXW; 
			
			spinup_SXW[i].transpTotal = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealD), "_init_spinup_globals()");
			spinup_SXW[i].swc = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealF), "_init_spinup_globals()");
			
			memcpy(spinup_SXW[i].transpTotal, SXW.transpTotal, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
			memcpy(spinup_SXW[i].swc, SXW.swc, SXW.NPds * SXW.NSoLyrs * sizeof(RealF));
			
			spinup_SW_Soilwat[i] = SW_Soilwat;
			spinup_SW_Site[i] = SW_Site;
			spinup_SW_VegProd[i] = SW_VegProd;
			
        		spinup_SW_Site[i].lyr = Mem_Calloc(SW_Site.n_layers, sizeof(SW_LAYER_INFO *), "_init_grid_globals()");
        		for(j = 0; j < SW_Site.n_layers; j++) {
        			spinup_SW_Site[i].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_init_grid_globals()");
        			memcpy(spinup_SW_Site[i].lyr[j], SW_Site.lyr[j], sizeof(SW_LAYER_INFO));
        		}
		}
	}
	if(UseSoils && UseSoilwat) ChDir(".."); //get back to our previous directory
	
}

/***********************************************************/
static void _free_grid_globals( void ) {
	//frees memory allocated in _load_grid_globals() function.
	int i, j;
	GrpIndex c;
	SppIndex s;
	
	for( i = 0; i < grid_Cells; i++ ) {
	
		ForEachSpecies(s) {
			if(!Species[s]->use_me) continue;
			Mem_Free(grid_Species[s][i].kills);
			Mem_Free(grid_Species[s][i].seedprod);
			_free_head(grid_Species[s][i].IndvHead);
		}
		
		ForEachGroup(c)
			if(RGroup[c]->use_me) Mem_Free(grid_RGroup[c][i].kills);
			
		if(UseSoilwat) {
			Mem_Free(grid_SXW[i].transpTotal);
			Mem_Free(grid_SXW[i].swc);
			if(UseSoils) {
				Mem_Free(grid_SXW_ptrs[i].roots_max);
   				Mem_Free(grid_SXW_ptrs[i].rootsXphen);
				Mem_Free(grid_SXW_ptrs[i].roots_active);
				Mem_Free(grid_SXW_ptrs[i].roots_active_rel);
				Mem_Free(grid_SXW_ptrs[i].roots_active_sum);
				Mem_Free(grid_SXW_ptrs[i].phen);
			}
			
			for(j = 0; j < grid_SW_Site[i].n_layers; j++)
				Mem_Free(grid_SW_Site[i].lyr[j]);
			Mem_Free(grid_SW_Site[i].lyr);
		}
			
	}
	
}

/***********************************************************/
static void _free_spinup_globals( void ) {
	//frees memory allocated in _load_spinup_globals() function.
	int i, j;
	GrpIndex c;
	SppIndex s;
	
	for( i = 0; i < nSoilTypes; i++ ) {
	
		ForEachSpecies(s) {
			if(!Species[s]->use_me) continue;
			Mem_Free(spinup_Species[s][i].kills);
			Mem_Free(spinup_Species[s][i].seedprod);
			_free_head(spinup_Species[s][i].IndvHead);
		}
		
		ForEachGroup(c)
			if(RGroup[c]->use_me) Mem_Free(spinup_RGroup[c][i].kills);
			
		if(UseSoilwat) {
			Mem_Free(spinup_SXW[i].transpTotal);
			Mem_Free(spinup_SXW[i].swc);
			if(UseSoils) {
				Mem_Free(spinup_SXW_ptrs[i].roots_max);
   				Mem_Free(spinup_SXW_ptrs[i].rootsXphen);
				Mem_Free(spinup_SXW_ptrs[i].roots_active);
				Mem_Free(spinup_SXW_ptrs[i].roots_active_rel);
				Mem_Free(spinup_SXW_ptrs[i].roots_active_sum);
				Mem_Free(spinup_SXW_ptrs[i].phen);
			}
			
			for(j = 0; j < grid_SW_Site[i].n_layers; j++)
				Mem_Free(spinup_SW_Site[i].lyr[j]);
			Mem_Free(spinup_SW_Site[i].lyr);
		}
			
	}
	
}

/***********************************************************/
static void _free_grid_memory( void ) { 
	//frees all the memory allocated in this file ST_Grid.c (most of it is dynamically allocated in _init_grid_globals() & _load_grid_globals() functions)
	
	int i;
	GrpIndex c;
	SppIndex s;
	
	_free_grid_globals();
	if(sd_Option2a || sd_Option2b)
		_free_spinup_memory();
	
	ForEachSpecies(s)
		if(Species[s]->use_me) Mem_Free(grid_Species[s]);
	ForEachGroup(c)
		if(RGroup[c]->use_me) Mem_Free(grid_RGroup[c]);
	
	Mem_Free(grid_Succulent);
	Mem_Free(grid_Env);
	Mem_Free(grid_Plot);
	Mem_Free(grid_Globals);
	if(UseSoilwat) {
		Mem_Free(grid_SXW);
		Mem_Free(grid_SW_Soilwat);
		Mem_Free(grid_SW_Site);
		Mem_Free(grid_SW_VegProd);
	}
	
	if(UseSoils && UseSoilwat) {
		free_all_sxw_memory();
		Mem_Free(grid_SXW_ptrs);
		for( i=0; i < grid_Cells; i++)
			Mem_Free(grid_Soils[i].lyr);
		Mem_Free(grid_Soils);
	}
	if(UseDisturbances)
		Mem_Free(grid_Disturb);
	if(UseSeedDispersal) {
		ForEachSpecies(s) 
			if(Species[s]->use_me && Species[s]->use_dispersal) { 
				for(i = 0; i < grid_Cells; i++) {
					Mem_Free(grid_SD[s][i].cells);
					Mem_Free(grid_SD[s][i].prob);
					grid_SD[s][i].size = 0;
				}
				Mem_Free(grid_SD[s]);
			}

		for(i = 0; i < grid_Cells; i++)
			Mem_Free(grid_initSpecies[i].species_seed_avail);
		Mem_Free(grid_initSpecies);
	}

	stat_Free_Accumulators(); //free our memory we allocated for all the accumulators now that they're unnecessary to have
	
	for(i = 0; i < N_GRID_DIRECTORIES; i++) //frees the strings allocated in _init_grid_files()
    		Mem_Free(grid_directories[i]);
    	for(i = 0; i < N_GRID_FILES; i++)
    		Mem_Free(grid_files[i]);
    	
    	// freeing random memory that other parts of steppe/soilwat allocate... this isn't quite everything but it's a good start
    	parm_free_memory(); //frees memory allocated in ST_params.c
    
	ForEachSpecies(s) {
		if(!Species[s]->use_me) continue;
		_free_head(Species[s]->IndvHead);
		Mem_Free(Species[s]->kills);
		Mem_Free(Species[s]->seedprod);
	}
	    
	ForEachGroup(c)
		if(RGroup[c]->use_me) 
			Mem_Free(RGroup[c]->kills);
	    
	if(UseSoilwat) {
		Mem_Free(SXW.swc);
		Mem_Free(SXW.debugfile);
		for(i = 0; i < SW_Site.n_layers; i++)
			Mem_Free(SW_Site.lyr[i]);
		Mem_Free(SW_Site.lyr);
	}
    	
}

/***********************************************************/
static void _free_spinup_memory( void ) { 
	// frees spinup memory	
	
	GrpIndex c;
	SppIndex s;
	
	_free_spinup_globals();
	
	ForEachSpecies(s)
		if(Species[s]->use_me) Mem_Free(spinup_Species[s]);
	ForEachGroup(c)
		if(RGroup[c]->use_me) Mem_Free(spinup_RGroup[c]);
	
	Mem_Free(spinup_Succulent);
	Mem_Free(spinup_Env);
	Mem_Free(spinup_Plot);
	Mem_Free(spinup_Globals);
	if(UseSoilwat) {
		Mem_Free(spinup_SXW);
		Mem_Free(spinup_SW_Soilwat);
		Mem_Free(spinup_SW_Site);
		Mem_Free(spinup_SW_VegProd);
	}
	
	if(UseSoils && UseSoilwat) {
		Mem_Free(spinup_SXW_ptrs);
	}
    
	//if(sd_Option2a || sd_Option2b) {
		Mem_Free(soilTypes_Array);
		Mem_Free(grid_SoilTypes);
	//}
}

/***********************************************************/
static void _load_cell( int row, int col, int year, Bool useAccumulators ) {	
	// loads the specified cell into the global variables
	
	int cell = col + ( (row-1) * grid_Cols) - 1;  // converts the row/col into an array index
	int j;
	GrpIndex c;
	SppIndex s;
	//fprintf(stderr, " loading cell: %d; ", cell);
	if(useAccumulators)	
		stat_Load_Accumulators(cell, year);

	ForEachSpecies(s) {
		if(!Species[s]->use_me) continue;

		Mem_Free(Species[s]->kills);
		Mem_Free(Species[s]->seedprod);
		_free_head(Species[s]->IndvHead); //free_head() frees the memory allocated by the head and the memory allocated by each part of the linked list
		
		*Species[s] = grid_Species[s][cell];

		Species[s]->kills = Mem_Calloc(grid_Species[s][cell].max_age, sizeof(IntUS), "_load_cell(Species[s]->kills)");
		Species[s]->seedprod = Mem_Calloc(grid_Species[s][cell].viable_yrs, sizeof(RealF), "_load_cell(Species[s]->seedprod)");
		
		memcpy(Species[s]->kills, grid_Species[s][cell].kills, grid_Species[s][cell].max_age * sizeof(IntUS));
		memcpy(Species[s]->seedprod, grid_Species[s][cell].seedprod, grid_Species[s][cell].viable_yrs * sizeof(RealF));
		Species[s]->IndvHead = _copy_head(grid_Species[s][cell].IndvHead); //copy_head() deep copies the linked list structure (allocating memory when needed)... it will even allocate memory for the head of the list
	
	}
		
	ForEachGroup(c) {
		if(!RGroup[c]->use_me) continue;
		Mem_Free(RGroup[c]->kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)... we need to free it then reallocate it then memcpy it to get the deep copy we want
		
		*RGroup[c] = grid_RGroup[c][cell]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (ie copy the values in the pointers instead of the addresses) the pointers.  A shallow copy will copy over the values for every non-pointer (C itself does not inherently know how to deep copy, so we must code this behaviour).
		
		RGroup[c]->kills = Mem_Calloc(grid_RGroup[c][cell].max_age, sizeof(IntUS), "_load_cell(RGroup[c]->kills");
		memcpy(RGroup[c]->kills, grid_RGroup[c][cell].kills, grid_RGroup[c][cell].max_age * sizeof(IntUS));
	}
		
	Succulent = grid_Succulent[cell];
	Env = grid_Env[cell];
	Plot = grid_Plot[cell];
	Globals = grid_Globals[cell];
	
	if(UseSoilwat) {
		Mem_Free(SXW.transpTotal);
		if(SXW.swc != NULL) Mem_Free(SXW.swc);
		for(j = 0; j < SW_Site.n_layers; j++)
			Mem_Free(SW_Site.lyr[j]);
		Mem_Free(SW_Site.lyr);
			
		SXW = grid_SXW[cell];
		SW_Site = grid_SW_Site[cell];
		SW_Soilwat = grid_SW_Soilwat[cell];
		SW_VegProd = grid_SW_VegProd[cell];
		
		SXW.transpTotal = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealD), "_load_cell(SXW.transp)");
		SXW.swc = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealF), "_load_cell(SXW.swc)");
		
		memcpy(SXW.transpTotal, grid_SXW[cell].transpTotal, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
		memcpy(SXW.swc, grid_SXW[cell].swc, SXW.NPds * SXW.NSoLyrs * sizeof(RealF));
		
		SW_Site.lyr = Mem_Calloc(grid_SW_Site[cell].n_layers, sizeof(SW_LAYER_INFO *), "_load_cell(SW_Site.lyr)");
        for(j = 0; j < grid_SW_Site[cell].n_layers; j++) {
        	SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_load_cell(SW_Site.lyr[j])");
        	memcpy(SW_Site.lyr[j], grid_SW_Site[cell].lyr[j], sizeof(SW_LAYER_INFO));
        }
        
        if(UseSoils) load_sxw_memory(grid_SXW_ptrs[cell].roots_max, grid_SXW_ptrs[cell].rootsXphen, grid_SXW_ptrs[cell].roots_active, grid_SXW_ptrs[cell].roots_active_rel, grid_SXW_ptrs[cell].roots_active_sum, grid_SXW_ptrs[cell].phen, grid_SXW_ptrs[cell].prod_bmass, grid_SXW_ptrs[cell].prod_pctlive);
	}
}

/***********************************************************/
static void _load_spinup_cell( int cell ) {	
	// loads the specified cell into the global variables (from the spinup)
	
	int j;
	GrpIndex c;
	SppIndex s;

	ForEachSpecies(s) {
		if(!Species[s]->use_me) continue;

		Mem_Free(Species[s]->kills);
		Mem_Free(Species[s]->seedprod);
		_free_head(Species[s]->IndvHead); //free_head() frees the memory allocated by the head and the memory allocated by each part of the linked list
		
		*Species[s] = spinup_Species[s][cell];

		Species[s]->kills = Mem_Calloc(spinup_Species[s][cell].max_age, sizeof(IntUS), "_load_spinup_cell(Species[s]->kills)");
		Species[s]->seedprod = Mem_Calloc(spinup_Species[s][cell].viable_yrs, sizeof(RealF), "_load_spinup_cell(Species[s]->seedprod)");
		
		memcpy(Species[s]->kills, spinup_Species[s][cell].kills, spinup_Species[s][cell].max_age * sizeof(IntUS));
		memcpy(Species[s]->seedprod, spinup_Species[s][cell].seedprod, spinup_Species[s][cell].viable_yrs * sizeof(RealF));
		Species[s]->IndvHead = _copy_head(spinup_Species[s][cell].IndvHead); //copy_head() deep copies the linked list structure (allocating memory when needed)... it will even allocate memory for the head of the list
	
	}
		
	ForEachGroup(c) {
		if(!RGroup[c]->use_me) continue;
		Mem_Free(RGroup[c]->kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)... we need to free it then reallocate it then memcpy it to get the deep copy we want
		
		*RGroup[c] = spinup_RGroup[c][cell]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (ie copy the values in the pointers instead of the addresses) the pointers.  A shallow copy will copy over the values for every non-pointer (C itself does not inherently know how to deep copy, so we must code this behaviour).
		
		RGroup[c]->kills = Mem_Calloc(spinup_RGroup[c][cell].max_age, sizeof(IntUS), "_load_spinup_cell(RGroup[c]->kills");
		memcpy(RGroup[c]->kills, spinup_RGroup[c][cell].kills, spinup_RGroup[c][cell].max_age * sizeof(IntUS));
	}
		
	Succulent = spinup_Succulent[cell];
	Env = spinup_Env[cell];
	Plot = spinup_Plot[cell];
	Globals = spinup_Globals[cell];
	
	if(UseSoilwat) {
		Mem_Free(SXW.transpTotal);
		if(SXW.swc != NULL) Mem_Free(SXW.swc);
		for(j = 0; j < SW_Site.n_layers; j++)
			Mem_Free(SW_Site.lyr[j]);
		Mem_Free(SW_Site.lyr);
			
		SXW = spinup_SXW[cell];
		SW_Site = spinup_SW_Site[cell];
		SW_Soilwat = spinup_SW_Soilwat[cell];
		SW_VegProd = spinup_SW_VegProd[cell];
		
		SXW.transpTotal = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealD), "_load_spinup_cell(SXW.transp)");
		SXW.swc = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealF), "_load_spinup_cell(SXW.swc)");
		
		memcpy(SXW.transpTotal, spinup_SXW[cell].transpTotal, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
		memcpy(SXW.swc, spinup_SXW[cell].swc, SXW.NPds * SXW.NSoLyrs * sizeof(RealF));
		
		SW_Site.lyr = Mem_Calloc(spinup_SW_Site[cell].n_layers, sizeof(SW_LAYER_INFO *), "_load_spinup_cell(SW_Site.lyr)");
        for(j = 0; j < spinup_SW_Site[cell].n_layers; j++) {
        	SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_load_spinup_cell(SW_Site.lyr[j])");
        	memcpy(SW_Site.lyr[j], spinup_SW_Site[cell].lyr[j], sizeof(SW_LAYER_INFO));
        }
        
        if(UseSoils) load_sxw_memory(spinup_SXW_ptrs[cell].roots_max, spinup_SXW_ptrs[cell].rootsXphen, spinup_SXW_ptrs[cell].roots_active, spinup_SXW_ptrs[cell].roots_active_rel, spinup_SXW_ptrs[cell].roots_active_sum, spinup_SXW_ptrs[cell].phen, spinup_SXW_ptrs[cell].prod_bmass, spinup_SXW_ptrs[cell].prod_pctlive);
	}
}

/***********************************************************/
static void _save_cell( int row, int col, int year, Bool useAccumulators ) {	
	// saves the specified cell into the grid variables

	int cell = col + ( (row-1) * grid_Cols) - 1;  // converts the row/col into an array index
	int j;
	GrpIndex c;
	SppIndex s;
	//fprintf(stderr, "saving cell: %d\n", cell);

	if(useAccumulators)	
		stat_Save_Accumulators(cell, year);

	ForEachSpecies(s) {
		if(!Species[s]->use_me) continue;

		Mem_Free(grid_Species[s][cell].kills);
		Mem_Free(grid_Species[s][cell].seedprod);
		_free_head(grid_Species[s][cell].IndvHead);
		
		grid_Species[s][cell] = *Species[s];
		
		grid_Species[s][cell].kills = Mem_Calloc(Species[s]->max_age, sizeof(IntUS), "_save_cell(grid_Species[cell][s].kills)");
		grid_Species[s][cell].seedprod = Mem_Calloc(Species[s]->viable_yrs, sizeof(RealF), "_save_cell(grid_Species[cell][s].seedprod)");
		
		memcpy(grid_Species[s][cell].kills, Species[s]->kills, Species[s]->max_age * sizeof(IntUS));
		memcpy(grid_Species[s][cell].seedprod, Species[s]->seedprod, Species[s]->viable_yrs * sizeof(RealF));
		grid_Species[s][cell].IndvHead = _copy_head(Species[s]->IndvHead);
	}
		
	ForEachGroup(c) {
		if(!RGroup[c]->use_me) continue;
		Mem_Free(grid_RGroup[c][cell].kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)
		
		grid_RGroup[c][cell] = *RGroup[c]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (i.e. copy the values in the pointers instead of the addresses) the pointers
		
		grid_RGroup[c][cell].kills = Mem_Calloc(RGroup[c]->max_age, sizeof(IntUS), "_save_cell(grid_RGroup[cell][c].kills)");
		memcpy(grid_RGroup[c][cell].kills, RGroup[c]->kills, RGroup[c]->max_age * sizeof(IntUS));
	}
		
	grid_Succulent[cell] = Succulent;
	grid_Env[cell] = Env;
	grid_Plot[cell] = Plot;
	grid_Globals[cell] = Globals;
	
	if(UseSoilwat) {
		Mem_Free(grid_SXW[cell].transpTotal);
		Mem_Free(grid_SXW[cell].swc);
		for(j = 0; j < SW_Site.n_layers; j++)
			Mem_Free(grid_SW_Site[cell].lyr[j]);
		Mem_Free(grid_SW_Site[cell].lyr);	
	
		grid_SXW[cell] = SXW;
		grid_SW_Site[cell] = SW_Site;
		grid_SW_Soilwat[cell] = SW_Soilwat;
		grid_SW_VegProd[cell] = SW_VegProd;
		
		grid_SXW[cell].transpTotal = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealD), "_save_cell(grid_SXW[cell].transp)");
		grid_SXW[cell].swc = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealF), "_save_cell(grid_SXW[cell].swc)");
		
		memcpy(grid_SXW[cell].transpTotal, SXW.transpTotal, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
		memcpy(grid_SXW[cell].swc, SXW.swc, SXW.NPds * SXW.NSoLyrs * sizeof(RealF));
		
		grid_SW_Site[cell].lyr = Mem_Calloc(SW_Site.n_layers, sizeof(SW_LAYER_INFO *), "_save_cell(grid_SW_Site[cell].lyr[j])");
       		for(j = 0; j < SW_Site.n_layers; j++) {
        		grid_SW_Site[cell].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_save_cell(grid_SW_Site[cell].lyr[j])");
        		memcpy(grid_SW_Site[cell].lyr[j], SW_Site.lyr[j], sizeof(SW_LAYER_INFO));
        	}
        
        if(UseSoils) save_sxw_memory(grid_SXW_ptrs[cell].roots_max, grid_SXW_ptrs[cell].rootsXphen, grid_SXW_ptrs[cell].roots_active, grid_SXW_ptrs[cell].roots_active_rel, grid_SXW_ptrs[cell].roots_active_sum, grid_SXW_ptrs[cell].phen, grid_SXW_ptrs[cell].prod_bmass, grid_SXW_ptrs[cell].prod_pctlive);
	}
}

/***********************************************************/
static void _save_spinup_cell( int cell ) {	
	// saves the specified cell into the grid variables (from the spinup)

	int j;
	GrpIndex c;
	SppIndex s;
	
	ForEachSpecies(s) {
		if(!Species[s]->use_me) continue;

		Mem_Free(spinup_Species[s][cell].kills);
		Mem_Free(spinup_Species[s][cell].seedprod);
		_free_head(spinup_Species[s][cell].IndvHead);
		
		spinup_Species[s][cell] = *Species[s];
		
		spinup_Species[s][cell].kills = Mem_Calloc(Species[s]->max_age, sizeof(IntUS), "_save_cell(grid_Species[cell][s].kills)");
		spinup_Species[s][cell].seedprod = Mem_Calloc(Species[s]->viable_yrs, sizeof(RealF), "_save_cell(grid_Species[cell][s].seedprod)");
		
		memcpy(spinup_Species[s][cell].kills, Species[s]->kills, Species[s]->max_age * sizeof(IntUS));
		memcpy(spinup_Species[s][cell].seedprod, Species[s]->seedprod, Species[s]->viable_yrs * sizeof(RealF));
		spinup_Species[s][cell].IndvHead = _copy_head(Species[s]->IndvHead);
	}
		
	ForEachGroup(c) {
		if(!RGroup[c]->use_me) continue;
		Mem_Free(spinup_RGroup[c][cell].kills); //kills is the only pointer in the resourcegroup_st struct (which is what RGroup is defined as)
		
		spinup_RGroup[c][cell] = *RGroup[c]; //does a shallow copy, we have to do the freeing/malloc/memcpy to deep copy (i.e. copy the values in the pointers instead of the addresses) the pointers
		
		spinup_RGroup[c][cell].kills = Mem_Calloc(RGroup[c]->max_age, sizeof(IntUS), "_save_cell(grid_RGroup[cell][c].kills)");
		memcpy(spinup_RGroup[c][cell].kills, RGroup[c]->kills, RGroup[c]->max_age * sizeof(IntUS));
	}
		
	spinup_Succulent[cell] = Succulent;
	spinup_Env[cell] = Env;
	spinup_Plot[cell] = Plot;
	spinup_Globals[cell] = Globals;
	
	if(UseSoilwat) {
		Mem_Free(spinup_SXW[cell].transpTotal);
		Mem_Free(spinup_SXW[cell].swc);
		for(j = 0; j < SW_Site.n_layers; j++)
			Mem_Free(spinup_SW_Site[cell].lyr[j]);
		Mem_Free(spinup_SW_Site[cell].lyr);	
	
		spinup_SXW[cell] = SXW;
		spinup_SW_Site[cell] = SW_Site;
		spinup_SW_Soilwat[cell] = SW_Soilwat;
		spinup_SW_VegProd[cell] = SW_VegProd;
		
		spinup_SXW[cell].transpTotal = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealD), "_save_cell(grid_SXW[cell].transp)");
		spinup_SXW[cell].swc = Mem_Calloc(SXW.NPds * SXW.NSoLyrs, sizeof(RealF), "_save_cell(grid_SXW[cell].swc)");
		
		memcpy(spinup_SXW[cell].transpTotal, SXW.transpTotal, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
		memcpy(spinup_SXW[cell].swc, SXW.swc, SXW.NPds * SXW.NSoLyrs * sizeof(RealF));
		
		spinup_SW_Site[cell].lyr = Mem_Calloc(SW_Site.n_layers, sizeof(SW_LAYER_INFO *), "_save_cell(grid_SW_Site[cell].lyr[j])");
       		for(j = 0; j < SW_Site.n_layers; j++) {
        		spinup_SW_Site[cell].lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_save_cell(grid_SW_Site[cell].lyr[j])");
        		memcpy(spinup_SW_Site[cell].lyr[j], SW_Site.lyr[j], sizeof(SW_LAYER_INFO));
        	}
        
        if(UseSoils) save_sxw_memory(spinup_SXW_ptrs[cell].roots_max, spinup_SXW_ptrs[cell].rootsXphen, spinup_SXW_ptrs[cell].roots_active, spinup_SXW_ptrs[cell].roots_active_rel, spinup_SXW_ptrs[cell].roots_active_sum, spinup_SXW_ptrs[cell].phen, grid_SXW_ptrs[cell].prod_bmass, grid_SXW_ptrs[cell].prod_pctlive);
	}
}

/**************************************************************/
static Bool GetALine2( FILE *f, char buf[], int limit) {
  	//this is similar to the getaline function in filefuncs.c, except this one checks for carriage return characters and doesn't deal with whitespace/... (since excel writes them into .csv files for some aggravating reason)... this one is probably less efficient overall though.
  	//this is treating '\r', '\n', and '\r\n' all like they are valid line feed characters... in reality '\r' by itself should never be, but it's not like I can stop excel from outputting .csv files however the heck it feels like...
  	//only read limit characters
  	int i=0, aChar;
  	aChar = getc(f);
  	if(aChar == EOF) return FALSE;
  	while( (i < (limit-1)) && aChar != EOF && aChar != '\r' && aChar != '\n') {
  		buf[i++] = (char) aChar;
  		aChar = getc(f);
  	}
  	if(aChar == '\r') //this part handles the '\r\n' case
  		if(getc(f) != '\n')
  			fseek(f, -1, SEEK_CUR); //back up one character in the file because we didn't find a new-line character
  
  	buf[i] = '\0';
  	return TRUE;
}

/***********************************************************/
static void _read_disturbances_in( void ) {
	// reads the grid disturbances input file
	// the file should be something like: "cell,use_fecal_pats,use_ant_mounds,use_animal_burrows,kill_yr"
	// there should be no spaces in between, just commas separating the values
	// kill_yr will overwrite the kill year for each RGroup in the cell (0 means don't use, a # > 0 means kill everything at this year)

	FILE *f;
	char buf[1024];
	int i, cell, num;
	
    	f = OpenFile(grid_files[2], "r");
    
    	GetALine2(f, buf, 1024); // gets rid of the first line (since it just defines the columns)
    	for(i = 0; i < grid_Cells; i++) {
    		if(!GetALine2(f, buf, 1024)) break;

    		num = sscanf( buf, "%d,%d,%d,%d,%d,%d,%d", &cell, &grid_Disturb[i].choices[0], &grid_Disturb[i].choices[1], &grid_Disturb[i].choices[2], &grid_Disturb[i].kill_yr, &grid_Disturb[i].killfrq, &grid_Disturb[i].extirp);
		if(num != 7)
			LogError(logfp, LOGFATAL, "Invalid %s file line %d wrong", grid_files[2], i+2);
	}
    	if(i != grid_Cells)
    		LogError(logfp, LOGFATAL, "Invalid %s file wrong number of cells", grid_files[2]);	
    
    	CloseFile(&f);
}

/***********************************************************/
static int _get_value_index( char* s, char seperator, int nSeperators ) {
	//pretty much this function goes through s until it find nSeperators worth of seperators and then returns the index of the next character
	//this is used to do most of the parsing in the _read_soils_in() function
	int i = 0, sep = 0;
	while(*s) {
		i++;
		if(*s++ == seperator) //checks if the current char equals the seperator... and then increments the char pointer
			if(++sep == nSeperators) //needs ++sep, not sep++
				break;
	}
	return i;
}

/***********************************************************/
static void _read_soils_in( void ) {
	// reads the grid soils input
	// the file should be something like: "cell,copy_cell,copy_which,num_layers,..."
	// there should be no spaces in between, just commas separating the values
	// this function reads in pretty much a .csv file, but it will not account for all of the possibilities that a .csv file could be written as (as accounting for all these possibilities would take a while to code and be unproductive) so keep that in mind
		
	FILE *f;
	char buf[4096];
	int i, j, k, cell, num, do_copy, copy_cell, num_layers, depth, depthMin;
	float d[11];

	if(sd_Option2a || sd_Option2b) 
		nSoilTypes = 0; //initialize our soil types counter

	f = OpenFile(grid_files[3], "r");
	
	GetALine2(f, buf, 4096); // gets rid of the first line (since it just defines the columns)... it's only there for user readability
	for(i = 0; i < grid_Cells; i++) {
		if(!GetALine2(f, buf, 4096)) break;
		
		num = sscanf( buf, "%d,%d,%d,%d", &cell, &do_copy, &copy_cell, &num_layers );
		if(num != 4)
			LogError(logfp, LOGFATAL, "Invalid %s file", grid_files[3]);
		int stringIndex = _get_value_index(buf, ',', 4); //gets us the index of the string that is right after what we just parsed in
		
		if(num_layers > MAX_LAYERS)
			LogError(logfp, LOGFATAL, "Invalid %s file line %d num_layers (%d) exceeds MAX_LAYERS (%d)", grid_files[3], i+2, num_layers, MAX_LAYERS); 

		if(do_copy == 1 && copy_cell > -1 && copy_cell < grid_Cells && cell != 0 && copy_cell < cell) { //copy this cells values from a previous cell's
			grid_Soils[i].lyr = Mem_Calloc(grid_Soils[copy_cell].num_layers, sizeof(Grid_Soil_Lyr), "_read_soils_in()");
			for(j = 0; j < grid_Soils[copy_cell].num_layers; j++)
				grid_Soils[i].lyr[j] = grid_Soils[copy_cell].lyr[j];
			grid_Soils[i].num_layers = grid_Soils[copy_cell].num_layers;

			if(sd_Option2a || sd_Option2b)
				grid_SoilTypes[cell] = grid_SoilTypes[copy_cell];

			continue;
		} else if(do_copy == 1)
			LogError(logfp, LOGFATAL, "Invalid %s file line %d invalid copy_cell attempt", grid_files[3], i+2);

		if(sd_Option2a || sd_Option2b) {
			grid_SoilTypes[cell] = nSoilTypes;
			soilTypes_Array[nSoilTypes] = cell;
			nSoilTypes++;
		}

		depthMin = 0;
		grid_Soils[i].num_layers = num_layers;
		grid_Soils[i].lyr = Mem_Calloc(num_layers, sizeof(Grid_Soil_Lyr), "_read_soils_in()");
		for(j = 0; j < num_layers; j++) {
			 //the idea behind using &buf[stringIndex] is that we start scanning at the point in the string that is right after what we just parsed... the & is there because we have to send sscanf the pointer that points to that location
			num = sscanf( &buf[stringIndex], "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &depth, &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10] );
			if(num != 12)
				LogError(logfp, LOGFATAL, "Invalid %s file line %d invalid soil layer input", grid_files[3], i+2);
			
			k = stringIndex;
			stringIndex += _get_value_index(&buf[stringIndex], ',', 12); //updates the index of the string that we are at
			if(k == stringIndex)
				LogError(logfp, LOGFATAL, "Invalid %s file line %d not enough soil layers", grid_files[3], i+2);
			
			for(k = 0; k < 11; k++) 
				grid_Soils[i].lyr[j].data[k] = d[k];
			grid_Soils[i].lyr[j].width = depth-depthMin;
			depthMin = depth;
		}
	}
	
	if(i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file, not enough cells", grid_files[3]);
		
	/*for(i = 0; i < grid_Cells; i++) {
		printf("cell %d:\n", i);
		for(j = 0; j < grid_Soils[i].num_layers; j++) {
			printf("layer %d : %d", j, grid_Soils[i].lyr[j].width);
			for(k = 0; k < 11; k++) printf(" %f", grid_Soils[i].lyr[j].data[k]);
			printf("\n");
		}
	}*/
	
	CloseFile(&f);
}

/***********************************************************/
static void _init_soil_layers(int cell, int isSpinup) {
	// initializes the soilwat soil layers for the cell correctly based upon the input gathered from our grid_soils input file
	// pretty much takes the data from grid_Soils (read in in _read_soils_in()) and converts it to what SW_Site needs...
	// this function does generally the same things that the _read_layers() function in SW_Site.c does, except that it does it in a way that lets us use it in the grid...
	int i, j;
	i = cell;
	
	Bool evap_ok = TRUE, transp_ok_forb = TRUE, transp_ok_tree = TRUE, transp_ok_shrub = TRUE, transp_ok_grass = TRUE; /* mitigate gaps in layers */
	
	for(j = 0; j < SW_Site.n_layers+SW_Site.deepdrain; j++)
		Mem_Free(SW_Site.lyr[j]);
	Mem_Free(SW_Site.lyr);
			
	SW_Site.n_layers = grid_Soils[i].num_layers;
	SW_Site.n_evap_lyrs = SW_Site.n_transp_lyrs_forb = SW_Site.n_transp_lyrs_tree = SW_Site.n_transp_lyrs_shrub = SW_Site.n_transp_lyrs_grass = 0;
			
	SW_Site.lyr = Mem_Calloc(SW_Site.n_layers+SW_Site.deepdrain, sizeof(SW_LAYER_INFO *), "_init_grid_globals()");
    	for(j = 0; j < SW_Site.n_layers; j++) {
        	SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_init_grid_globals()");
        		
        	//indexes (for grid_Soils[i].lyr[j].data):
        	//0		   1				2		3	  		4			5			6			7	   8		9		10
        	//matricd	gravel_content  evco  	trco_grass  trco_shrub  trco_tree  	trco_forb	%sand  %clay imperm soiltemp
        	SW_Site.lyr[j]->width = grid_Soils[i].lyr[j].width;
        	SW_Site.lyr[j]->soilBulk_density = grid_Soils[i].lyr[j].data[0];
        	SW_Site.lyr[j]->fractionVolBulk_gravel = grid_Soils[i].lyr[j].data[1];
        	SW_Site.lyr[j]->evap_coeff = grid_Soils[i].lyr[j].data[2];
        	SW_Site.lyr[j]->transp_coeff_grass = grid_Soils[i].lyr[j].data[3];
        	SW_Site.lyr[j]->transp_coeff_shrub = grid_Soils[i].lyr[j].data[4];
        	SW_Site.lyr[j]->transp_coeff_tree = grid_Soils[i].lyr[j].data[5];
        	SW_Site.lyr[j]->transp_coeff_forb = grid_Soils[i].lyr[j].data[6];
        	SW_Site.lyr[j]->fractionWeightMatric_sand = grid_Soils[i].lyr[j].data[7];
        	SW_Site.lyr[j]->fractionWeightMatric_clay = grid_Soils[i].lyr[j].data[8];
        	SW_Site.lyr[j]->impermeability = grid_Soils[i].lyr[j].data[9];
        	SW_Site.lyr[j]->my_transp_rgn_tree = 0;
        	SW_Site.lyr[j]->my_transp_rgn_shrub = 0;
        	SW_Site.lyr[j]->my_transp_rgn_grass = 0;
        	SW_Site.lyr[j]->sTemp = grid_Soils[i].lyr[j].data[10];
        
		if ( evap_ok ) {
			if ( GT(SW_Site.lyr[j]->evap_coeff, 0.0) )
				SW_Site.n_evap_lyrs++;
			else
				evap_ok = FALSE;
		}
		if ( transp_ok_tree ) {
			if ( GT(SW_Site.lyr[j]->transp_coeff_tree, 0.0) )
				SW_Site.n_transp_lyrs_tree++;
			else
				transp_ok_tree = FALSE;
		}
		if ( transp_ok_shrub ) {
			if ( GT(SW_Site.lyr[j]->transp_coeff_shrub, 0.0) )
				SW_Site.n_transp_lyrs_shrub++;
			else
				transp_ok_shrub = FALSE;
		}
		if ( transp_ok_grass ) {
			if ( GT(SW_Site.lyr[j]->transp_coeff_grass, 0.0) )
				SW_Site.n_transp_lyrs_grass++;
			else
				transp_ok_grass = FALSE;
		}
		if (transp_ok_forb) {
			if (GT(SW_Site.lyr[j]->transp_coeff_forb, 0.0))
				SW_Site.n_transp_lyrs_forb++;
			else
				transp_ok_forb = FALSE;
		}
		water_eqn(SW_Site.lyr[j]->fractionVolBulk_gravel, SW_Site.lyr[j]->fractionWeightMatric_sand, SW_Site.lyr[j]->fractionWeightMatric_clay, j);//in SW_Site.c, called to initialize some layer data...
		SW_Site.lyr[j]->swcBulk_fieldcap = SW_SWPmatric2VWCBulk(SW_Site.lyr[j]->fractionVolBulk_gravel, 0.333, j) * SW_Site.lyr[j]->width;
		SW_Site.lyr[j]->swcBulk_wiltpt = SW_SWPmatric2VWCBulk(SW_Site.lyr[j]->fractionVolBulk_gravel, 15, j) * SW_Site.lyr[j]->width;
		//From calculate_soilBulkDensity in SW_Site.c
		SW_Site.lyr[j]->soilBulk_density = SW_Site.lyr[j]->soilMatric_density * (1 - SW_Site.lyr[j]->fractionVolBulk_gravel) + (SW_Site.lyr[j]->fractionVolBulk_gravel * 2.65);
		//already checked for max_layers condition
	}
    if(SW_Site.deepdrain) {
    	SW_Site.n_layers++;
    	SW_Site.lyr[j] = Mem_Calloc(1, sizeof(SW_LAYER_INFO), "_init_grid_globals()");
    	SW_Site.lyr[j]->width = 1.0;
    }

	init_site_info(); //in SW_Site.c, called to initialize layer data...
	
	free_all_sxw_memory();
	_init_SXW_inputs(FALSE); //we call this so that SXW can set the correct sizes/values up for the memory dynamically allocated in sxw.c

	if(!isSpinup) {

		grid_SXW_ptrs[i].roots_max = Mem_Calloc(SXW.NGrps * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].rootsXphen = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].roots_active = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].roots_active_rel = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].roots_active_sum = Mem_Calloc(SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].phen = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].prod_bmass = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].prod_pctlive = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "_init_soil_layers()");

		save_sxw_memory(grid_SXW_ptrs[i].roots_max, grid_SXW_ptrs[i].rootsXphen, grid_SXW_ptrs[i].roots_active, grid_SXW_ptrs[i].roots_active_rel, grid_SXW_ptrs[i].roots_active_sum, grid_SXW_ptrs[i].phen, grid_SXW_ptrs[i].prod_bmass, grid_SXW_ptrs[i].prod_pctlive);
	} else {

		spinup_SXW_ptrs[i].roots_max = Mem_Calloc(SXW.NGrps * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].rootsXphen = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].roots_active = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].roots_active_rel = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].roots_active_sum = Mem_Calloc(SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "_init_soil_layers()");
		spinup_SXW_ptrs[i].phen = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].prod_bmass = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "_init_soil_layers()");
		grid_SXW_ptrs[i].prod_pctlive = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "_init_soil_layers()");

		save_sxw_memory(spinup_SXW_ptrs[i].roots_max, spinup_SXW_ptrs[i].rootsXphen, spinup_SXW_ptrs[i].roots_active, spinup_SXW_ptrs[i].roots_active_rel, spinup_SXW_ptrs[i].roots_active_sum, spinup_SXW_ptrs[i].phen, grid_SXW_ptrs[i].prod_bmass, grid_SXW_ptrs[i].prod_pctlive);
	}
}

/***********************************************************/
static float _read_a_float(FILE *f, char *buf, const char *filename, const char *descriptor) {
	//small function to reduce code duplication in the _read_seed_dispersal_in() function...
	//f should be already open, and all of the character arrays should be pre-allocated before calling the function...
	float result;

	if(!GetALine(f, buf))
		LogError(logfp, LOGFATAL, "Invalid %s file: %s", filename, descriptor);
	if(sscanf(buf, "%f", &result) != 1) 
		LogError(logfp, LOGFATAL, "Invalid %s file: %s", filename, descriptor); 
		
	return result;
} 

/***********************************************************/
static float _cell_dist(int row1, int row2, int col1, int col2, float cellLen) {
	//returns the distance between the two grid cells
	if(row1 == row2) {
		return (abs(col1-col2) * cellLen);
	} else if(col1 == col2) {
		return (abs(row1-row2) * cellLen);
	} else { // row1 != row2 && col1 != col2
		//the problem can be thought of in terms of a right triangle...
		//using the pythagorean theorem: c = sqrt(a^2 + b^2)... c (the hypotenuse) represents the distance that we need.  a is the distance between columns and b is the distance between rows.
		return sqrt(pow(abs(col1-col2)*cellLen, 2.0) + pow(abs(row1-row2)*cellLen, 2.0));
	}
}

/***********************************************************/
static void _read_seed_dispersal_in( void ) {
	// reads the grid seed dispersal input file and sets up grid_SD with the correct values and probabilities

	FILE *f;
	char buf[1024];
	float sd_Rate, H, VW, VT, MAXD, plotLength, d, pd;
	int maxCells, i, j, k, MAXDP, row, col, cell;
	SppIndex s;
	
    	// read in the seed dispersal input file to get the constants that we need
    	f = OpenFile(grid_files[4], "r");
    
    	VW = _read_a_float(f, buf, grid_files[4], "VW line");
	
	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_DoOutput) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: seed dispersal output line\n", grid_files[4]);
	   
	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_MakeHeader) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: seed dispersal make header line\n", grid_files[4]);

	GetALine(f, buf);
	if(sscanf(buf, "%c", &sd_Sep) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: seed dispersal seperator line\n", grid_files[4]);

	if(sd_Sep == 't') //dealing with tab and space special cases...
		sd_Sep = '\t';
	else if(sd_Sep == 's')
		sd_Sep = ' ';

	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_NYearsSeedsAvailable) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 1 line\n", grid_files[4]);

	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_Option1a) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 1a line\n", grid_files[4]);

	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_Option1b) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 1b line\n", grid_files[4]);

	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_Option2a) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 2a line\n", grid_files[4]);

	GetALine(f, buf);
	if(sscanf(buf, "%d", &sd_Option2b) != 1)
		LogError(logfp, LOGFATAL, "Invalid %s file: option 2b line\n", grid_files[4]);

	if( (sd_Option1a && (sd_Option1b || sd_Option2a || sd_Option2b)) || (sd_Option2a && (sd_Option1a || sd_Option1b || sd_Option2b)) )
		LogError(logfp, LOGFATAL, "Invalid %s file: conflicting options chosen\n", grid_files[4]);
		

    	CloseFile(&f);
	
	ForEachSpecies(s) {
		// set up grid_SD with the seed dispersal probabilities needed later on...
		H = Species[s]->sd_H;
		VT = Species[s]->sd_VT;
		MAXD = ((H * VW) / VT) / 100.0; // divide by 100.0 because we want it in meters, not centimeters
		sd_Rate = -(log(0.005)/MAXD); //sd_Rate is the seed dispersal rate... 0.005 = exp(-RATE*MAXD) => RATE = -(ln(0.005)/MAXD)
		
		plotLength = sqrt(Globals.plotsize);
		MAXDP = (int) ceil(MAXD / plotLength); //MAXD in terms of plots... rounds up to the nearest integer
		maxCells = (int) pow((MAXDP*2) + 1.0, 2.0); //gets the maximum number of cells that a grid cell can possibly disperse seeds to... it ends up being more then the maximum actually...
		if(grid_Cells < maxCells)
			maxCells = grid_Cells;
		if(! (Species[s]->use_me && Species[s]->use_dispersal) ) continue;
		
		for(i = 0; i < grid_Cells; i++) {
			grid_SD[s][i].cells = Mem_Calloc(maxCells, sizeof(int), "_read_seed_dispersal_in()"); //the cell number
			grid_SD[s][i].prob = Mem_Calloc(maxCells, sizeof(float), "_read_seed_dispersal_in()"); //the probability that the cell will disperse seeds to this distance
			grid_SD[s][i].size = 0; //refers to the number of cells reachable...
		}

		for(row = 1; row <= grid_Rows; row++)
			for(col = 1; col <= grid_Cols; col++) {

				cell = col + ( (row-1) * grid_Cols) - 1;
				k = 0;

				for(i = 1; i <= grid_Rows; i++)
					for(j = 1; j <= grid_Cols; j++) {
						if(i == row && j == col) continue;

						d = _cell_dist(i, row, j, col, plotLength); //distance
						pd = (d > MAXD) ? (0.0) : (exp(-sd_Rate*d)); //dispersal probability

						if(!ZRO(pd)) {
							grid_SD[s][cell].cells[k] = i + ( (j-1) * grid_Cols) - 1;
							grid_SD[s][cell].prob[k] = pd;
							grid_SD[s][cell].size++;
							k++;
							//fprintf(stderr, "cell: %d; i: %d; j: %d; dist: %f; pd: %f %d %d\n", i + ( (j-1) * grid_Cols) - 1, i, j, d, pd, row, col); 
						}
					}
				//fprintf(stderr, "size %d index %d maxsize %d\n", grid_SD[cell].size, cell, maxCells);
			}

		for(i = 0; i < grid_Cells; i++)
			if(grid_SD[s][i].size > 0) {
				grid_SD[s][i].cells = Mem_ReAlloc(grid_SD[s][i].cells, grid_SD[s][i].size * sizeof(int));
				grid_SD[s][i].prob = Mem_ReAlloc(grid_SD[s][i].prob, grid_SD[s][i].size * sizeof(float));
			}
	}
}

/***********************************************************/
static void _do_seed_dispersal(void ) {
	float biomass, randomN, LYPPT, presentProb, receivedProb;
	int i, j, germ, sgerm, year;
	SppIndex s;

	if(Globals.currYear == 1 && !sd_Option1a && !sd_Option1b) { //since we have no previous data to go off of, use the current years...
		for(i = 0; i < grid_Cells; i++) 
			ForEachSpecies(s) {
				if( ! (Species[s]->use_me && Species[s]->use_dispersal) ) continue;	
				grid_Species[s][i].allow_growth = grid_Species[s][i].sd_sgerm = 1;	// since it's the first year, we have to allow growth...
				if(UseDisturbances)
					if(1 == grid_Disturb[i].kill_yr)
						grid_Species[s][i].allow_growth = 0;
				grid_SD[s][i].lyppt = grid_Env[i].ppt;
			}
	} else {
		// figure out whether or not to allow growth for the current year... based upon whether the species already has plants or germination allowed this year and seeds received last year...
		ForEachSpecies(s) {
			
			if(! (Species[s]->use_me && Species[s]->use_dispersal) ) continue;
			
			// germination probability
			randomN = RandUni(); 
			germ = LE(randomN, Species[s]->seedling_estab_prob);

			year = Globals.currYear - 1;
			
			for(i = 0; i < grid_Cells; i++) {

				if(sd_Option1a && Globals.currYear <= sd_NYearsSeedsAvailable) {
					grid_SD[s][i].seeds_present = 1;	
				} else if(sd_Option1b && Globals.currYear <= sd_NYearsSeedsAvailable && grid_initSpecies[i].species_seed_avail[s]) {
					grid_SD[s][i].seeds_present = 1;
				}

				sgerm = (grid_SD[s][i].seeds_present || grid_SD[s][i].seeds_received) && germ; //refers to whether the species has seeds available from the previous year and conditions are correct for germination this year
				grid_Species[s][i].allow_growth = FALSE;
				biomass = grid_Species[s][i].relsize * grid_Species[s][i].mature_biomass;

				if(UseDisturbances) {
					if((sgerm ||  year < grid_Disturb[i].kill_yr || grid_Disturb[i].kill_yr <= 0 ||  GT(biomass, 0.0)) && (year != grid_Disturb[i].kill_yr))
						grid_Species[s][i].allow_growth = TRUE;
				} else if(sgerm || GT(biomass, 0.0))
					grid_Species[s][i].allow_growth = TRUE;
				grid_Species[s][i].sd_sgerm = sgerm; //based upon whether we have received/produced seeds that germinated
				//if(grid_Species[s][i].allow_growth == TRUE &&  i == 52 && s == 0 && Globals.currIter == 1)
				//	printf("%s allow_growth:%d year:%d sgerm:%d iter:%d\n", grid_Species[s][i].name, grid_Species[s][i].allow_growth, year, sgerm, Globals.currIter);
			}
		}

	}

	// calculate whether or not seeds were received/produced this year, this data is used the next time the function is called
	ForEachSpecies(s) {
		if(! (Species[s]->use_me && Species[s]->use_dispersal) ) continue;

		IndivType* indiv;

		// figure out which species in each cell produced seeds...
		for(i = 0; i < grid_Cells; i++) {
			grid_SD[s][i].seeds_present = grid_SD[s][i].seeds_received = grid_Species[s][i].received_prob = 0;
			
			biomass = 0;	//getting the biggest individual in the species...
			ForEachIndiv(indiv, &grid_Species[s][i])
				if(indiv->relsize * grid_Species[s][i].mature_biomass > biomass)
					biomass = indiv->relsize * grid_Species[s][i].mature_biomass;  

			if(GE(biomass, grid_Species[s][i].mature_biomass * grid_Species[s][i].sd_Param1)) {
				randomN = RandUni();

				LYPPT = grid_SD[s][i].lyppt;
				float PPTdry = grid_Species[s][i].sd_PPTdry, PPTwet = grid_Species[s][i].sd_PPTwet;
				float Pmin = grid_Species[s][i].sd_Pmin, Pmax = grid_Species[s][i].sd_Pmax;

				//p3 = Pmin, if LYPPT < PPTdry 
				//p3 = 1 - (1-Pmin) * exp(-d * (LYPPT - PPTdry)) with d = - ln((1 - Pmax)/(1 - Pmin)) / (PPTwet - PPTdry), if PPTdry <= LYPPT <= PPTwet
				//p3 = Pmax, if LYPPT > PPTwet
			
				presentProb = 0.0;
				if(PPTdry <= LYPPT && LYPPT <= PPTwet) {
					float d = -log(((1 - Pmax)/(1 - Pmin))) / (PPTwet - PPTdry); //log is the natural log in STD c's math.h
					presentProb = 1 - (1-Pmin) * exp((-d * (LYPPT - PPTdry)));
				} else if(LYPPT < PPTdry)
					presentProb = Pmin;
				else if(LYPPT > PPTwet)
					presentProb = Pmax;

				if(LE(randomN, presentProb))
					grid_SD[s][i].seeds_present = 1;
			}
			//if(i == 0) printf("cell: %d lyppt: %f\n", i, grid_SD[i].lyppt);
		}

		// figure out which species in each cell received seeds...
		for(i = 0; i < grid_Cells; i++) {
			if(grid_SD[s][i].seeds_present) continue;
			receivedProb = 0;
			
			for(j = 0; j < grid_SD[s][i].size; j++)
				if(grid_SD[s][grid_SD[s][i].cells[j]].seeds_present)
					receivedProb += grid_SD[s][i].prob[j];

			randomN = RandUni();
			if(LE(randomN, receivedProb) && !ZRO(receivedProb)) 
				grid_SD[s][i].seeds_received = 1;
			else
				grid_SD[s][i].seeds_received = 0;

			grid_Species[s][i].received_prob = receivedProb;
		}
	}
} 

/***********************************************************/
static void _set_sd_lyppt(int row, int col) {
	int cell = col + ( (row-1) * grid_Cols) - 1;
	SppIndex s;
	
	ForEachSpecies(s)
		if(Species[s]->use_me && Species[s]->use_dispersal) 
			grid_SD[s][cell].lyppt = grid_Env[cell].ppt;
}

/***********************************************************/
static void _kill_groups_and_species( void ) {
	// kills all groups and species, leaving the ability to regenerate at a later time...
	// I first assumed that killing all of the RGroups would also kill all of the Species, but for some reason it wasn't getting everything in very rare cases (which is weird, b/c I thought that a Species should always belong to an RGroup)... so instead we are now just killing all of the species individually
	GrpIndex c;
	SppIndex sp;
	
	ForEachGroup(c)
		RGroup[c]->pr = 0.0; // reset the pr, so our output doesn't look weird
									
	ForEachSpecies(sp) {
		if(!Species[sp]->use_me) continue;
		Species_Kill(sp);
		Species[sp]->relsize = 0.0;
	}
}

/***********************************************************/
static int _do_grid_disturbances( int row, int col ) {
	// return 1 if a disturbance occurs, else return 0
	if(UseDisturbances) {
		int cell = col + ( (row-1) * grid_Cols) - 1;
		if(Globals.currYear == grid_Disturb[cell].kill_yr) {
			//basically if a disturbance occurs, we kill everything and then don't allow any species to grow for the year
			_kill_groups_and_species();
			SppIndex s;
			ForEachSpecies(s)
				if(Species[s]->use_me && Species[s]->use_dispersal)
					Species[s]->allow_growth = FALSE;
			return 1;
		}
	}
	return 0;
}

/***********************************************************/
static void _read_init_species( void ) {
	// reads the grid init species input
	// the file should be something like: "cell,copy_cell,copy_which,use_SpinUp,(all the species names seperated by a comma)"
	// there should be no spaces in between, just commas separating the values (ie it should be a .csv file, but does not account for all of the possibilities that a .csv file could be)

	FILE *f;
	int i, j, num, cell, do_copy, copy_cell, use_SpinUp, seeds_Avail;
	char buf[4096];
	
	//open the file/do the reading	
	f = OpenFile(grid_files[5], "r"); //grid_files[5] is the grid_initSpecies.csv file

	GetALine2(f, buf, 4096); // gets rid of the first line (since it just defines the columns)... it's only there for user readability
	for(i = 0; i < grid_Cells; i++) {
		if(!GetALine2(f, buf, 4096)) break;

		num = sscanf( buf, "%d,%d,%d,%d", &cell, &do_copy, &copy_cell, &use_SpinUp );
		if(num != 4)
			LogError(logfp, LOGFATAL, "Invalid %s file", grid_files[5]);

		grid_initSpecies[i].use_SpinUp = use_SpinUp;
		
		int stringIndex = _get_value_index(buf, ',', 4); //gets us the index of the string that is right after what we just parsed in
		
		if(do_copy == 1 && copy_cell > -1 && copy_cell < grid_Cells && cell != 0 && copy_cell < cell) { //copy this cells values from a previous cell's
			for(j = 0; j < Globals.sppCount; j++)
				grid_initSpecies[i].species_seed_avail[j] = grid_initSpecies[copy_cell].species_seed_avail[j];
			grid_initSpecies[i].use_SpinUp = grid_initSpecies[copy_cell].use_SpinUp;
			continue;
		} else if(do_copy == 1)
			LogError(logfp, LOGFATAL, "Invalid %s file line %d invalid copy_cell attempt", grid_files[5], i+2);

		//going through each species
		SppIndex s;
		ForEachSpecies(s) {
			num = sscanf( &buf[stringIndex], "%d,", &seeds_Avail);
			if(num != 1)
				LogError(logfp, LOGFATAL, "Invalid %s file line %d invalid species input", grid_files[5], i+2);

			grid_initSpecies[i].species_seed_avail[s] = seeds_Avail;
			stringIndex += _get_value_index(&buf[stringIndex], ',', 1);
		}
	}

	if(i != grid_Cells)
		LogError(logfp, LOGFATAL, "Invalid %s file, not enough cells", grid_files[5]);


	CloseFile(&f);	
}
