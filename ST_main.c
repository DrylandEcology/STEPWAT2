/**
 * \file ST_main.c
 * \brief Main program loop and argument processing
 * 
 *  History:
 *     (6/15/2000) -- INITIAL CODING - cwb
 *     15-Apr-02 (cwb) -- added code to interface with SOILWAT
 *	   5-24-2013 (DLM) -- added gridded option to program... see ST_grid.c
 *               source file for the rest of the gridded code 
 * 
 * \author CWB (initial programming)
 * \author DLM (added gridded mode option)
 * \author Kyle Palmquist
 * \author Chandler Haukap
 * \date 15 April 2000 (initial programming)
 * \ingroup STEPPE
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "ST_steppe.h"
#include "sw_src/include/generic.h"
#include "sw_src/include/filefuncs.h"
#include "sw_src/include/myMemory.h"
#include "sw_src/include/SW_Main_lib.h"
#include "sw_src/include/SW_VegProd.h"
#include "sw_src/include/SW_Control.h"
#include "sw_src/include/SW_Defines.h"
#include "sw_src/include/SW_Domain.h"

#include "sxw_funcs.h"
#include "sxw.h"
#include "sw_src/include/SW_Output.h"
#include "sw_src/include/SW_Output_outtext.h"
#include "sw_src/include/SW_Output_outarray.h"
#include "sw_src/include/rands.h"
#include "ST_functions.h" // externs `environs_rng`, `resgroups_rng`, `species_rng`
#include "ST_spinup.h"
#include "ST_stats.h"
#include "ST_progressBar.h"
#include "ST_seedDispersal.h" // externs `UseSeedDispersal`, `dispersal_rng`
#include "ST_mortality.h" // externs `mortality_rng`, `*_SomeKillage`
#include "ST_grid.h" // externs `grid_rng`




/************* External Function Declarations **************/
/***********************************************************/

  void rgroup_Grow( void);
  void rgroup_Establish( void) ;
  void rgroup_IncrAges( void);
  void rgroup_PartResources( void);

  void parm_Initialize(void);
  void parm_SetFirstName( char *s);

  void output_Bmass_Yearly( Int year );
  void output_Mort_Yearly( void );

  void save_annual_species_relsize(void);

  void files_init(void);
  void maxrgroupspecies_init(void);

/*************** Local Function Declarations ***************/
/***********************************************************/
void Plot_Initialize( void);

#ifndef STDEBUG
/** \brief Prints a description of flags then exits the program.
 *  
 * Meant to be used when the user inputs an undefined flag.
 * 
 * \sa init_args(int argc, char **argv)
 */
static void usage(void) {
  char *s ="STEPPE plant community dynamics (SGS-LTER Jan-04).\n"
           "   Usage : steppe [-d startdir] [-f files.in] [-q] [-e] [-o] [-g]\n"
           "      -d : supply working directory (default=.)\n"
           "      -f : supply list of input files (default=files.in)\n"
           "      -q : quiet mode, don't print message to check logfile.\n"
           "      -p : prints progress bar\n"
           "      -e : echo initialization results to logfile\n"
           "      -o : write SOILWAT output to output files. Contains average over all iterations and standard deviation.\n"
           "      -g : use gridded mode\n"
           "      -i : write SOILWAT output to output files for each iteration\n" // dont need to set -o flag to use this flag
		   "-STdebug : generate sqlite database with STEPWAT information\n";
  fprintf(stderr,"%s", s);
  exit(0);
}

static void init_args(int argc, char **argv);
static void check_log(void);

void allocate_Globals(void);
void deallocate_Globals(Bool isGriddedMode);
#endif

/* a couple of debugging routines */
void check_sizes(const char *);



/************ External Variable Definitions  ***************/
/*              see ST_globals.h                       */
/***********************************************************/


/** \brief Global struct holding species-specific variables. */
SpeciesType  **Species;
/** \brief Global struct holding rgroup-specific variables. */
GroupType    **RGroup;
/** \brief Global struct holding succulent-specific constants. */
SucculentType  *Succulent;
/** \brief Global struct holding environment-specific variables. */
EnvType        *Env;
/** \brief Global struct holding plot-specific variables. */
PlotType       *Plot;
/** \brief Global struct holding global variables. */
ModelType      *Globals;
/** \brief Global struct holding biomass output flags. */
GlobalType     SuperGlobals;
/** \brief Global struct holding variables describing the domain */
SW_DOMAIN      SoilWatDomain;
/** \brief Global struct holding SOILWAT2 variables*/
SW_RUN         SoilWatRun;
/** \brief Global struct holding log information (used by STEPWAT2) */
LOG_INFO       LogInfo;
/** \brief Global struct holding log information (used by SOILWAT2) */
LOG_INFO       LogInfoSW;
/** \brief Local booleans to echo inputs/any output (used by SOILWAT2) */
Bool           EchoInits;

BmassFlagsType BmassFlags;
/** \brief Global struct holding biomass quantile mapping information */
BmassQMType BmassQM;
/** \brief Global struct holding mortality output flags. */
MortFlagsType  MortFlags;

Bool UseGrid;
Bool UseProgressBar;


/**
 * \brief If TRUE an SQL database wil be generated containing debug information.
 * This flag can be set by adding "-STDebug" as an option when calling the program.
 * \ingroup SQL
 */
Bool STdebug_requested;
GrpIndex rg;
SppIndex sp;
IndivType *ndv; /* shorthand for the current indiv */


#ifndef STDEBUG

/** \brief Runs the program.
 *
 * Initializes flags and parameters, and runs the non-gridded mode. If the
 * user requests gridded mode this function calls RunGrid.
 */
int main(int argc, char **argv) {
  IntS year, iter;
	Bool killedany;

	sw_init_logs(stdout, &LogInfo);
	atexit(check_log);
	sw_init_logs(stdout, &LogInfoSW);
	/* provides a way to inform user that something
	 * was logged.  see generic.h */
    
    SW_DOM_init_ptrs(&SoilWatDomain);
    SW_CTL_init_ptrs(&SoilWatRun);

  SuperGlobals.prepare_IterationSummary = FALSE; // dont want to get soilwat output unless -o flag
  SuperGlobals.storeAllIterations = FALSE; // dont want to store all soilwat output iterations unless -i flag
  STdebug_requested = FALSE;

	init_args(argc, argv); // read input arguments and intialize proper flags

	printf("STEPWAT  init_args() executed successfully \n");

	if (UseGrid) {
        writeSOILWAT2Output = SuperGlobals.prepare_IterationSummary;
		runGrid();
		sw_write_warnings("", &LogInfo);
		sw_write_warnings("(SW) ", &LogInfoSW);
		return 0;
	}

    /* Read files.in and maxrgroupspecies.in before everything else */
    files_init();
    maxrgroupspecies_init();

	allocate_Globals();

	parm_Initialize();
        
	SXW_Init(TRUE, NULL); // allocate SOILWAT2-memory
	SW_OUT_set_ncol(SoilWatDomain.nMaxSoilLayers, SoilWatDomain.nMaxEvapLayers,
                    SoilWatRun.VegEstabIn.count, SoilWatDomain.OutDom.ncol_OUT,
                    SoilWatDomain.OutDom.nvar_OUT, SoilWatDomain.OutDom.nsl_OUT,
                    SoilWatDomain.OutDom.npft_OUT, &LogInfo); // set number of output columns
 	SW_OUT_set_colnames(SoilWatRun.RunIn.SiteRunIn.n_layers, SoilWatRun.VegEstabIn.parms,
 						SoilWatDomain.OutDom.ncol_OUT,
 						SoilWatDomain.OutDom.colnames_OUT, &LogInfo); // set column names for output files

 	if (SuperGlobals.prepare_IterationSummary) {
 		SW_OUT_create_summary_files(&SoilWatDomain.OutDom, &SoilWatRun.SW_PathOutputs,
                                    SoilWatDomain.SW_PathInputs.txtInFiles,
                                    SoilWatRun.RunIn.SiteRunIn.n_layers, &LogInfo);
        SW_OUT_construct_outarray(1, &SoilWatDomain.OutDom, &SoilWatRun.OutRun, &LogInfo);
 	}
        
	/* Connect to ST db and insert static data */
	if(STdebug_requested){
	    ST_connect("Output/stdebug");
	}

	/* --- Begin a new iteration ------ */
	for (iter = 1; iter <= SuperGlobals.runModelIterations; iter++) {
		Plot_Initialize();

		Globals->currIter = SoilWatRun.OutRun.currIter = iter;

		if (SuperGlobals.storeAllIterations) {
 			SW_OUT_create_iteration_files(&SoilWatDomain.OutDom, &SoilWatRun.SW_PathOutputs,
                                          iter, SoilWatDomain.SW_PathInputs.txtInFiles,
                                          SoilWatRun.RunIn.SiteRunIn.n_layers, &LogInfo);
		}

		if (SuperGlobals.prepare_IterationSummary) {
 			SoilWatDomain.OutDom.print_IterationSummary =
				(Bool) (Globals->currIter == SuperGlobals.runModelIterations);
		}

		/* ------  Begin running the model ------ */
		for (year = 1; year <= SuperGlobals.runModelYears; year++) {
            if(UseProgressBar){
                logProgress(iter, year, SIMULATION);
            }

			//printf("------------------------Repetition/year = %d / %d\n", iter, year);

			set_all_rngs(SuperGlobals.randseed, iter, year, 0);

			Globals->currYear = year;

			rgroup_Establish();

			Env_Generate();

			rgroup_PartResources();

			if (!isnull(SXW->debugfile) ) SXW_PrintDebug(0);

			rgroup_Grow();

			mort_Main(&killedany);

			rgroup_IncrAges();

			// Added functions for Grazing and mort_end_year as proportional killing effect before exporting biomass end of the year
			grazing_EndOfYear();

			save_annual_species_relsize();

      		mort_EndOfYear();

			stat_Collect(year);

			if (BmassFlags.yearly)
				output_Bmass_Yearly(year);

			// Moved kill annual and kill extra growth after we export biomass, and recovery of biomass after fire before the next year
			killAnnuals();
                        
			killMaxage();

			proportion_Recovery();

			killExtraGrowth();
			
			// if the user requests the stdebug.sqlite3 file to be generated
			// it is populated here.
			if(STdebug_requested){
				//species info
				ForEachSpecies(sp) {
					//individual info
					insertSpecieYearInfo(sp);
					for ((ndv) = Species[sp]->IndvHead; (ndv) != NULL; (ndv) = (ndv)->Next) {
						insertIndivYearInfo(ndv);
						insertIndiv(ndv);
					}
				}
				//Rgroup info
				ForEachGroup(rg){
					insertRGroupYearInfo(rg);
				}
			}
		} /* end model run for this year*/

		if (MortFlags.summary) {
			stat_Collect_GMort();
			stat_Collect_SMort();
		}

		if (MortFlags.yearly)
			output_Mort_Yearly(); // writes yearly file

		// dont need to restart if last iteration finished
		// this keeps it from re-writing the output folder and overwriting output files
		if (Globals->currIter != SuperGlobals.runModelIterations)
		{
			// don't reset in last iteration because we need to close files
			// before clearing/de-allocated SOILWAT2-memory
			SXW_Reset(SXW->f_watin, FALSE);
		}
	} /* end model run for this iteration*/

    if(UseProgressBar){
        logProgress(0, 0, OUTPUT);
    }

	/*------------------------------------------------------*/
	if (MortFlags.summary)
		stat_Output_AllMorts();
	if (BmassFlags.summary)
		stat_Output_AllBmass();
        
    /* Disconnect from the database */
	if(STdebug_requested){
		ST_disconnect();
	}

  if (!isnull(SXW->debugfile)){
    printf("entering debugfile\n");
    SXW_PrintDebug(1);
  }

  SW_OUT_close_files(&SoilWatRun.SW_PathOutputs, &SoilWatDomain.OutDom, &LogInfo);
  SW_DOM_deconstruct(&SoilWatDomain);
  SW_CTL_clear_model(TRUE, &SoilWatRun); // de-allocate all memory
  free_all_sxw_memory();
  freeMortalityMemory();

	deallocate_Globals(FALSE);

    sw_write_warnings("", &LogInfo);
    sw_write_warnings("(SW) ", &LogInfoSW);

    // This isn't wrapped in an if statement on purpose.
    // We should print "Done" either way.
    logProgress(0, 0, DONE);

	return 0;
}
/* END PROGRAM */
#endif

/** \brief (re)initializes the plot.
 * 
 * Zeros out Species and RGroup and kills all individuals.
 * Finally this function resets sxw.
 * 
 * \sa SXW_InitPlot(void)
 */
void Plot_Initialize(void) {
	GrpIndex rg;
	SppIndex sp;

	/* Clear remaining individuals and
	 resource counters */
	ForEachSpecies(sp)
	{

		if (!Species[sp]->use_me)
			continue;

		/* reset extirpated RGroups' species, if any */
		if (RGroup[Species[sp]->res_grp]->extirpated) {
			Species[sp]->seedling_estab_prob =
					Species[sp]->seedling_estab_prob_old;
		}

		/* clear estab and kills information */
		if (!isnull(Species[sp]->kills))
			Mem_Set(Species[sp]->kills, 0, sizeof(IntUS) * (SppMaxAge(sp)));

		/* Kill all individuals of each species.
		 This should zero everything necessary (inc. estab&kilz) */
		Species_Kill(sp, 0);

		/* This should no longer occur following the resolution of issue #209 on GitHub */
		if (!ZRO(getSpeciesRelsize(sp))) {
			LogError(&LogInfo, LOGWARN, 
							 "%s relsize = %f in Plot_Initialize. This indicates that some individuals in this species were not killed.",
							 Species[sp]->name, getSpeciesRelsize(sp));
		}
		if (Species[sp]->est_count) {
			LogError(&LogInfo, LOGWARN, "%s est_count (%d) forced "
					"in Plot_Initialize", Species[sp]->name,
					Species[sp]->est_count);
			Species[sp]->est_count = 0;
		}
	}

	ForEachGroup(rg)
	{

		if (!RGroup[rg]->use_me)
			continue;

		/* Clearing kills-accounting for survival data */
		if (!isnull(RGroup[rg]->kills))
			Mem_Set(RGroup[rg]->kills, 0, sizeof(IntUS) * GrpMaxAge(rg));

    /* This should no longer occur following the resolution of issue #209 on GitHub */
		if (!ZRO(getRGroupRelsize(rg))) {
			LogError(&LogInfo, LOGWARN, 
							 "%s relsize = %f in Plot_Initialize. This indicates that some individuals in this RGroup were not killed.",
							 RGroup[rg]->name, getRGroupRelsize(rg));
			/*printf("in plot_initialize before forcing, Rgroup = %s, relsize = %f, est_count= %d\n",
			RGroup[rg]->name, RGroup[rg]->relsize, RGroup[rg]->est_count); */
		}
                
		/* THIS NEVER SEEMS TO OCCUR */
		if (RGroup[rg]->est_count) {
			LogError(&LogInfo, LOGWARN, "%s est_count (%d) forced "
					"in Plot_Initialize", RGroup[rg]->name,
					RGroup[rg]->est_count);
			RGroup[rg]->est_count = 0;
		}
		RGroup[rg]->yrs_neg_pr = 0;
		RGroup[rg]->extirpated = FALSE;
	}

	/* Note, currently initCheatgrassPrecip is not being used, but retained for the time being */
	initCheatgrassPrecip();
	initWildfireClimate();
	SXW_InitPlot();
}


/** Set up all random number generators

  STEPWAT2 expects random number generators to produce sequences of random
  numbers that are reproducible (if user-provided "seed" is non-zero) and
      * unique among RNGs, iterations, years, and grid cells (most RNGs)
      * unique among RNGs, iterations, and years but identical among grid cells
        (weather generator RNG).

  A user-provided "seed" of zero produces non-reproducible random number
  sequences which are non-coinciding among RNGs, iterations, and grid cells.

  The set up of RNGs retains its characteristics even if `iter`, `year`, and/or
  `cell_id` are zero.

  See #RNG_INITSEQ and SOILWAT2's `RandSeed()` for further details.

  \param initstate The initial state of the system,
         i.e., the user provided "seed".
  \param iter The iteration identification number.
  \param year The simulated year.
  \param cell_id The cell identification number.

  \sideeffect STEPWAT2's random number generators are
              initialized with state and sequence.
 */
void set_all_rngs(
	unsigned long initstate,
	int iter,
	int year,
	int cell_id
) {
	/* Set up RNGs with seed/state and sequence identifier that is
		 reproducible and unique among RNGs, iterations, years and grid cells */
	RandSeed(initstate, RNG_INITSEQ(1, iter, year, cell_id), &environs_rng);
	RandSeed(initstate, RNG_INITSEQ(2, iter, year, cell_id), &mortality_rng);
	RandSeed(initstate, RNG_INITSEQ(3, iter, year, cell_id), &resgroups_rng);
	RandSeed(initstate, RNG_INITSEQ(4, iter, year, cell_id), &species_rng);
	RandSeed(initstate, RNG_INITSEQ(5, iter, year, cell_id), &grid_rng);
	RandSeed(initstate, RNG_INITSEQ(6, iter, year, cell_id), &dispersal_rng);
	RandSeed(initstate, RNG_INITSEQ(7, iter, year, cell_id), &resource_rng);

	/* Initialize RNGs with seed/state and sequence identifier that is
		 reproducible and unique among RNGs, iterations, and year
		 but not grid cells */
	RandSeed(initstate, RNG_INITSEQ(8, iter, year, 0), &SoilWatRun.MarkovIn.markov_rng);
}



#ifndef STDEBUG

/**************************************************************/
/* Allocates memory for any global variables defined as pointers */
void allocate_Globals(void){
	Env = (EnvType*) Mem_Calloc(1, sizeof(EnvType), "allocate_Globals: Env", &LogInfo);
	Succulent = (SucculentType*) Mem_Calloc(1, sizeof(SucculentType), "allocate_Globals: Succulent", &LogInfo);
	Globals = (ModelType*) Mem_Calloc(1, sizeof(ModelType), "allocate_Globals: Globals", &LogInfo);
	Plot = (PlotType*) Mem_Calloc(1, sizeof(PlotType), "allocate_Globals: Plot", &LogInfo);
	_SomeKillage = (Bool*) Mem_Calloc(1, sizeof(Bool), "allocate_Globals: _SomeKillage", &LogInfo);
}

/* Deallocates the global variables */
void deallocate_Globals(Bool isGriddedMode){
	GrpIndex rg;
	SppIndex sp;
	
	double **bMassQMFreeArray[] = {
		&BmassQM.rap_annual_points,
		&BmassQM.rap_perennial_points,
		&BmassQM.stepwat_annual_points,
		&BmassQM.stepwat_perennial_points
	};
	const int numBmassFreeElems = 4;
	int bMassIndex;

	/* Free Species */
	ForEachSpecies(sp){
		/* Start by freeing any pointers in the Species struct */
		free(Species[sp]->kills);
		free(Species[sp]->seedprod);
		free(Species[sp]->name);
		IndivType *indv = Species[sp]->IndvHead, *next;
		/* Next free the linked list of individuals */
		while(indv){
			next = indv->Next;
			free(indv);
			indv = next;
		}
		free(indv);
		/* Finally free the actual species */
		free(Species[sp]);
	}
	/* Then free the entire array */
	free(Species); 

	/* Free RGroup */
	ForEachGroup(rg){
		/* Free all pointers in the RGroup struct */
		free(RGroup[rg]->est_spp);
		free(RGroup[rg]->kills);
		free(RGroup[rg]->name);
		free(RGroup[rg]->species);
		free(RGroup[rg]);
	}
	/* Then free the entire array */
	free(RGroup);

	if(!isGriddedMode){
			free(Env);
			free(Succulent);
			free(Globals);
			free(Plot);
			free(_SomeKillage);
		}
		
	/* Free BmassQM */
	for (bMassIndex = 0; bMassIndex < numBmassFreeElems; bMassIndex++) {
		if (!isnull(*bMassQMFreeArray[bMassIndex])) {
			free(*bMassQMFreeArray[bMassIndex]);
			*bMassQMFreeArray[bMassIndex] = NULL;
		}
	}
}

/** \brief Translates the input flags to in program flags.
 * 
 * The recognised flags are -d, -f, -q, -e, -p, -g, -o, -i, -s and -S.
 * Note that flags are case sensitive. 
 * 
 * When the -f flag is uses this function looks next for the name of the file.
 */
static void init_args(int argc, char **argv) {
  /* to add an option:
   *  - include it in opts[]
   *  - set a flag in valopts indicating no value (0),
   *    value required (1), or value optional (-1),
   *  - then tell us what to do in the switch statement
   *
   * 3/1/03 - cwb - Current options are
   *                -d=chg to work dir <opt=dir_name>
   *                -f=chg deflt first file <opt=file.in>
   *                -q=quiet, noprint "Check logfile" at end of program
   *                -e=echo init values to logfile.
   * 06/27/16 -AKT  -o= Print all the Soilwat output as well while running with STEPWAT
   *                  This option is required to have soilwat_input_files is next path after this
   *                  like  -o ../../sw_src/testing/files_step_soilwat_grid.in
   * 1/8/04 - cwb - Added -p option to help the GUI with a progress bar.
   *         This is another "secret" option, insofar as the
   *         command-line user doesn't need it.  The option directs
   *         the program to write progress info (iter) to stdout.
   *         Without the option, progress info (dots) is written to
   *         stderr.
   * 8/16/17 - BEB  Updated option for -o flag. Now if this flag is set the output files from
   *                files_v30.in are written to.
   * 10/9/17 - BEB Added -i flag for writing SOILWAT output for every iteration
   */
  char str[1024],
       *opts[]  = {"-d","-f","-q","-e", "-p", "-g", "-o", "-i", "-s", "-S"};  /* valid options */
  int valopts[] = {  1,   1,   0,  -1,   0,    0,    0,   0,   0,   0};  /* indicates options with values */
                 /* 0=none, 1=required, -1=optional */
  int i, /* looper through all cmdline arguments */
      a, /* current valid argument-value position */
      op, /* position number of found option */
      nopts=sizeof(opts)/sizeof(char *);
  Bool lastop_noval = FALSE;

  /* Defaults */
  parm_SetFirstName( DFLT_FIRSTFILE);
  LogInfo.QuietMode = EchoInits = UseSeedDispersal = FALSE;
  LogInfo.logfp = stderr;


  a=1;
	for (i = 1; i <= nopts; i++)
	{
		if (a >= argc)
			break;

		/* figure out which option by its position 0-(nopts-1) */
		for (op = 0; op < nopts; op++)
		{
			if (strncmp(opts[op], argv[a], 2) == 0)
				break; /* found it, move on */
		}
		if (op == nopts)
		{
			fprintf(stderr, "Invalid option %s\n", argv[a]);
			usage();
			exit(-1);
		}
		if (a == argc - 1 && strlen(argv[a]) == 2)
			lastop_noval = TRUE;

		*str = '\0';
		/* extract value part of option-value pair */
		if (valopts[op])
		{
			if (lastop_noval && valopts[op] < 0)
			{
				/* break out, optional value not available */
				/* avoid checking past end of array */

			}
			else if (lastop_noval && valopts[op] > 0)
			{
				fprintf(stderr, "Incomplete option %s\n", opts[op]);
				usage();
				exit(-1);

			}
			else if ('\0' == argv[a][2] && valopts[op] < 0)
			{
				/* break out, optional value not available */

			}
			else if ('\0' != argv[a][2])
			{ /* no space betw opt-value */
				strcpy(str, (argv[a] + 2));

			}
			else if ('-' != *argv[a + 1])
			{ /* space betw opt-value */
				strcpy(str, argv[++a]);

			}
			else if (0 < valopts[op])
			{ /* required opt-val not found */
				fprintf(stderr, "Incomplete option %s\n", opts[op]);
				usage();
				exit(-1);
			} /* opt-val not required */
		}

		/* set indicators/variables based on results */
		switch (op)
		{
		case 0: /* -d */
			if (!ChDir(str))
			{
				LogError(&LogInfo, LOGERROR, "Invalid project directory (%s)",
						str);
			}
			break;
		case 1:
			parm_SetFirstName(str);
			break; /* -f */

		case 2:
			LogInfo.QuietMode = TRUE;
			break; /* -q */

		case 3:
			EchoInits = TRUE;
			break; /* -e */

		case 4:
			LogInfo.logfp = stdout; /* -p */
			UseProgressBar = TRUE;
			break;

		case 5:
			UseGrid = TRUE;
			break; /* -g */

		case 6:
      		printf("storing SOILWAT output aggregated across-iterations (-o flag)\n");
      		SuperGlobals.prepare_IterationSummary = TRUE;
			break; /* -o */

    	case 7: // -i
      		printf("storing SOILWAT output for each iteration (-i flag)\n");
      		SuperGlobals.storeAllIterations = TRUE;
      		break;

		case 8: // -s
			if (strlen(argv[a]) > 1){
				printf("Generating SXW debug file\n");
				SXW->debugfile = Str_Dup(&argv[a][1], &LogInfo);
			}
			break;
	  
	  	case 9: // -S
		    if(!strncmp("-STdebug", argv[a], 8)){ //  -STdebug
				printf("Generating STdebug.sqlite database (-STdebug flag)\n");
				STdebug_requested = TRUE;
				break;
			} else {
				printf("Invalid option. (Did you mean -STdebug?)\n");
				usage();
			}

		default:
			LogError(&LogInfo, LOGERROR,
					"Programmer: bad option in main:init_args:switch");
		}

		a++; /* move to next valid option-value position */

	} /* end for(i) */

  LogInfoSW.logfp = LogInfo.logfp;

}

/** \brief Prints a warning if there is an entry in the logfile
 * 
 * The warning is printed to the progress file, which is usually the 
 * same as stdout.
 * 
 * check_log is registered to run automatically at exit.
 */
static void check_log(void) {
/* =================================================== */

  if (LogInfo.logfp != stdout) {
    if ((LogInfo.stopRun || LogInfo.numWarnings > 0) && !LogInfo.QuietMode)
      fprintf(LogInfo.logfp, "\nCheck logfile for error messages.\n");

    CloseFile(&LogInfo.logfp, &LogInfo);
  }

}
#endif

/** \brief Compares the getRelsize funcitons to calculated values.
 * 
 * Used for debugging. 
 * 
 * \sa getSpeciesRelsize
 * \sa getRGroupRelsize
 */
void check_sizes(const char *chkpt) {
    /* =================================================== */
    /* Use this for debugging to check that the sum of the individual
     * sizes add up to the RGroup and Species relsize registers.
     * The chkpt is a string that gets output to help you know where
     * the difference was found.
     */
    GrpIndex rg;
    SppIndex sp;
    IndivType *ndv;
    int i;
    RealF spsize, rgsize,
            diff = .000005; /* amount of difference allowed */

    ForEachGroup(rg) {

        rgsize = 0.0;

        ForEachEstSpp(sp, rg, i) {

            spsize = 0.0;
            ForEachIndiv(ndv, Species[sp]) spsize += ndv->relsize;
            rgsize += spsize;

            if (LT(diff, fabs(spsize - getSpeciesRelsize(sp)))) {
                LogError(&LogInfo, LOGWARN, "%s (%d:%d): SP: \"%s\" size error: "
                        "SP=%.7f, ndv=%.7f",
                        chkpt, Globals->currIter, Globals->currYear,
                        Species[sp]->name, getSpeciesRelsize(sp), spsize);
            }
        }

        if (LT(diff, fabs(rgsize - getRGroupRelsize(rg)))) {
            LogError(&LogInfo, LOGWARN, "%s (%d:%d): RG \"%s\" size error: "
                    "RG=%.7f, ndv=%.7f",
                    chkpt, Globals->currIter, Globals->currYear,
                    RGroup[rg]->name, getRGroupRelsize(rg), rgsize);
        }
    }

}

#ifdef DEBUG_GROW
/**************************************************************/
void Debug_AddByIter( Int iter) {

     Species_Add_Indiv(1,1);
   Species_Add_Indiv(2,4);
   Species_Add_Indiv(12,4);

}

void Debug_AddByYear( Int year) {

  if (year == 1 && 1) {
    if(1) Species_Add_Indiv(1,3);
    if(1) Species_Add_Indiv(2,10);
    if(1) Species_Add_Indiv(12,20);
  }

  if (year == 20 && 1) {
    if(1) Species_Add_Indiv(1,1);
  }
  if (year > 3 && !(year & 5) && 0) {
    if(1) Species_Add_Indiv(1,1);
    if(1) Species_Add_Indiv(2,1);
    if(1) Species_Add_Indiv(12,2);
  }

  if (year == 140 && 1) {
    Species_Add_Indiv(2,1);
  }

}

#endif
