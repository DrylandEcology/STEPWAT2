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
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sw_src/SW_Control.h"
#include "sxw_funcs.h"
#include "sw_src/SW_Output.h"
#include "sw_src/SW_Output_outtext.h"
#include "sw_src/SW_Output_outarray.h"
#include "sw_src/rands.h"
#include "ST_spinup.h"
#include "ST_progressBar.h"
#include "ST_mortality.h"

extern Bool prepare_IterationSummary; // defined in `SOILWAT2/SW_Output.c`
extern Bool print_IterationSummary; // defined in `SOILWAT2/SW_Output_outtext.c`
extern Bool storeAllIterations; // defined in `SOILWAT2/SW_Output.c`
extern SW_VEGPROD SW_VegProd;
extern Bool *_SomeKillage;				// From ST_mortality.c 
SW_FILE_STATUS SW_File_Status;

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

#ifdef DEBUG_MEM
  #define chkmem_f CheckMemoryIntegrity(FALSE);
  #define chkmem_t CheckMemoryIntegrity(TRUE);
  void CheckMemoryIntegrity(Bool flag); /* local */
  void Stat_NoteMemoryRefs(void) ;
#else
  #define chkmem_f
  #define chkmem_t
#endif

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

Bool QuietMode;

/************ External Variable Definitions  ***************/
/*              see ST_globals.h                       */
/***********************************************************/
char errstr[1024];
char inbuf[1024];
/** \brief The log file */
FILE *logfp;
/** \brief optional place to put progress info */
FILE *progfp;
/** \brief indicator that err file was written to */
int logged;
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

BmassFlagsType BmassFlags;
/** \brief Global struct holding mortality output flags. */
MortFlagsType  MortFlags;

Bool UseGrid;
Bool EchoInits;
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

pcg32_random_t environs_rng;
pcg32_random_t resgroups_rng;
pcg32_random_t species_rng;
pcg32_random_t grid_rng;
extern pcg32_random_t markov_rng;

#ifndef STDEBUG

/** \brief Runs the program.
 * 
 * Initializes flags and parameters, and runs the non-gridded mode. If the
 * user requests gridded mode this function calls RunGrid.
 */
int main(int argc, char **argv) {
  IntS year, iter;
	Bool killedany;

	logged = FALSE;
	atexit(check_log);
	/* provides a way to inform user that something
	 * was logged.  see generic.h */

  prepare_IterationSummary = FALSE; // dont want to get soilwat output unless -o flag
  storeAllIterations = FALSE; // dont want to store all soilwat output iterations unless -i flag
  STdebug_requested = FALSE;

	init_args(argc, argv); // read input arguments and intialize proper flags

	printf("STEPWAT  init_args() executed successfully \n");

	if (UseGrid) {
        writeSOILWAT2Output = prepare_IterationSummary;
		runGrid();
		return 0;
	}

    /* Read files.in and maxrgroupspecies.in before everything else */
    files_init();
    maxrgroupspecies_init();

	allocate_Globals();

	parm_Initialize();
        
	SXW_Init(TRUE, NULL); // allocate SOILWAT2-memory
	SW_OUT_set_ncol(); // set number of output columns
	SW_OUT_set_colnames(); // set column names for output files
	if (prepare_IterationSummary) {
		SW_OUT_create_summary_files();
		// allocate `p_OUT` and `p_OUTsd` arrays to aggregate SOILWAT2 output across iterations
		setGlobalSTEPWAT2_OutputVariables();
	}
        
	/* Connect to ST db and insert static data */
	if(STdebug_requested){
	    ST_connect("Output/stdebug");
	}

	/* --- Begin a new iteration ------ */
	for (iter = 1; iter <= SuperGlobals.runModelIterations; iter++) {
		Plot_Initialize();

		RandSeed(SuperGlobals.randseed, &environs_rng);
		RandSeed(SuperGlobals.randseed, &mortality_rng);
		RandSeed(SuperGlobals.randseed, &resgroups_rng);
		RandSeed(SuperGlobals.randseed, &species_rng);
		RandSeed(SuperGlobals.randseed, &grid_rng);
		RandSeed(SuperGlobals.randseed, &markov_rng);

		Globals->currIter = iter;
                
		if (storeAllIterations) {
			SW_OUT_create_iteration_files(Globals->currIter);
		}

		if (prepare_IterationSummary) {
			print_IterationSummary = (Bool) (Globals->currIter == SuperGlobals.runModelIterations);
		}

		/* ------  Begin running the model ------ */
		for (year = 1; year <= SuperGlobals.runModelYears; year++) {
            if(UseProgressBar){
                logProgress(iter, year, SIMULATION);
            }

			//printf("------------------------Repetition/year = %d / %d\n", iter, year);

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
			SXW_Reset(SXW->f_watin);
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

  SW_OUT_close_files();
  SW_CTL_clear_model(TRUE); // de-allocate all memory
  free_all_sxw_memory();
  freeMortalityMemory();

	deallocate_Globals(FALSE);

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
			LogError(logfp, LOGNOTE, 
							 "%s relsize = %f in Plot_Initialize. This indicates that some individuals in this species were not killed.",
							 Species[sp]->name, getSpeciesRelsize(sp));
		}
		if (Species[sp]->est_count) {
			LogError(logfp, LOGNOTE, "%s est_count (%d) forced "
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
			LogError(logfp, LOGNOTE, 
							 "%s relsize = %f in Plot_Initialize. This indicates that some individuals in this RGroup were not killed.",
							 RGroup[rg]->name, getRGroupRelsize(rg));
                        /*printf("in plot_initialize before forcing, Rgroup = %s, relsize = %f, est_count= %d\n",
                        RGroup[rg]->name, RGroup[rg]->relsize, RGroup[rg]->est_count); */
		}
                
		/* THIS NEVER SEEMS TO OCCUR */
		if (RGroup[rg]->est_count) {
			LogError(logfp, LOGNOTE, "%s est_count (%d) forced "
					"in Plot_Initialize", RGroup[rg]->name,
					RGroup[rg]->est_count);
			RGroup[rg]->est_count = 0;
		}
		RGroup[rg]->yrs_neg_pr = 0;
		RGroup[rg]->extirpated = FALSE;
	}

	initCheatgrassPrecip();
	SXW_InitPlot();
}

#ifndef STDEBUG

/**************************************************************/
/* Allocates memory for any global variables defined as pointers */
void allocate_Globals(void){
	Env = (EnvType*) Mem_Calloc(1, sizeof(EnvType), "allocate_Globals: Env");
	Succulent = (SucculentType*) Mem_Calloc(1, sizeof(SucculentType), "allocate_Globals: Succulent");
	Globals = (ModelType*) Mem_Calloc(1, sizeof(ModelType), "allocate_Globals: Globals");
	Plot = (PlotType*) Mem_Calloc(1, sizeof(PlotType), "allocate_Globals: Plot");
	_SomeKillage = (Bool*) Mem_Calloc(1, sizeof(Bool), "allocate_Globals: _SomeKillage");
}

/* Deallocates the global variables */
void deallocate_Globals(Bool isGriddedMode){
	GrpIndex rg;
	SppIndex sp;

	if(!isGriddedMode){
		Mem_Free(Env);
		Mem_Free(Succulent);
		Mem_Free(Globals);
		Mem_Free(Plot);
		Mem_Free(_SomeKillage);
	}
	
	/* Free Species */
	ForEachSpecies(sp){
		/* Start by freeing any pointers in the Species struct */
		Mem_Free(Species[sp]->kills);
		Mem_Free(Species[sp]->seedprod);
		Mem_Free(Species[sp]->name);
		IndivType *indv = Species[sp]->IndvHead, *next;
		/* Next free the linked list of individuals */
		while(indv){
			next = indv->Next;
			Mem_Free(indv);
			indv = next;
		}
		Mem_Free(indv);
		/* Finally free the actual species */
		Mem_Free(Species[sp]);
	}
	/* Then free the entire array */
	Mem_Free(Species); 

	/* Free RGroup */
	ForEachGroup(rg){
		/* Free all pointers in the RGroup struct */
		Mem_Free(RGroup[rg]->est_spp);
		Mem_Free(RGroup[rg]->kills);
		Mem_Free(RGroup[rg]->name);
		Mem_Free(RGroup[rg]->species);
		Mem_Free(RGroup[rg]);
	}
	/* Then free the entire array */
	Mem_Free(RGroup);
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
  QuietMode = EchoInits = UseSeedDispersal = FALSE;
  progfp = stderr;


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
				LogError(stderr, LOGFATAL, "Invalid project directory (%s)",
						str);
			}
			break;
		case 1:
			parm_SetFirstName(str);
			break; /* -f */

		case 2:
			QuietMode = TRUE;
			break; /* -q */

		case 3:
			EchoInits = TRUE;
			break; /* -e */

		case 4:
			progfp = stdout; /* -p */
			UseProgressBar = TRUE;
			break;

		case 5:
			UseGrid = TRUE;
			break; /* -g */

		case 6:
      		printf("storing SOILWAT output aggregated across-iterations (-o flag)\n");
      		prepare_IterationSummary = TRUE;
			break; /* -o */

    	case 7: // -i
      		printf("storing SOILWAT output for each iteration (-i flag)\n");
      		storeAllIterations = TRUE;
      		break;

		case 8: // -s
			if (strlen(argv[a]) > 1){
				printf("Generating SXW debug file\n");
				SXW->debugfile = Str_Dup(&argv[a][1]);
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
			LogError(logfp, LOGFATAL,
					"Programmer: bad option in main:init_args:switch");
		}

		a++; /* move to next valid option-value position */

	} /* end for(i) */


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

  if (logfp != stdout) {
    if (logged && !QuietMode)
      fprintf(progfp, "\nCheck logfile for error messages.\n");

    CloseFile(&logfp);
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
                LogError(stdout, LOGWARN, "%s (%d:%d): SP: \"%s\" size error: "
                        "SP=%.7f, ndv=%.7f",
                        chkpt, Globals->currIter, Globals->currYear,
                        Species[sp]->name, getSpeciesRelsize(sp), spsize);
            }
        }

        if (LT(diff, fabs(rgsize - getRGroupRelsize(rg)))) {
            LogError(stdout, LOGWARN, "%s (%d:%d): RG \"%s\" size error: "
                    "RG=%.7f, ndv=%.7f",
                    chkpt, Globals->currIter, Globals->currYear,
                    RGroup[rg]->name, getRGroupRelsize(rg), rgsize);
        }
    }

}

#ifdef DEBUG_MEM
void CheckMemoryIntegrity(Bool flag) {
// DLM - 6/6/2013 : NOTE - The dynamically allocated variables added for the new grid option are not accounted for here.  I haven't bothered to add them because it would add a ton of code and I haven't been using the code to facilitate debugging of memory.  Instead I have been using the valgrind program to debug my memory, as I feel like it is a better solution.  If wanting to use this code to debug memory, memory references would probably have to be set up for the dynamically allocated variables in ST_grid.c as well as the new dynamically allocated grid_Stat variable in ST_stats.c.


  ClearMemoryRefs();

  if (flag || Globals.currIter > 1 || Globals.currYear > 1)
    Stat_SetMemoryRefs();

  RGroup_SetMemoryRefs();
  Species_SetMemoryRefs();
  Parm_SetMemoryRefs();

  SXW_SetMemoryRefs();
  CheckMemoryRefs();

}
#endif



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
