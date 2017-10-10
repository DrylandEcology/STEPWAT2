/********************************************************/
/********************************************************/
/*  Source file: ST_main.c
/*  Type: module
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: Main program loop and initializations.
/*  History:
/*     (6/15/2000) -- INITIAL CODING - cwb
 *     15-Apr-02 (cwb) -- added code to interface with SOILWAT
 *	   5-24-2013 (DLM) -- added gridded option to program... see ST_grid.c source file for the rest of the gridded code
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ST_steppe.h"
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"


#ifdef STEPWAT
  #include "sxw_funcs.h"
  #include "sxw.h"
  extern SXW_t SXW;
#endif

extern Bool isPartialSoilwatOutput;
extern Bool storeAllIterations;

/************* External Function Declarations **************/
/***********************************************************/

  void rgroup_Grow( void);
  void rgroup_Establish( void) ;
  void rgroup_IncrAges( void);
  void rgroup_PartResources( void);

  void mort_Main( Bool *killed);
  void mort_EndOfYear( void);
  void proportion_Recovery(void);
  void grazing_EndOfYear( void);

  void parm_Initialize( Int);
  void parm_SetFirstName( char *s);

  void output_Bmass_Yearly( Int year );
  void output_Mort_Yearly( void );

  void stat_Collect( Int year ) ;
  void stat_Collect_GMort ( void ) ;
  void stat_Collect_SMort ( void ) ;

  void stat_Output_AllMorts( void) ;
  void stat_Output_AllBmass(void) ;

  void stat_Output_AllSoilwatVariables(void);

  void runGrid( void ); //for the grid... declared in ST_grid.c

  void _kill_annuals(void);
  void _kill_extra_growth(void);

#ifdef DEBUG_MEM
  #define chkmem_f CheckMemoryIntegrity(FALSE);
  #define chkmem_t CheckMemoryIntegrity(TRUE);
  void CheckMemoryIntegrity(Bool flag); /* local */
  void Stat_NoteMemoryRefs(void) ;
#else
  #define chkmem_f
  #define chkmem_t
#endif


#ifndef STEPWAT
 void SXW_Init( Bool init_SW, char *f_roots ) {}
 void SXW_Run_SOILWAT(void) {}
 void SXW_InitPlot( void ) {}
#endif

/*************** Local Function Declarations ***************/
/***********************************************************/
void Plot_Initialize( void);

static void usage(void) {
  char *s ="STEPPE plant community dynamics (SGS-LTER Jan-04).\n"
           "Usage: steppe [-d startdir] [-f files.in] [-q] [-s] [-e] [-o] [-g]\n"
           "  -d : supply working directory (default=.)\n"
           "  -f : supply list of input files (default=files.in)\n"
           "  -q : quiet mode, don't print message to check logfile.\n"
           "  -s : use SOILWAT model for resource partitioning.\n"
           "  -e : echo initialization results to logfile\n"
           "  -o : print all the soilwat output\n"
           "  -g : use gridded mode\n"
           "  -i : print SOILWAT output for each iteration\n"; // dont need to set -o flag to use this flag
  fprintf(stderr,"%s", s);
  exit(0);
}
static void init_args(int argc, char **argv);
static void check_log(void);

/* a couple of debugging routines */
void check_sizes(const char *);

Bool QuietMode;

/************ External Variable Definitions  ***************/
/*              see ST_globals.h                       */
/***********************************************************/
char errstr[1024];
char inbuf[1024];
FILE *logfp,   /* used everywhere by LogError */
     *progfp;  /* optional place to put progress info */
int logged;  /* indicator that err file was written to */
SpeciesType   *Species[MAX_SPECIES];
GroupType     *RGroup [MAX_RGROUPS];
SucculentType  Succulent;
EnvType        Env;
PlotType       Plot;
ModelType      Globals;
BmassFlagsType BmassFlags;
MortFlagsType  MortFlags;

Bool UseSoilwat; /* not used in every module */
Bool UseGrid;
Bool UseSeedDispersal;
Bool DuringSpinup;
Bool EchoInits;
Bool UseProgressBar;

/******************** Begin Model Code *********************/
/***********************************************************/
int main(int argc, char **argv) {
  IntS year, iter, incr, g, s;
	IndivType *i;
	Bool killedany;
  int k;

	logged = FALSE;
	atexit(check_log);
	/* provides a way to inform user that something
	 * was logged.  see generic.h */

  isPartialSoilwatOutput = TRUE; // dont want to get soilwat output unless -o flag
  storeAllIterations = FALSE; // dont want to store all soilwat output iterations unless -i flag

	init_args(argc, argv); // read input arguments and intialize proper flags

	printf("STEPWAT  init_args() executed successfully \n ");

	if (UseGrid == TRUE) {
		runGrid();
		return 0;
	}

	parm_Initialize(0);

	if (UseSoilwat)
		SXW_Init(TRUE, NULL);

	incr = (IntS) ((float) Globals.runModelIterations / 10);
	if (incr == 0)
		incr = 1;

  /*----------------------------------------------------------
    Get proper order for rank_SWPcrits
  ----------------------------------------------------------*/
  int outerLoop, innerLoop;
  float key;
  RealF tempArray[4];
  tempArray[0] = SXW.critSoilWater[0];
  tempArray[1] = SXW.critSoilWater[1];
  tempArray[2] = SXW.critSoilWater[2];
  tempArray[3] = SXW.critSoilWater[3];

  // printing for testing purposes
  /*printf("%f\n", tempArray[0]);
  printf("%f\n", tempArray[1]);
  printf("%f\n", tempArray[2]);
  printf("%f\n\n", tempArray[3]);*/

  // insertion sort to rank the veg types and store them in their proper order
  for (outerLoop = 1; outerLoop < 4; outerLoop++)
   {
       key = tempArray[outerLoop];
       innerLoop = outerLoop-1;
       while (innerLoop >= 0 && tempArray[innerLoop] < key)
       {
         // code to switch values
         SXW.rank_SWPcrits[innerLoop+2] = innerLoop;
         innerLoop = innerLoop-1;
       }
       tempArray[innerLoop+1] = key;
       SXW.rank_SWPcrits[innerLoop+2] = outerLoop;
   }
   SXW.rank_SWPcrits[0] = 0; // setting to 0 since that is was equation calls for

   // printing for testing purposes
   /*printf("%d = %f\n", SXW.rank_SWPcrits[1], tempArray[0]);
   printf("%d = %f\n", SXW.rank_SWPcrits[2], tempArray[1]);
   printf("%d = %f\n", SXW.rank_SWPcrits[3], tempArray[2]);
   printf("%d = %f\n\n", SXW.rank_SWPcrits[4], tempArray[3]);*/

   /*----------------------------------------------------------
     End of rank_SWPcrits
   ----------------------------------------------------------*/

	/* --- Begin a new iteration ------ */
	for (iter = 1; iter <= Globals.runModelIterations; iter++) {
		if (progfp == stderr) {
			if (iter % incr == 0)
				fprintf(progfp, ".");
		} else {
			fprintf(progfp, "%d\n", iter);
		}

    // set these to 0 for use with -i flag (need to create column headers for every iteration file)
    SXW.col_status_dy = 0;
    SXW.col_status_wk = 0;
    SXW.col_status_mo = 0;
    SXW.col_status_yr = 0;

		if (BmassFlags.yearly || MortFlags.yearly)
			parm_Initialize(iter);

		Plot_Initialize();
		RandSeed(Globals.randseed);
		Globals.currIter = iter;

    memset(SXW.transp_SWA,0,sizeof(SXW.transp_SWA)); // set transp_SWA to 0; needs to be reset each iteration

		/* ------  Begin running the model ------ */
		for (year = 1; year <= Globals.runModelYears; year++) {
			Globals.currYear = year;

			rgroup_Establish();

			Env_Generate();

			rgroup_PartResources();

#ifdef STEPWAT
			if (!isnull(SXW.debugfile) ) SXW_PrintDebug(0);
#endif

			rgroup_Grow();

			mort_Main(&killedany);

			rgroup_IncrAges();

      // Added functions for Grazing and mort_end_year as proportional killing effect before exporting biomass end of the year
			grazing_EndOfYear();

			mort_EndOfYear();

			stat_Collect(year);

			if (BmassFlags.yearly)
				output_Bmass_Yearly(year);

       // Moved kill annual and kill extra growth after we export biomass, and recovery of biomass after fire before the next year
			_kill_annuals();
			 proportion_Recovery();
			_kill_extra_growth();
		} /* end model run for this year*/

		if (MortFlags.summary) {
			stat_Collect_GMort();
			stat_Collect_SMort();
		}

		if (MortFlags.yearly)
			output_Mort_Yearly(); // writes yearly file

		if (UseSoilwat)
			{
		     	//stat_Output_AllSoilwatVariables();
          // dont need to restart if last iteration finished
          // this keeps it from re-writing the output folder and overwriting output files
          if(Globals.currIter != Globals.runModelIterations)
            SXW_Reset();
			}
	} /* end model run for this iteration*/

	/*------------------------------------------------------*/
	if (MortFlags.summary)
		stat_Output_AllMorts();
	if (BmassFlags.summary)
		stat_Output_AllBmass();

#ifdef STEPWAT
			if (!isnull(SXW.debugfile)){
        printf("entering debugfile\n");
        SXW_PrintDebug(1);
      }
#endif

  printf("\nend program\n");
	fprintf(progfp, "\n");
	return 0;
}
/* END PROGRAM */



/**************************************************************/
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

		/* programmer alert: INVESTIGATE WHY THIS OCCURS */
		if (!ZRO(Species[sp]->relsize)) {
			LogError(logfp, LOGNOTE, "%s relsize (%f) forced "
					"in Plot_Initialize", Species[sp]->name,
					Species[sp]->relsize);
			Species[sp]->relsize = 0.0;
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

	if (UseSoilwat)
		SXW_InitPlot();
}


/**************************************************************/
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
   *                -s=soilwat model-derived resource
   *                   optional parm debugfile for pgmr testing; see code.
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
       *opts[]  = {"-d","-f","-q","-s","-e", "-p", "-g", "-o", "-i"};  /* valid options */
  int valopts[] = {  1,   1,   0,  -1,   0,    0,    0,   0,   0};  /* indicates options with values */
                 /* 0=none, 1=required, -1=optional */
  int i, /* looper through all cmdline arguments */
      a, /* current valid argument-value position */
      op, /* position number of found option */
      nopts=sizeof(opts)/sizeof(char *);
  Bool lastop_noval = FALSE;

  /* Defaults */
  parm_SetFirstName( DFLT_FIRSTFILE);
  UseSoilwat = QuietMode = EchoInits = UseSeedDispersal = FALSE;
  SXW.debugfile = NULL;
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
			UseSoilwat = TRUE; /* -s */
			if (strlen(str) > 1)
				SXW.debugfile = Str_Dup(str);
			break;

		case 4:
			EchoInits = TRUE;
			break; /* -e */

		case 5:
			progfp = stdout; /* -p */
			UseProgressBar = TRUE;
			break;

		case 6:
			UseGrid = TRUE;
			break; /* -g */

		case 7:
      printf("storing SOILWAT output (flag -o)\n");
      isPartialSoilwatOutput = FALSE;
			break; /* -o    also get all the soilwat output*/

    case 8: // -i
      printf("storing SOILWAT output for all iterations\n");
      storeAllIterations = TRUE;
      break;

		default:
			LogError(logfp, LOGFATAL,
					"Programmer: bad option in main:init_args:switch");
		}

		a++; /* move to next valid option-value position */

	} /* end for(i) */


}


static void check_log(void) {
/* =================================================== */

  if (logfp != stdout) {
    if (logged && !QuietMode)
      fprintf(progfp, "\nCheck logfile for error messages.\n");

    CloseFile(&logfp);
  }

}


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
        diff=.000005;    /* amount of difference allowed */


  ForEachGroup(rg) {

    rgsize=0.0;
    ForEachEstSpp(sp, rg, i) {

      spsize = 0.0;
      ForEachIndiv(ndv, Species[sp]) spsize += ndv->relsize;
      rgsize += spsize;

      if (LT(diff, abs(spsize - Species[sp]->relsize)) ) {
        LogError(stdout, LOGWARN, "%s (%d:%d): SP: \"%s\" size error: "
                                  "SP=%.9f, ndv=%.9f\n",
                chkpt, Globals.currIter, Globals.currYear,
                Species[sp]->name, Species[sp]->relsize, spsize);
      }
    }

    if ( LT(diff, abs(rgsize -RGroup[rg]->relsize)) ) {
      LogError(stdout, LOGWARN, "%s (%d:%d): RG \"%s\" size error: "
                                "RG=%.9f, ndv=%.9f\n",
              chkpt, Globals.currIter, Globals.currYear,
              RGroup[rg]->name, RGroup[rg]->relsize, rgsize);
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

#ifdef STEPWAT
  SXW_SetMemoryRefs();
#endif

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
