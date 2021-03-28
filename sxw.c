/**
 *  \file sxw.c
 *  \brief Interface module for the STEPWAT2 to SOILWAT
 *         data flow.
 * 
 *  Application: STEPWAT2 - plant community dynamics simulator
 *  coupled with the  SOILWAT2 model.
 *  History 
 *     (9-May-2002) -- INITIAL CODING - cwb
 *     28-Feb-02 - cwb - The model runs but plants die
 *         soon after establishment in a way that suggests
 *         chronic stretching of resources.  At this time
 *         I'm setting the input to soilwat to be based on
 *         full-sized plants to provide maximum typical
 *         transpiration and interpreting that as "available"
 *         resource. Requires the addition of a new array to
 *         hold the summation of the maximum (mature) biomass
 *         of each group _Grp_BMass[].  The only affected routines
 *         in this file are SXW_Init() and SXW_Run_SOILWAT().
 *         But see also sxw_soilwat.c.
 *      18-Jun-03 - cwb - Solved the basic problem above but
 *         we're still getting rapid increase in PR.  Stepping
 *         through the code and comparing with the spreadsheet
 *         suggests rounding errors in the internal matrices that
 *         perform the decomposition of transpiration to per-group
 *         transpiration values, so I'm making those arrays
 *         double precision with the RealD typedef.  This has
 *         the potential to cause bugs in dynamic memory
 *         routines, so pay attention to which arrays are
 *         defined with double vs single precision and take
 *         appropriate casting measures.
 *	07-16-12 (DLM) - made a ton of changes to try and
 *          get it to compile with the new updated version of soilwat (version 23) 
 *  11-15-19 - Chandler Haukap - The functionality described by cwb on February
 *         28 2002 has been entirely deprecated. I removed the last reference
 *         to SXW_BYMAXSIZE and _Grp_BMass today.
 * 
 * \author CWB (initial coding)
 * \author Chandler Haukap
 * \date 9 May 2002 (initial coding)
 * \ingroup SXW_PRIVATE
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sw_src/generic.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/SW_Defines.h"
#include "sxw.h"
#include "sxw_funcs.h"
#include "sxw_module.h"
#include "sw_src/SW_Control.h"
#include "sw_src/SW_Model.h"
#include "sw_src/SW_VegProd.h"
#include "sw_src/SW_Carbon.h"
#include "sw_src/SW_Site.h"
#include "sw_src/SW_SoilWater.h"
#include "sw_src/SW_Files.h"
#include "sw_src/SW_Weather.h"
#include "sw_src/SW_Markov.h"
#include "sw_src/SW_Output.h"
#include "sw_src/rands.h"
#include "sw_src/pcg/pcg_basic.h"


/*************** Global Variable Declarations ***************/
/***********************************************************/
SXW_t* SXW;

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_VEGPROD SW_VegProd;
extern SW_WEATHER SW_Weather;
extern SW_MARKOV SW_Markov;
//extern SW_SOILWAT SW_Soilwat;

// defined in `SW_Output.c`:
extern Bool prepare_IterationSummary;
extern Bool storeAllIterations;

/*************** Module/Local Variable Declarations ***************/
/***********************************************************/
/* these are initialized and maybe populated here but are used
 * in sxw_resource.c so they aren't declared static.
 */
// Window of transpiration used by _transp_contribution_by_group() in sxw_resource.c
// "Window" refers to the number of years over which transpiration data is averaged.
transp_t* transp_window;
pcg32_random_t resource_rng; //rng for swx_resource.c functions.
SXW_resourceType* SXWResources;

/* These are only used here so they are static.  */
// static char inbuf[FILENAME_MAX];   /* reusable input buffer */
static char _swOutDefName[FILENAME_MAX];
static char *MyFileName;
static char **_sxwfiles[SXW_NFILES];
static char _debugout[256];
static TimeInt _debugyrs[100], _debugyrs_cnt;

/*************** Local Function Declarations ***************/
/***********************************************************/

static void _allocate_memory(void);
static void _read_files( void );
static void _read_roots_max(void);
static void _read_phen(void);
static void _read_prod(void);
static void _read_watin(void);
static void _make_arrays(void);
static void _make_roots_arrays(void);
static void _make_phen_arrays(void);
static void _make_prod_arrays(void);
static void _make_transp_arrays(void);
static void _read_debugfile(void);
void _print_debuginfo(void);
void debugCleanUp(void);
static void _make_swc_array(void);
static void SXW_SW_Setup_Echo(void);
static void SXW_Reinit(char* SOILWAT_file);

void save_sxw_memory( RealD * grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive );
SXW_t* getSXW(void);
SXW_resourceType* getSXWResources(void);
transp_t* getTranspWindow(void);

void copy_sxw_variables(SXW_t* newSXW, SXW_resourceType* newSXWResources, transp_t* newTransp_window);

/****************** Begin Function Code ********************/

/**
 * \brief Read [SOILWAT2](\ref sw_src) input files and initialize variables.
 * 
 * \ingroup SXW
 */
void SXW_Init( Bool init_SW, char *f_roots ) {
  /* read SOILWAT2's input files and initialize some variables */
  /* The shorthand table dimensions are set at the point
   * at which they become defined in _read_files().
   *
   * 2/14/03 - cwb - Added the summation for each group's biomass
   *     as per the notes at the top of the file.
   *
   * 2/21/15 - rjm - new param f_roots is a replacement for the regular
   * 	roots file. This is for the gridded version where soils needs to
   * 	match sxwroots.in
   */
	char roots[MAX_FILENAMESIZE] = { '\0' };

	RandSeed(SuperGlobals.randseed, &resource_rng);

	_allocate_memory(); //Allocate memory for all local pointers

   _sxwfiles[0] = &SXW->f_roots;
   _sxwfiles[1] = &SXW->f_phen;
   _sxwfiles[2] = &SXW->f_prod;
   _sxwfiles[3] = &SXW->f_watin;

  SXW->debugfile = NULL;
  SXW->NGrps = Globals->grpCount;

  _read_files();
  if(f_roots != NULL) {
	  //Copy the directory for the sxwroots files
	  strcpy(roots, DirName(*_sxwfiles[0]));
	  strcat(roots, f_roots);
	  Mem_Free(*_sxwfiles[0]);
	  _sxwfiles[0] = &SXW->f_roots;
	  *_sxwfiles[0] = Str_Dup(roots);
  }
  SXW->NPds = MAX_MONTHS;
  _read_watin();

  if (SXW->debugfile)
	  _read_debugfile();


  if (init_SW)
  {
		SXW_Reinit(SXW->f_watin);
  }

  SXW->NTrLyrs = SW_Site.n_transp_lyrs[0];
  if(SW_Site.n_transp_lyrs[1] > SXW->NTrLyrs)
  	SXW->NTrLyrs = SW_Site.n_transp_lyrs[1];
  if(SW_Site.n_transp_lyrs[3] > SXW->NTrLyrs)
  	SXW->NTrLyrs = SW_Site.n_transp_lyrs[3];
  if(SW_Site.n_transp_lyrs[2] > SXW->NTrLyrs)
	  SXW->NTrLyrs = SW_Site.n_transp_lyrs[2];

  SXW->NSoLyrs = SW_Site.n_layers;

  /* Print general information to stdout. 
     If we are using gridded mode this functionallity will be handled in ST_grid.c */
  if(!UseGrid){
    printf("Number of iterations: %d\n", SuperGlobals.runModelIterations);
    printf("Number of years: %d\n", SuperGlobals.runModelYears);
  }

  _make_arrays();

  _read_roots_max();
  _read_phen();
  _read_prod();

  _sxw_root_phen();

#ifdef TESTING
  _sxw_test();
  exit(0);
#endif
}

/**
 * @brief This function initializes and allocates SOILWAT2 structures,
 *		  and reads SOILWAT2 inputs.
 *
 * \ingroup SXW
 */
static void SXW_Reinit(char* SOILWAT_file) {
	char *temp;

	// setup and construct model (independent of inputs)
	temp = strdup(SOILWAT_file);
	SW_CTL_setup_model(temp);

	// read user inputs
	SW_CTL_read_inputs_from_disk();

	// initialize simulation run (based on user inputs)
	SW_CTL_init_run();

	// initialize output: transfer between STEPPE and SOILWAT2
	SW_OUT_set_SXWrequests();
    free(temp);
}


/**
 * @brief This function resets the model to default conditions
 *
 * It clears SOILWAT2-memory, initializes and allocates SOILWAT2 structures,
 * and reads SOIILWAT2 inputs.
 *
 * However, it does **not** reset memory allocated by
 * `setGlobalSTEPWAT2_OutputVariables` because those variables are carrying
 * over from one STEPWAT2 iteration to the next. They are only de-allocated
 * at the end of an entire STEPWAT2 run (see `ST_main.c/main()`).
 * 
 * \ingroup SXW
 */
void SXW_Reset(char* SOILWAT_file) {
	SW_CTL_clear_model(FALSE); // don't reset output arrays
	SXW_Reinit(SOILWAT_file);
}

/**
 * \brief Resets the SOILWAT2 transpiration tables. This should be called 
 *        every iteration.
 * 
 * \sa Plot_Init() where this function is called.
 * 
 * \ingroup SXW
 */
void SXW_InitPlot (void) {
	_sxw_sw_clear_transp();
	_sxw_update_resource();
}

/**
 * \brief Executes SOILWAT2 which generates soil water resources for plants to 
 *        utilize this year. 
 * 
 * \sa Env_Generate() where this function is called.
 * 
 * \ingroup SXW
 */
void SXW_Run_SOILWAT(void) {
    GrpIndex g;
    Int j;
    SppIndex sp;
    RealF *sizes;

        sizes = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF), "SXW_Run_SOILWAT");

    /* Compute current STEPPE biomass which represents last year's biomass and biomass due to establishment this year (for perennials) and biomass due to establishment this year (for annuals) */
    ForEachGroup(g) {
        sizes[g] = RGroup_GetBiomass(g);

        //printf("First call to sizes: RGroup = %s, sizes[g] = %f\n", RGroup[g]->name, sizes[g]);

        ForEachEstSpp(sp, g, j) {
            
            /* For annual species, increment the biomass that is passed into SOILWAT2 to also include last year's biomass, in addition to biomass due to establishment this year */
            if (Species[sp]->max_age == 1) {

                sizes[g] += Species[sp]->lastyear_relsize * Species[sp]->mature_biomass;
            }
        }
        //printf("Second call to sizes: RGroup = %s, sizes[g] = %f\n", RGroup[g]->name, sizes[g]);

    }

    _sxw_sw_setup(sizes);

    // Initialize `SXW` values for current year's run:
	SXW->aet = 0.; /* used to be in sw_setup() but it needs clearing each run */

    //SXW_SW_Setup_Echo();
    _sxw_sw_run();

    /* Now compute resource availability for each STEPPE functional group */
    _sxw_update_resource();

    /* Set annual precipitation and annual temperature */
    _sxw_set_environs();

    Mem_Free(sizes);
}

void SXW_SW_Setup_Echo(void) {
	char name[256] = {0};
	strcat(name, _debugout);
	FILE *f = OpenFile(strcat(name, ".input.out"), "a");
	int i;
	fprintf(f, "\n================== %d ==============================\n", SW_Model.year);
	fprintf(f,"Fractions Grass:%f Shrub:%f Tree:%f Forb:%f BareGround:%f\n", SW_VegProd.veg[3].cov.fCover, SW_VegProd.veg[1].cov.fCover, SW_VegProd.veg[0].cov.fCover, SW_VegProd.veg[2].cov.fCover, SW_VegProd.bare_cov.fCover);
	fprintf(f,"Monthly Production Values\n");
	fprintf(f,"Grass\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.veg[3].litter[i],
				SW_VegProd.veg[3].biomass[i], SW_VegProd.veg[3].pct_live[i],
				SW_VegProd.veg[3].lai_conv[i]);
	}

	fprintf(f,"Shrub\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.veg[1].litter[i],
				SW_VegProd.veg[1].biomass[i], SW_VegProd.veg[1].pct_live[i],
				SW_VegProd.veg[1].lai_conv[i]);
	}

	fprintf(f,"Tree\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.veg[0].litter[i],
				SW_VegProd.veg[0].biomass[i], SW_VegProd.veg[0].pct_live[i],
				SW_VegProd.veg[0].lai_conv[i]);
	}

	fprintf(f,"Forb\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.veg[2].litter[i],
				SW_VegProd.veg[2].biomass[i], SW_VegProd.veg[2].pct_live[i],
				SW_VegProd.veg[2].lai_conv[i]);
	}

	SW_SITE *s = &SW_Site;
	fprintf(f,"Soils Transp_coeff\n");
	fprintf(f,"Forb\tTree\tShrub\tGrass\n");
	ForEachSoilLayer(i)
	{// %u %u %u %u s->lyr[i]->my_transp_rgn_forb, s->lyr[i]->my_transp_rgn_tree, s->lyr[i]->my_transp_rgn_shrub, s->lyr[i]->my_transp_rgn_grass
		fprintf(f,"%6.2f %6.2f %6.2f %6.2f\n", s->lyr[i]->transp_coeff[2], s->lyr[i]->transp_coeff[0], s->lyr[i]->transp_coeff[1], s->lyr[i]->transp_coeff[3]);
	}

  // adding values to sxw structure for use in ST_stats.c
  /*SXW->grass_cover = SW_VegProd.grass.conv_stcr;
  SXW->shrub_cover = SW_VegProd.shrub.conv_stcr;
  SXW->tree_cover = SW_VegProd.tree.conv_stcr;
  SXW->forbs_cover = SW_VegProd.forb.conv_stcr;*/


	fprintf(f, "\n");
	CloseFile(&f);
}

/**
 * \brief Obtain transpiration (resource availability) for each resource group.
 * 
 * \param rg is the [index](\ref GrpIndex) in \ref RGroup of the resource group.
 * 
 * \return The transpiration availible for \ref RGroup[rg].
 */
RealF SXW_GetTranspiration( GrpIndex rg) {
	return SXWResources->_resource_cur[rg];
}

void SXW_PrintDebug(Bool cleanup) {
/*======================================================*/
	TimeInt i;
	static Bool beenhere = FALSE;

	if(cleanup) {
		debugCleanUp();
	} else {
		for (i = 0; i < _debugyrs_cnt; i++) {
			if (SW_Model.year == _debugyrs[i]) {
				SXW_SW_Setup_Echo();
				_print_debuginfo();
				break;
			}
		}
		if (!beenhere) {
			beenhere = TRUE;
			insertInfo();
			insertSXWPhen();
			insertSXWProd();
			insertRootsXphen(SXWResources->_rootsXphen);
		}
		insertInputVars();
		insertInputProd();
		insertInputSoils();
		insertOutputVars(SXWResources->_resource_cur, transp_window->added_transp);
		insertRgroupInfo(SXWResources->_resource_cur);
		insertOutputProd(&SW_VegProd);
		insertRootsSum(SXWResources->_roots_active_sum);
		insertRootsRelative(SXWResources->_roots_active_rel);
		insertTranspiration();
		insertSWCBulk();
	}
}

/** \brief Allocate memory for sxw local variables. 
 * 
 * Local variables allocated are SXW, transp_window, and SXWResources.
 * However, memory allocation in SXW is not modularized well. Check
 * \ref _make_arrays(), \ref _make_roots_arrays(),
 * \ref _make_phen_arrays(), \ref _make_prod_arrays(),
 * \ref _make_transp_arrays(), and \ref _make_swc_array() before assuming that
 * something isn't allocated.
 * 
 * \ingroup SXW_private
 */
static void _allocate_memory(void){
	transp_window = (transp_t*) Mem_Calloc(1, sizeof(transp_t), "_allocate_memory: transp_window");

	/* Set the size of our transpiration window. */
	if(Globals->transp_window > MAX_WINDOW){
		LogError(logfp, LOGWARN, "Requested transp_window (%d) > %d. Setting transp_window to %d.", 
		         Globals->transp_window, MAX_WINDOW, MAX_WINDOW);
		transp_window->size = MAX_WINDOW;
	} else {
		transp_window->size = Globals->transp_window;
	}

	transp_window->ratios = (RealF*) Mem_Calloc(transp_window->size, sizeof(RealF), "_allocate_memory: transp_window->ratios");
	transp_window->transp = (RealF*) Mem_Calloc(transp_window->size, sizeof(RealF), "_allocate_memory: transp_window->transp");
	transp_window->SoS_array = (RealF*) Mem_Calloc(transp_window->size, sizeof(RealF), "_allocate_memory: transp_window->SoS_array");

	SXW = (SXW_t*) Mem_Calloc(1, sizeof(SXW_t), "_allocate_memory: SXW");
	SXWResources = (SXW_resourceType*) Mem_Calloc(1, sizeof(SXW_resourceType), "_allocate_memory: SXWResources");

	SXWResources->_resource_cur = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF), "_allocate_memory: _resource_cur");
}

/* Returns a pointer to the local SXW variable. */
SXW_t* getSXW(void){
	return SXW;
}

/* Returns a pointer to the local SXWResources variable. */
SXW_resourceType* getSXWResources(void){
	return SXWResources;
}

/* Returns a pointer to the local transp_window variable. */
transp_t* getTranspWindow(void){
	return transp_window;
}

/* Shallow copy variables into the local sxw variables. */
void copy_sxw_variables(SXW_t* newSXW, SXW_resourceType* newSXWResources, transp_t* newTransp_window){
	SXW = newSXW;
	SXWResources = newSXWResources;
	transp_window = newTransp_window;
}

static void  _read_files( void ) {
/*======================================================*/
  /* read list of input files.   */
  FILE *fin;
  int i, nfiles = SXW_NFILES;

  SXW->f_files = Parm_name(F_SXW);  /* aliased */
  MyFileName = SXW->f_files;
  fin = OpenFile(MyFileName,"r");

  for(i=0; i < nfiles; i++) {
    if (!GetALine(fin, inbuf)) break;
    *_sxwfiles[i] = Str_Dup(Str_TrimLeftQ(inbuf));
  }

  if (i < nfiles) {
    LogError(logfp, LOGFATAL, "STEPWAT: %s: Insufficient files found",
            MyFileName);
  }

  CloseFile(&fin);


}

static void  _read_roots_max(void) {
/*======================================================*/
	GrpIndex g;
	int cnt = 0, lyr;
	char *p;
	char *name;
	FILE *fp;
        
        name = (char *)Mem_Calloc(SuperGlobals.max_groupnamelen + 1, sizeof(char), "_read_roots_max");

	MyFileName = SXW->f_roots;
	fp = OpenFile(MyFileName, "r");

	while (GetALine(fp, inbuf)) {
		p = strtok(inbuf, " \t"); /* g'teed to work via GetALine() */

		if ((g = RGroup_Name2Index(p)) < 0) {
			LogError(logfp, LOGFATAL, "%s: Invalid group name (%s) found.",
					MyFileName, p);
		}
		strcpy(name, p);
		cnt++;
		lyr = 0;
		while ((p = strtok(NULL, " \t"))) {
			SXWResources->_roots_max[Ilg(lyr, g)] = atof(p);
			lyr++;
		}
		if (lyr != SXW->NTrLyrs) {
			LogError(logfp, LOGFATAL,
					"%s: Group : %s : Missing layer values. Match up with soils.in file. Include zeros if necessary. Layers needed %u. Layers defined %u",
					MyFileName, name, SXW->NTrLyrs, lyr);
		}
	}

	if (cnt < Globals->grpCount) {
		LogError(logfp, LOGFATAL, "%s: Not enough valid groups found.",
				MyFileName);
	}

	CloseFile(&fp);
        
        Mem_Free(name);
}

static void _read_phen(void) {
/*======================================================*/

  GrpIndex g;
  IntUS cnt=0;
  TimeInt m;
  char *p;
  FILE *fp;

  MyFileName = SXW->f_phen;
  fp = OpenFile(MyFileName,"r");


  while( GetALine(fp, inbuf) ) {
    p = strtok(inbuf," \t"); /* g'teed to work via GetALine() */

    if ( (g=RGroup_Name2Index(p)) <0 ) {
      LogError(logfp, LOGFATAL,
               "%s: Invalid group name (%s) found.", MyFileName, p);
    }
    cnt++;
    m = Jan;
    while ((p=strtok(NULL," \t")) ) {
      if (m > Dec) {
        LogError(logfp, LOGFATAL,
                 "%s: More than 12 months of data found.", MyFileName);
      }
      SXWResources->_phen[Igp(g,m)] = atof(p);
      m++;
    }

  }

  if (cnt < Globals->grpCount) {
    LogError(logfp, LOGFATAL,
             "%s: Not enough valid groups found.", MyFileName);
  }

  CloseFile(&fp);

}

/** \brief Read the values from the SXW prod file.
 * 
 * These values are LITTER, BIOMASS, and PCTLIVE. All three values are input
 * for each group for each month in separate two dimensional tables. If a table
 * if poorly formated, i.e. incorrect \ref RGroup names, incorrect number of
 * months, or incorrect number of \ref RGroups, this function will throw a 
 * fatal error.
 * 
 * \sideeffect
 *         Populates multiple one dimensional and two dimensional arrays.
 * 
 * \ingroup SXW_private
 */
static void _read_prod(void) {
	FILE *fp;
	TimeInt month = Jan;

	GrpIndex g;
	IntUS count = 0;
	char *p;
	char * pch;
	MyFileName = SXW->f_prod;
	fp = OpenFile(MyFileName, "r");

    /* Read LITTER values for each group for each month */
	while (GetALine(fp, inbuf)) {
		pch = strstr(inbuf, "[end]");
		if(pch != NULL) break;

		p = strtok(inbuf, " \t"); /* guaranteed to work via GetALine() */
		if ((g = RGroup_Name2Index(p)) < 0) {
			LogError(logfp, LOGFATAL, "%s: Invalid group name for LITTER (%s).",
					MyFileName, p);
		}

		month = Jan;
		while ((p = strtok(NULL, " \t"))) {
			if (month > Dec) {
				LogError(logfp, LOGFATAL,
						"%s: More than 12 months of data found for LITTER.", MyFileName);
			}
			SXWResources->_prod_litter[g][month] = atof(p);
			month++;
		}

		count++;
		if(count == SXW->NGrps)
			break;
	}

	if (count < Globals->grpCount) {
		LogError(logfp, LOGFATAL, "%s: Not enough valid groups found for LITTER values.",
				 MyFileName);
	}

    /* Read BIOMASS values for each group for each month */
	GetALine(fp,inbuf); /* toss [end] keyword */
	month = Jan;
    count = 0;

	while (GetALine(fp, inbuf)) {
		pch = strstr(inbuf, "[end]");
		if(pch != NULL)
			break;
		p = strtok(inbuf, " \t"); /* guaranteed to work via GetALine() */

		if ((g = RGroup_Name2Index(p)) < 0) {
			LogError(logfp, LOGFATAL, "%s: Invalid group name for biomass (%s) found.",
					MyFileName, p);
		}
		month = Jan;
		while ((p = strtok(NULL, " \t"))) {
			if (month > Dec) {
				LogError(logfp, LOGFATAL,
						"%s: More than 12 months of data found.", MyFileName);
			}
			SXWResources->_prod_bmass[Igp(g, month)] = atof(p);
			month++;
		}
		count++;
		if(count == SXW->NGrps)
			break;
	}

	if (count < Globals->grpCount) {
		LogError(logfp, LOGFATAL, "%s: Not enough valid groups found.",
				MyFileName);
	}

    /* Read PCTLIVE values for each group for each month */
	GetALine(fp, inbuf); /* toss [end] keyword */
	month = Jan;
	count=0;

	while (GetALine(fp, inbuf)) {
		pch = strstr(inbuf, "[end]");
		if (pch != NULL)
			break;
		p = strtok(inbuf, " \t"); /* guaranteed to work via GetALine() */

		if ((g = RGroup_Name2Index(p)) < 0) {
			LogError(logfp, LOGFATAL,
					"%s: Invalid group name for pctlive (%s) found.",
					MyFileName, p);
		}
		month = Jan;
		while ((p = strtok(NULL, " \t"))) {
			if (month > Dec) {
				LogError(logfp, LOGFATAL,
						"%s: More than 12 months of data found.", MyFileName);
			}
			SXWResources->_prod_pctlive[Igp(g, month)] = atof(p);
			month++;
		}
		count++;
		if (count == SXW->NGrps)
			break;
	}

	if (count < Globals->grpCount) {
		LogError(logfp, LOGFATAL, "%s: Not enough valid groups found.",
				MyFileName);
	}

	GetALine(fp, inbuf); /* toss [end] keyword */
	CloseFile(&fp);
}

static void _read_watin(void) {
/*======================================================*/
/* get the name of the soilwat output definition file.  */
/* assume that there is no path prepended to the file
 * specification in the SOILWAT input files.  It might be
 * nice to allow relative paths, but there's really no need
 * for it as far as the STEPWAT program is concerned because
 * all the SOILWAT input files should be in one directory
 * and that is defined in sxw.in.  Thus, we'll treat the
 * outsetup.in filename as though it has no dirname and
 * append it to the soilwat input dirname.
 */
   FILE *f;
   int lineno = 0;
   Bool found = FALSE;

   MyFileName = SXW->f_watin;
   f = OpenFile(MyFileName, "r");


   while( GetALine(f, inbuf) ) {
     if (++lineno == (eOutput + 2)) {
       strcpy(_swOutDefName, DirName(SXW->f_watin));
       strcat(_swOutDefName, inbuf);
       found = TRUE;
       break;
     }
   }
   CloseFile(&f);

   if (!found) {
     LogError(logfp, LOGFATAL,
              "%s: Too few files (%d)", MyFileName, lineno);
   }

}

static void _make_arrays(void) {
/*======================================================*/
/* central point to make all dynamically allocated arrays
 * now that the dimensions are known.
 */
	_make_roots_arrays();
	_make_phen_arrays();
	_make_prod_arrays();
	_make_transp_arrays();
	_make_swc_array();
}

static void _make_roots_arrays(void) {
/*======================================================*/
  int size;
  char *fstr = "_make_roots_array()";

  size  = SXW->NGrps * SXW->NTrLyrs;
  SXWResources->_roots_max     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

  size = SXW->NGrps * SXW->NPds * SXW->NTrLyrs;
  SXWResources->_rootsXphen       = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  SXWResources->_roots_active     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  SXWResources->_roots_active_rel = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

  //4 - Grass,Frob,Tree,Shrub
  size = NVEGTYPES * SXW->NPds * SXW->NTrLyrs;
  SXWResources->_roots_active_sum = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
}

static void _make_phen_arrays(void) {
/*======================================================*/
  int size;
  char *fstr = "_make_phen_arrays()";

  size = SXW->NGrps * MAX_MONTHS;
  SXWResources->_phen = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

}

/** \brief Allocate the "prod" arrays. 
 * 
 * _prod_bmass, _prod_pctlive, and _prod_litter from the \ref SXWResources
 * struct are all allocated here.
 * 
 * \sideeffect two arrays and one 2D array will be allocated.
 * 
 * \ingroup
 */
static void _make_prod_arrays(void) {
	int size, i;
	char *fstr = "_make_phen_arrays()";

	size = SXW->NGrps * MAX_MONTHS;
	SXWResources->_prod_bmass = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	SXWResources->_prod_pctlive = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

    /* Allocate the _prod_litter 2D array, where rows are rgroups and columns are months. */
    SXWResources->_prod_litter = (RealD**) Mem_Calloc(SXW->NGrps, sizeof(RealD*), fstr);
    for(i = 0; i < SXW->NGrps; ++i){
        SXWResources->_prod_litter[i] = (RealD*) Mem_Calloc(MAX_MONTHS, sizeof(RealD), fstr);
    }
}

static void _make_transp_arrays(void) {
/*======================================================*/
/* SXW->transp holds year-to-year values and is populated in SOILWAT.
 * both are indexed by the macro Ilp().
 */
	char *fstr = "_make_transp_array()";
	int size, k;

  size = SXW->NPds * SXW->NSoLyrs; // monthly values for each soil layer
	SXW->transpTotal = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

	ForEachVegType(k) {
		SXW->transpVeg[k] = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	}
}

static void _make_swc_array(void) {
/*======================================================*/
/* SXW->swc holds year-to-year values and is populated in SOILWAT.
 * it is indexed by the macro Ilp(). only used if soilwat option
 * specified with debugfile
 */
	char *fstr = "_make_swc_array()";
	int size = SXW->NPds * SXW->NSoLyrs; // monthly values for each soil layer

	SXW->swc = (RealF *) Mem_Calloc(size, sizeof(RealF *), fstr);
}

static void _read_debugfile(void) {
/*======================================================*/
  /* provides a way to specify optional debug information
   * to be printed at various points in the model run.
   *
   * 7-Jan-03 - format of file is as follows:
   *   debugfilename #name of file to write debuginfo to (len<256 chars)
   *   yyyy yyyy yyyy ... # 4-digit years to write output (see years.in)
   *   yyyy yyyy yyyy ... # possible multi-line continuation of above
   *
   *  Notes:
   *  - LIMIT of 100 4-digit years
   *
   *  - You'll probably want to limit the number of model iterations to 1.
   *
   *  - Start year is defined in SOILWAT's model parameter file
   *    (probably years.in).
   *
   *  -If you want to add more sections, use something like [end] after
   *   each and test to break out of the current section's loop. Eg,
   *     while (GetALine(f,inbuf) && !strcmp(inbuf, "[end]") ) {...}
   *   with the above example, don't forget to test for end of file (EOF)
   *   before continuing.
   */

	FILE *f;
	char *date, str[102];
	int cnt = 0;
	TimeInt i;
	char name[256] = {0};

	f = OpenFile(SXW->debugfile, "r");

	/* get name of output file */
	if (!GetALine(f, inbuf)) {
		CloseFile(&f);
		return;
	}
	strcpy(_debugout, inbuf);

	/* get output years */
	while (GetALine(f, inbuf)) {
		_debugyrs[cnt++] = atoi(strtok(inbuf, " \t")); /* g'teed via getaline() */
		while (NULL != (date = strtok(NULL, " \t"))) {
			_debugyrs[cnt++] = atoi(date);
		}
	}
	_debugyrs_cnt = cnt;

	sprintf(errstr, "Debugging Transpiration turned on.\n%s will contain"
			" %d years of output:\n", _debugout, _debugyrs_cnt);

	for (i = 0; i < _debugyrs_cnt; i++) {
		sprintf(str, "%d\n", _debugyrs[i]);
		strcat(errstr, str);
	}
	strcat(errstr, "Note that data will always be appended,\n");
	strcat(errstr, "so clear file contents before re-use.\n");
	LogError(logfp, LOGNOTE, errstr);

	CloseFile(&f);

	/* now empty the file prior to the run */
	strcat(name, _debugout);
	f = OpenFile(strcat(name, ".output.out"), "w");
	CloseFile(&f);

	name[0] = 0;
	strcat(name, _debugout);
	f = OpenFile(strcat(name, ".input.out"), "w");
	CloseFile(&f);

	connect(_debugout);
	createTables();
}

void debugCleanUp() {
  printf("in debugCleanUp\n");
	disconnect();
}

void _print_debuginfo(void) {
/*======================================================*/
	SW_VEGPROD *v = &SW_VegProd;
	TimeInt p;
	LyrIndex t;
	int l;
	FILE *f;
	GrpIndex r;
	RealF sum = 0.;
	RealF sum1 = 0.;
	RealF sum2 = 0.;
	RealF sum3 = 0.;
	static Bool beenhere = FALSE;
	char vegProdNames[NVEGTYPES][7];
	strcpy(vegProdNames[SW_TREES], "TREE");
	strcpy(vegProdNames[SW_SHRUB], "SHRUB");
	strcpy(vegProdNames[SW_GRASS], "GRASS");
	strcpy(vegProdNames[SW_FORBS], "FORB");
	char name[256] = {0};
	strcat(name, _debugout);
	f = OpenFile(strcat(name, ".output.out"), "a");

	if (!beenhere) {
		beenhere = TRUE;
		fprintf(f, "\n------ Roots X Phen Array -------\n");
		ForEachGroup(r)
		{
			fprintf(f, "         --- %s ---\n", RGroup[r]->name);

			fprintf(f, "Layer:");
			ForEachTrPeriod(p)
				fprintf(f, "\t%d", p + 1);
			fprintf(f, "\n");
			int nLyrs = getNTranspLayers(RGroup[r]->veg_prod_type);
			for (t = 0; t < nLyrs; t++) {
				//ForEachTranspLayer(t) {
				fprintf(f, "%d", t + 1);
				ForEachTrPeriod(p)
					fprintf(f, "\t%.4f", SXWResources->_rootsXphen[Iglp(r, t, p)]);
				fprintf(f, "\n");
			}
		}
	}


  /* sum actual total transpiration */
	ForEachTrPeriod(p)
	{
		for (t = 0; t < SXW->NSoLyrs; t++)
			sum += SXW->transpTotal[Ilp(t, p)];
	}

	fprintf(f, "\n================== %d =============================\n", SW_Model.year);
	fprintf(f, "MAP = %d(mm)\tMAT = %5.2f(C)\tAET = %5.4f(cm)\tT = %5.4f(cm)\tTADDED = %5.4f(cm)\tAT = %5.4f(cm)\n\n", Env->ppt, Env->temp, SXW->aet, sum, transp_window->added_transp, sum + transp_window->added_transp);

	fprintf(f, "Group     \tRelsize\tPR\tResource_cur\tResource_cur\n");
	fprintf(f, "-----     \t-------\t-----\t-no scaling-\t-with scaling-\n");
	ForEachGroup(r) {
		sum1 += getRGroupRelsize(r);
		sum2 += RGroup[r]->pr;
		sum3 += SXWResources->_resource_cur[r];
		fprintf(f, "%s\t%.4f\t%.4f\t%.4f\t\t%.4f\n", RGroup[r]->name, getRGroupRelsize(r),
                RGroup[r]->pr, SXWResources->_resource_cur[r]/RGroup[r]->_bvt, 
                SXWResources->_resource_cur[r]);
	}

	fprintf(f, "\n------ Production Values Daily Summed Across Types Monthly Averaged -------\n");
	fprintf(f, "Month\tBMass\tPctLive\tLAIlive\tLAItotal\tTotAGB\n");
	fprintf(f, "-----\t-----\t-------\t-------\t------\t------\n");

	int doy = 1;
	ForEachMonth(p)
	{
		int days = 31, i;
		double pct_live = 0, lai_live = 0, bLAI_total = 0, total_agb = 0, biomass = 0;

		if (p == Apr || p == Jun || p == Sep || p == Nov) //all these months have 30 days
			days = 30;
		else if (p == Feb) { //February has either 28 or 29 days
			days = 28;
			if (isleapyear(SW_Model.year))
				days = 29;
		} // all the other months have 31 days

		for (i = doy; i < (doy + days); i++) { //accumulating the monthly values...
			lai_live += (v->veg[0].lai_live_daily[i])
					+ (v->veg[1].lai_live_daily[i])
					+ (v->veg[3].lai_live_daily[i])
					+ (v->veg[2].lai_live_daily[i]);
			bLAI_total += (v->veg[0].bLAI_total_daily[i]) + (v->veg[1].bLAI_total_daily[i])
					+ (v->veg[3].bLAI_total_daily[i]) + (v->veg[2].bLAI_total_daily[i]);
			total_agb += (v->veg[0].total_agb_daily[i])
					+ (v->veg[1].total_agb_daily[i])
					+ (v->veg[3].total_agb_daily[i])
					+ (v->veg[2].total_agb_daily[i]);
			pct_live += (v->veg[0].pct_live_daily[i])
					+ (v->veg[1].pct_live_daily[i])
					+ (v->veg[3].pct_live_daily[i])
					+ (v->veg[2].pct_live_daily[i]);
			biomass += (v->veg[0].biomass_daily[i])
					+ (v->veg[1].biomass_daily[i])
					+ (v->veg[3].biomass_daily[i])
					+ (v->veg[2].biomass_daily[i]);
		}
		doy += days; //updating the doy
		//biomass = (v->tree.biomass[p]) + (v->shrub.biomass[p])
		//		+ (v->grass.biomass[p]) + (v->forb.biomass[p]);		//getting the monthly biomass
		//pct_live = (v->tree.pct_live[p]) + (v->shrub.pct_live[p])
		//		+ (v->grass.pct_live[p]) + (v->forb.pct_live[p]);	//getting the monthly pct_live

		pct_live /= days;
		biomass /= days;
		lai_live /= days; //getting the monthly averages...
		bLAI_total /= days;
		total_agb /= days;

		fprintf(f, "%4d\t%.2f\t%.3f\t%.3f\t%.2f\t%.2f\n", p + 1, biomass,
				pct_live, lai_live, bLAI_total, total_agb);
	}

	ForEachVegType(t) {
		fprintf(f, "\n------ Active Roots (sum) %s -------\n", vegProdNames[t]);
		fprintf(f, "Layer:");
		ForEachTrPeriod(p)
			fprintf(f, "\t%d", p + 1);
		fprintf(f, "\n");
		for (l = 0; l < SXW->NTrLyrs; l++) {
			fprintf(f, "%d", t + 1);
			ForEachTrPeriod(p)
				fprintf(f, "\t%.4f", SXWResources->_roots_active_sum[Itlp(t, l, p)]);
			fprintf(f, "\n");
		}
	}

	fprintf(f, "\n------ Active Roots (relative) -------\n");
	ForEachGroup(r)
	{
		fprintf(f, "         --- %s ---\n", RGroup[r]->name);

		fprintf(f, "Layer:");
		ForEachTrPeriod(p)
			fprintf(f, "\t%d", p + 1);
		fprintf(f, "\n");
		int nLyrs = getNTranspLayers(RGroup[r]->veg_prod_type);
		for (t = 0; t < nLyrs; t++) {
			fprintf(f, "%d", t + 1);
			ForEachTrPeriod(p)
				fprintf(f, "\t%.4f", SXWResources->_roots_active_rel[Iglp(r, t, p)]);
			fprintf(f, "\n");
		}
	}

	fprintf(f, "\n------ Transpiration Total Values -------\nPeriod:");
	for (t = 0; t < SXW->NSoLyrs; t++)
		fprintf(f, "\t\tL%d", t + 1);
	fprintf(f, "\n");

	ForEachTrPeriod(p)
	{
		fprintf(f, "%d : ", p + 1);
		for (t = 0; t < SXW->NSoLyrs; t++)
			fprintf(f, "\t%.4f", SXW->transpTotal[Ilp(t, p)]);
		fprintf(f, "\n");
	}

	fprintf(f, "Current Soil Water Content:\n");
	ForEachTrPeriod(p)
	{
		fprintf(f, "%d : ", p + 1);
		ForEachSoilLayer(t)
			fprintf(f, "\t%5.4f", SXW->swc[Ilp(t, p)]);
		fprintf(f, "\n");
	}

	fprintf(f, "\n");
	CloseFile(&f);
}


#ifdef DEBUG_MEM
#include "sw_src/myMemory.h"
/*======================================================*/
void SXW_SetMemoryRefs( void) {
/* when debugging memory problems, use the bookkeeping
   code in myMemory.c
 This routine sets the known memory refs so they can be
 checked for leaks, etc.  Includes malloc-ed memory from
 SXW as well as SOILWAT.  All refs will have been cleared
 by a call to ClearMemoryRefs() before this, and will be
 checked via CheckMemoryRefs() after this, most likely in
 the main() function.
*/
  TimeInt p;
  int i, last = SXW_NFILES-1;  /* recall we skipped the first file */

   for (i=0; i < last; i++) {
     NoteMemoryRef(*_sxwfiles[i]);
   }
   NoteMemoryRef(_roots_max);
   NoteMemoryRef(_roots_rel);
   NoteMemoryRef(_roots_current);
   NoteMemoryRef(_rootsXphen);
   NoteMemoryRef(_roots_active);
   NoteMemoryRef(_roots_phen_totals);
   NoteMemoryRef(_phen);
   NoteMemoryRef(_phen_grp_rel);
   NoteMemoryRef(_transp_grp_totals);
   NoteMemoryRef(SXW->transp);


   SW_CTL_SetMemoryRefs();
}

#endif

/** Convert STEPWAT2 indices of SOILWAT2's vegetation type into a
    SOILWAT2 index

    @param veg_prod_type 1 for tree, 2 for shrub, 3 for grass, 4 for forb
      (see comments for "veg_prod_type" in `rgroup.in`).

    @return One of the SOILWAT2 defined values `SW_TREES`, `SW_SHRUB`,
      `SW_FORBS`, `SW_GRASS` or -1 if no match.
      See `SW_Defines.h` for definitions.
*/
int get_SW2_veg_index(int veg_prod_type) {
  if (veg_prod_type == 1)
    return SW_TREES;
  else if (veg_prod_type == 2)
    return SW_SHRUB;
  else if (veg_prod_type == 3)
    return SW_GRASS;
  else if (veg_prod_type == 4)
    return SW_FORBS;

  return -1;
}



 /***********************************************************/
 //returns the number of transpiration layers correctly for each veg_prod_type
int getNTranspLayers(int veg_prod_type) {
  return SW_Site.n_transp_lyrs[veg_prod_type];
}

/** 
 * \brief Free the memory allocated to the SXW, transp_window, and 
 *         SXWResources structs.
 * 
 * This function will free the memory of all pointers in the structs
 * then free the structs themselves. This function should always be preceded
 * by a call to \ref SXW_Init().
 * 
 * \sideeffect 
 *         All three variables mentioned will be completely deallocated.
 * 
 * \author Chandler Haukap
 * 
 * \ingroup SXW
 */
void free_all_sxw_memory( void ) {
	int k;

	/* Free transp_window */
	Mem_Free(transp_window->ratios);
	Mem_Free(transp_window->transp);
	Mem_Free(transp_window->SoS_array);
	Mem_Free(transp_window);

    /* Free SXWResources */
	Mem_Free(SXWResources->_phen);
	Mem_Free(SXWResources->_prod_bmass);
	Mem_Free(SXWResources->_prod_pctlive);
	Mem_Free(SXWResources->_resource_cur);
	Mem_Free(SXWResources->_roots_active);
	Mem_Free(SXWResources->_roots_active_rel);
	Mem_Free(SXWResources->_roots_active_sum);
	Mem_Free(SXWResources->_roots_max);
	Mem_Free(SXWResources->_rootsXphen);
    for(k = 0; k < SXW->NGrps; ++k){
        Mem_Free(SXWResources->_prod_litter[k]);
    }
    Mem_Free(SXWResources->_prod_litter);
	Mem_Free(SXWResources);

	/* Free SXW */
	Mem_Free(SXW->f_roots);
	Mem_Free(SXW->f_phen);
	Mem_Free(SXW->f_prod);
	Mem_Free(SXW->f_watin);
	Mem_Free(SXW->transpTotal);
	ForEachVegType(k) {
		Mem_Free(SXW->transpVeg[k]);
	}
	Mem_Free(SXW->swc);
	Mem_Free(SXW);
}
