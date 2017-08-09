/********************************************************/
/********************************************************/
/*  Source file: sxw.c
 *
/*  Type: module
 *
/*  Purpose: Interface module for the STEPPE to SOILWAT
 *           data flow.  Oversees transformation of
 *           data from STEPPE to SOILWAT.
 *
/*  Calls:  sxw2wat.c
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *
/*  History:
/*     (9-May-2002) -- INITIAL CODING - cwb
 *
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
 *
 *		07-16-12 (DLM) - made a ton of changes to try and get it to compile with the new updated version of soilwat (version 23)
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "SW_Defines.h"
#include "sxw.h"
#include "sxw_funcs.h"
#include "sxw_module.h"
#include "sw_src/SW_Control.h"
#include "sw_src/SW_Model.h"
#include "sw_src/SW_VegProd.h"
#include "sw_src/SW_Site.h"
#include "sw_src/SW_SoilWater.h"
#include "sw_src/SW_Files.h"
#include "sw_src/SW_Weather.h"
#include "sw_src/SW_Markov.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
SXW_t SXW;

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_VEGPROD SW_VegProd;
extern SW_WEATHER SW_Weather;
extern SW_MARKOV SW_Markov;
//extern SW_SOILWAT SW_Soilwat;


/*************** Module/Local Variable Declarations ***************/
/***********************************************************/
/* these are initialized and maybe populated here but are used
 * in sxw_resource.c so they aren't declared static.
 */
/* ----- 3d arrays ------- */
RealD * _rootsXphen, /* relative roots X phen in each lyr,grp,pd */
      * _roots_active, /* "active" in terms of size and phenology */
      * _roots_active_rel;


/* ----- 2D arrays ------- */
/* malloc'ed here for consistency but only used */
/* in sxw_resource.c and sxw_soilwat.c */

     /* rgroup by layer, ie, group-level values */
RealD * _roots_max,     /* read from root distr. file */
      * _roots_active_sum, /* used in sxw_resource */

     /* rgroup by period */
      * _phen;          /* phenology read from file */

/* simple vectors hold the resource information for each group */
/* curr/equ gives the available/required ratio */
RealF _resource_cur[MAX_RGROUPS],  /* current resource utilization */
      _resource_pr[MAX_RGROUPS];   /* resource convertable to PR */

#ifdef SXW_BYMAXSIZE
/* addition to meet changes specified at the top of the file */
RealF _Grp_BMass[MAX_RGROUPS];
#endif

/* and one vector for the production constants */
RealD _prod_litter[MAX_MONTHS];
RealD * _prod_bmass;
RealD * _prod_pctlive;

RealF _bvt;  /* ratio of biomass/m2 / transp/m2 */

/* These are only used here so they are static.  */
// static char inbuf[FILENAME_MAX];   /* reusable input buffer */
static char _swOutDefName[FILENAME_MAX];
static char *MyFileName;
static char **_sxwfiles[SXW_NFILES];
static char _debugout[256];
static TimeInt _debugyrs[100], _debugyrs_cnt;


/*************** Local Function Declarations ***************/
/***********************************************************/

static void _read_files( void );
static void _read_roots_max(void);
static void _read_phen(void);
static void _read_prod(void);
static void _read_bvt(void);
static void _read_watin(void);
static void _make_arrays(void);
static void _make_roots_arrays(void);
static void _make_phen_arrays(void);
static void _make_prod_arrays(void);
static void _make_transp_arrays(void);
static void _write_sw_outin(void);
//static void _recover_names(void);
static void _read_debugfile(void);
void _print_debuginfo(void);
void debugCleanUp(void);
static void _make_swc_array(void);
static void SXW_SW_Setup_Echo(void);
//static void SXW_SW_Output_Echo(void);

//these last four functions are to be used in ST_grid.c
void load_sxw_memory( RealD * grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive );
void save_sxw_memory( RealD * grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive );
void free_sxw_memory( void );
void free_all_sxw_memory( void );

/****************** Begin Function Code ********************/
/***********************************************************/
void SXW_Init( Bool init_SW, char *f_roots ) {
  /* read SOILWAT's input files and initialize some variables */
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
	char * temp;
	char roots[512] = { '\0' };

#ifdef SXW_BYMAXSIZE
   GrpIndex rg; SppIndex sp;
   /* Sum each group's maximum biomass */
   ForEachGroup(rg) _Grp_BMass[rg] = 0.0;
   ForEachSpecies(sp)
     _Grp_BMass[Species[sp]->res_grp] += Species[sp]->mature_biomass;
   /* end code 2/14/03 */
#endif

   _sxwfiles[0] = &SXW.f_roots;
   _sxwfiles[1] = &SXW.f_phen;
   _sxwfiles[2] = &SXW.f_bvt;
   _sxwfiles[3] = &SXW.f_prod;
   _sxwfiles[4] = &SXW.f_watin;

  SXW.NGrps = Globals.grpCount;

  _read_files();
  if(f_roots != NULL) {
	  //Copy the directory for the sxwroots files
	  strcpy(roots, DirName(*_sxwfiles[0]));
	  strcat(roots, f_roots);
	  Mem_Free(*_sxwfiles[0]);
	  _sxwfiles[0] = &SXW.f_roots;
	  *_sxwfiles[0] = Str_Dup(roots);
  }
  SXW.NPds = MAX_MONTHS;
  _read_watin();

  if (SXW.debugfile)
	  _read_debugfile();
  //_write_sw_outin();


  if(init_SW) {
	  temp = strdup(SXW.f_watin);
	  SW_CTL_init_model(temp);
	  free(temp);
  }

  SXW.NTrLyrs = SW_Site.n_transp_lyrs_tree;
  if(SW_Site.n_transp_lyrs_shrub > SXW.NTrLyrs)
  	SXW.NTrLyrs = SW_Site.n_transp_lyrs_shrub;
  if(SW_Site.n_transp_lyrs_grass > SXW.NTrLyrs)
  	SXW.NTrLyrs = SW_Site.n_transp_lyrs_grass;
  if(SW_Site.n_transp_lyrs_forb > SXW.NTrLyrs)
	  SXW.NTrLyrs = SW_Site.n_transp_lyrs_forb;

  SXW.NSoLyrs = SW_Site.n_layers;
  printf("Number of layers: %d\n", SW_Site.n_layers);

  // TODO: need to see if need to have output for every run of STEPPE or only when -o parameter is called
  /*stat_Output_Daily_CSV_Summary();
  stat_Output_Weekly_CSV_Summary();
  stat_Output_Monthly_CSV_Summary();
  stat_Output_Yearly_CSV_Summary();
  */

  _make_arrays();

  _read_roots_max();
  _read_phen();
  _read_prod();
  _read_bvt();    /* 12/29/03 */

  _sxw_root_phen();

#ifdef TESTING
  _sxw_test();
  exit(0);
#endif

  //_recover_names();
}

/**
 * This function resets the model to default conditions
 */
void SXW_Reset(void) {
	char * temp;

	if (SW_Weather.use_markov) {
		free(SW_Markov.avg_ppt);
		free(SW_Markov.dryprob);
		free(SW_Markov.std_ppt);
		free(SW_Markov.wetprob);
	}

	SW_SIT_clear_layers();
	SW_WTH_clear_runavg_list();

	temp = strdup(SXW.f_watin);
	SW_CTL_init_model(temp);
	free(temp);
}

void SXW_InitPlot (void) {
/*======================================================*/
/* Call this from main::Plot_Init() after killing everything
 * so the sxw tables will be reset.
 */
#ifdef SXW_BYMAXSIZE
	GrpIndex g;
	RealF sizes[MAX_RGROUPS];
#endif

	_sxw_sw_clear_transp();
	_sxw_update_resource();

#ifdef SXW_BYMAXSIZE
	/* this stuff was taken from Run_Soilwat() but we're now trying
	 * to minimize the dynamic effect, so resources are always based
	 * on full-sized plants.  So we only need to do this at the
	 * beginning of each steppe-model iteration.
	 */
	ForEachGroup(g) sizes[g] = 1.0;
	_sxw_update_root_tables(sizes);
	_sxw_sw_setup(sizes);
#endif


}


void SXW_Run_SOILWAT (void) {
/*======================================================*/
/* need to update the resource vectors and set up the
 * ppt and temp (environs) here before calling
 * Env_Generate() and rgroup_Establish() in main().
 *
 * 2/28/2003 - cwb - adding changes mentioned at top of file.
 * 3/31/2003 - cwb - because we're always running soilwat to
 *             emulate full-size plants, computing roots etc
 *             gets done once during init plot.
 */

#ifndef SXW_BYMAXSIZE
	GrpIndex g;
	RealF sizes[MAX_RGROUPS];
	/* compute production values for transp based on current plant sizes */
	ForEachGroup(g)
		sizes[g] = RGroup[g]->relsize;
	_sxw_update_root_tables(sizes);// Issue #5
	_sxw_sw_setup(sizes);
#endif

	SXW.aet = 0.; /* used to be in sw_setup() but it needs clearing each run */

	//SXW_SW_Setup_Echo();
	_sxw_sw_run();

	/* now compute resource availability for the given plant sizes */
	_sxw_update_resource();

	/* and set environmental variables */
	_sxw_set_environs();

}

void SXW_SW_Setup_Echo(void) {
	char name[256] = {0};
	strcat(name, _debugout);
	FILE *f = OpenFile(strcat(name, ".input.out"), "a");
	int i;
	fprintf(f, "\n================== %d =============================\n", SW_Model.year);
	fprintf(f,"Fractions Grass:%f Shrub:%f Tree:%f Forb:%f BareGround:%f\n", SW_VegProd.fractionGrass, SW_VegProd.fractionShrub, SW_VegProd.fractionTree, SW_VegProd.fractionForb, SW_VegProd.fractionBareGround);
	fprintf(f,"Monthly Production Values\n");
	fprintf(f,"Grass\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.grass.litter[i],
				SW_VegProd.grass.biomass[i], SW_VegProd.grass.pct_live[i],
				SW_VegProd.grass.lai_conv[i]);
	}

	fprintf(f,"Shrub\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.shrub.litter[i],
				SW_VegProd.shrub.biomass[i], SW_VegProd.shrub.pct_live[i],
				SW_VegProd.shrub.lai_conv[i]);
	}

	fprintf(f,"Tree\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.tree.litter[i],
				SW_VegProd.tree.biomass[i], SW_VegProd.tree.pct_live[i],
				SW_VegProd.tree.lai_conv[i]);
	}

	fprintf(f,"Forb\n");
	fprintf(f,"Month\tLitter\tBiomass\tPLive\tLAI_conv\n");
	for (i = 0; i < 12; i++) {
		fprintf(f,"%u\t%f\t%f\t%f\t%f\n", i + 1, SW_VegProd.forb.litter[i],
				SW_VegProd.forb.biomass[i], SW_VegProd.forb.pct_live[i],
				SW_VegProd.forb.lai_conv[i]);
	}

	SW_SITE *s = &SW_Site;
	fprintf(f,"Soils Transp_coeff\n");
	fprintf(f,"Forb\tTree\tShrub\tGrass\n");
	ForEachSoilLayer(i)
	{// %u %u %u %u s->lyr[i]->my_transp_rgn_forb, s->lyr[i]->my_transp_rgn_tree, s->lyr[i]->my_transp_rgn_shrub, s->lyr[i]->my_transp_rgn_grass
		fprintf(f,"%6.2f %6.2f %6.2f %6.2f\n", s->lyr[i]->transp_coeff_forb, s->lyr[i]->transp_coeff_tree, s->lyr[i]->transp_coeff_shrub, s->lyr[i]->transp_coeff_grass);
	}

  // adding values to sxw structure for use in ST_stats.c
  SXW.grass_cover = SW_VegProd.grass.conv_stcr;
  SXW.shrub_cover = SW_VegProd.shrub.conv_stcr;
  SXW.tree_cover = SW_VegProd.tree.conv_stcr;
  SXW.forbs_cover = SW_VegProd.forb.conv_stcr;


	fprintf(f, "\n");
	CloseFile(&f);
}

//void SXW_SW_Output_Echo() {
//
//}

RealF SXW_GetPR( GrpIndex rg) {
/*======================================================*/
/* see _sxw_update_resource() for _resource_cur[]
This function is no longer utilized, SXW_GetTranspiration has replaced it
_resource_pr is no longer used as a parameter. We remain the code for the time being
KAP 7/20/2016
*/
	RealF pr = ZRO(_resource_pr[rg]) ? 0.0 : 1. / _resource_pr[rg];
	return pr;
	//return pr > 10 ? 10 : pr;
}

RealF SXW_GetTranspiration( GrpIndex rg) {
/*======================================================*/
/* see _sxw_update_resource() for _resource_cur[]
*/
	//printf("SXW_GetTranspiration _resource_cur[%d] = %.5f \n", rg, _resource_cur[rg]);
	return _resource_cur[rg];
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
			insertRootsXphen(_rootsXphen);
		}
		insertInputVars();
		insertInputProd();
		insertInputSoils();
		insertOutputVars(_resource_cur);
		insertRgroupInfo(_resource_cur);
		insertOutputProd(&SW_VegProd);
		insertRootsSum(_roots_active_sum);
		insertRootsRelative(_roots_active_rel);
		insertTranspiration();
		insertSWCBulk();
	}
}


static void  _read_files( void ) {
/*======================================================*/
  /* read list of input files.   */
  FILE *fin;
  int i, nfiles = SXW_NFILES;

  SXW.f_files = Parm_name(F_SXW);  /* aliased */
  MyFileName = SXW.f_files;
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
	char name[MAX_GROUPNAMELEN];
	FILE *fp;

	MyFileName = SXW.f_roots;
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
			_roots_max[Ilg(lyr, g)] = atof(p);
			lyr++;
		}
		if (lyr != SXW.NTrLyrs) {
			LogError(logfp, LOGFATAL,
					"%s: Group : %s : Missing layer values. Match up with soils.in file. Include zeros if necessary. Layers needed %u. Layers defined %u",
					MyFileName, name, SXW.NTrLyrs, lyr);
		}
	}

	if (cnt < Globals.grpCount) {
		LogError(logfp, LOGFATAL, "%s: Not enough valid groups found.",
				MyFileName);
	}

	CloseFile(&fp);
}

static void _read_phen(void) {
/*======================================================*/

  GrpIndex g;
  IntUS cnt=0;
  Months m;
  char *p;
  FILE *fp;

  MyFileName = SXW.f_phen;
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
      _phen[Igp(g,m)] = atof(p);
      m++;
    }

  }

  if (cnt < Globals.grpCount) {
    LogError(logfp, LOGFATAL,
             "%s: Not enough valid groups found.", MyFileName);
  }

  CloseFile(&fp);

}

static void _read_bvt(void) {
/*======================================================*/
/* read two numbers, biomass (g/m2) and
 * transpiration (cm/m2) for that biomass to construct a
 * simple linear relationship that gives g biomass / cm transp
 * per m2.
 */

  FILE *fp;
  RealF bmass, transp;

  MyFileName = SXW.f_bvt;
  fp = OpenFile(MyFileName,"r");

  GetALine(fp, inbuf);
  bmass = atof(inbuf);

  GetALine(fp, inbuf);
  transp = atof(inbuf);


  CloseFile(&fp);

  _bvt = bmass / transp;

}


static void _read_prod(void) {
	int x;
	FILE *fp;
	Months mon = Jan;

	GrpIndex g;
	IntUS cnt = 0;
	char *p;
	char * pch;
	MyFileName = SXW.f_prod;
	fp = OpenFile(MyFileName, "r");

	while (GetALine(fp, inbuf)) {
		x = sscanf(inbuf, "%lf", &_prod_litter[mon]);
		if (x < 1) {
			LogError(logfp, LOGFATAL, "%s: invalid record for litter %d.", MyFileName,
					mon + 1);
		}

		if (++mon > Dec)
			break;
	}

	if (mon <= Dec) {
		LogError(logfp, LOGWARN,
				"%s: No Veg Production values found after month %d", MyFileName,
				mon + 1);
	}

	GetALine(fp,inbuf); /* toss [end] keyword */
	mon = Jan;

	while (GetALine(fp, inbuf)) {
		pch = strstr(inbuf, "[end]");
		if(pch != NULL)
			break;
		p = strtok(inbuf, " \t"); /* g'teed to work via GetALine() */

		if ((g = RGroup_Name2Index(p)) < 0) {
			LogError(logfp, LOGFATAL, "%s: Invalid group name for biomass (%s) found.",
					MyFileName, p);
		}
		mon = Jan;
		while ((p = strtok(NULL, " \t"))) {
			if (mon > Dec) {
				LogError(logfp, LOGFATAL,
						"%s: More than 12 months of data found.", MyFileName);
			}
			_prod_bmass[Igp(g, mon)] = atof(p);
			mon++;
		}
		cnt++;
		if(cnt == SXW.NGrps)
			break;
	}

	if (cnt < Globals.grpCount) {
		LogError(logfp, LOGFATAL, "%s: Not enough valid groups found.",
				MyFileName);
	}


	GetALine(fp, inbuf); /* toss [end] keyword */
	mon = Jan;
	cnt=0;

	while (GetALine(fp, inbuf)) {
		pch = strstr(inbuf, "[end]");
		if (pch != NULL)
			break;
		p = strtok(inbuf, " \t"); /* g'teed to work via GetALine() */

		if ((g = RGroup_Name2Index(p)) < 0) {
			LogError(logfp, LOGFATAL,
					"%s: Invalid group name for pctlive (%s) found.",
					MyFileName, p);
		}
		mon = Jan;
		while ((p = strtok(NULL, " \t"))) {
			if (mon > Dec) {
				LogError(logfp, LOGFATAL,
						"%s: More than 12 months of data found.", MyFileName);
			}
			_prod_pctlive[Igp(g, mon)] = atof(p);
			mon++;
		}
		cnt++;
		if (cnt == SXW.NGrps)
			break;
	}

	if (cnt < Globals.grpCount) {
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

   MyFileName = SXW.f_watin;
   f = OpenFile(MyFileName, "r");


   while( GetALine(f, inbuf) ) {
     if (++lineno == (eOutput + 2)) {
       strcpy(_swOutDefName, DirName(SXW.f_watin));
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

static void _write_sw_outin(void) {
/*======================================================*/
/* make sure the outsetup file for soilwat contains only
 * the given information */
/* Note that there won't actually be any output.  These
 * keys are required to trigger the correct part of the
 * output accumulation routines.  Refer to the Output.c
 * module of SOILWAT for more.
 */
	FILE *fp;
	char pd[3];

	switch (SXW.NPds) {
	case MAX_WEEKS:
		strcpy(pd, "WK");
		break;
	case MAX_MONTHS:
		strcpy(pd, "MO");
		break;
	case MAX_DAYS:
		strcpy(pd, "DY");
		break;
	}

  // reading outsetup_v30.in file to get timestep values
  FILE *zq;
	zq = fopen(_swOutDefName, "r");
	int lineNum = 0;
	if(zq){
    char lineLen[100];
    char timestepVals[100];
    char *timeStep_final[4];
		while(fgets(lineLen, sizeof lineLen, zq) != NULL){
			if(lineNum == 5 || lineNum == 3){ // check both lines in case debug is turned off which would make line 3 the last line
        char *timeCheck = strstr(lineLen, "TIMESTEP"); // check if this line is the one that outlines the timesteps
        if(timeCheck){
          printf("%s\n", lineLen);
          strcpy(timestepVals, lineLen);
        }
      }
      lineNum++;
		}
    // done reading outsetup_v30.in file

    char ** res  = NULL;
    char *  p    = strtok (timestepVals, " ");
    int n_spaces = 0, i;

    /* split string and append tokens to 'res' */

    while (p) {
      res = realloc (res, sizeof (char*) * ++n_spaces);

      if (res == NULL)
        exit (-1); /* memory allocation failed */

      res[n_spaces-1] = p;

      p = strtok (NULL, " ");

    }

    /* realloc one extra element for the last NULL */

    res = realloc (res, sizeof (char*) * (n_spaces+1));
    res[n_spaces] = 0;

    /* print the result */
    for (i = 0; i < (n_spaces); ++i){
      printf ("res[%d] = %s\n", i, res[i]);

      if(strcmp(res[i],"dy")){
        //printf("DAY\n");
        timeStep_final[0] = "dy";
        //printf("timeStep_final[0]: %s\n", timeStep_final[0]);
      }
      else if(strcmp(res[i],"wk")){
        //printf("WEEK\n");
        timeStep_final[1] = "wk";
        //printf("timeStep_final[1]: %s\n", timeStep_final[1]);
      }
      else if(strcmp(res[i],"mo")){
        //printf("MO\n");
        timeStep_final[2] = "mo";
        //printf("timeStep_final[2]: %s\n", timeStep_final[2]);
      }
      else if(strcmp(res[i],"yr")){
        //printf("YR\n");
        timeStep_final[3] = "yr";
        //printf("timeStep_final[3]: %s\n", timeStep_final[3]);
      }
      else printf("\nelse\n\n");
    }


    /* free the memory allocated */

    free (res);
    //

		fclose(zq);
	}

	fp = OpenFile(_swOutDefName, "w");
	fprintf(fp, "TRANSP  SUM  %s  1  end  transp\n", pd);
	fprintf(fp, "PRECIP  SUM  YR  1  end  precip\n");
	fprintf(fp, "TEMP    AVG  YR  1  end  temp\n");
	if (SXW.debugfile) {
		fprintf(fp, "AET     SUM  YR  1  end  aet\n");
		fprintf(fp, "SWCBULK     FIN  MO  1  end  swc_bulk\n");
	}
  fprintf(fp, "TIMESTEP dy wk mo yr\n");

	CloseFile(&fp);
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
	if (SXW.debugfile || UseGrid)
		_make_swc_array();
}

static void _make_roots_arrays(void) {
/*======================================================*/
  int size;
  char *fstr = "_make_roots_array()";

  size  = SXW.NGrps * SXW.NTrLyrs;
  _roots_max     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

  size = SXW.NGrps * SXW.NPds * SXW.NTrLyrs;
  _rootsXphen       = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _roots_active     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _roots_active_rel = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

  //4 - Grass,Frob,Tree,Shrub
  size = 4 * SXW.NPds * SXW.NTrLyrs;
  _roots_active_sum = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
}

static void _make_phen_arrays(void) {
/*======================================================*/
  int size;
  char *fstr = "_make_phen_arrays()";

  size = SXW.NGrps * MAX_MONTHS;
  _phen = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

}

static void _make_prod_arrays(void) {
	int size;
	char *fstr = "_make_phen_arrays()";

	size = SXW.NGrps * MAX_MONTHS;
	_prod_bmass = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	_prod_pctlive = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
}

static void _make_transp_arrays(void) {
/*======================================================*/
/* SXW.transp holds year-to-year values and is populated in SOILWAT.
 * both are indexed by the macro Ilp().
 */
	char *fstr = "_make_transp_array()";
	int size;

	size = SXW.NPds * SXW.NSoLyrs;
	SXW.transpTotal = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	SXW.transpTrees = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	SXW.transpShrubs = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	SXW.transpForbs = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
	SXW.transpGrasses = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
}


static void _make_swc_array(void) {
/*======================================================*/
/* SXW.swc holds year-to-year values and is populated in SOILWAT.
 * it is indexed by the macro Ilp(). only used if soilwat option
 * specified with debugfile
 */
	char *fstr = "_make_swc_array()";
	int size = SXW.NPds * SXW.NSoLyrs;
	SXW.swc = (RealF *) Mem_Calloc(size, sizeof(RealF *), fstr);
}


/*static void _recover_names(void) {
======================================================
	int i;

	for (i = 0; i < SXW_NFILES; i++) {
		Mem_Free(*_sxwfiles[i]);
	}
}*/

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

	f = OpenFile(SXW.debugfile, "r");

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
	char vegProdNames[4][7];
	strcpy(vegProdNames[0], "TREE");
	strcpy(vegProdNames[1], "SHRUB");
	strcpy(vegProdNames[2], "GRASS");
	strcpy(vegProdNames[3], "FORB");
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
					fprintf(f, "\t%.4f", _rootsXphen[Iglp(r, t, p)]);
				fprintf(f, "\n");
			}
		}
	}


  /* sum actual total transpiration */
	ForEachTrPeriod(p)
	{
		//ForEachTranspLayer(t) sum += SXW.transp[Ilp(t,p)];
		for (t = 0; t < SXW.NSoLyrs; t++)
			sum += SXW.transpTotal[Ilp(t, p)];
	}

	fprintf(f, "\n================== %d =============================\n", SW_Model.year);
	fprintf(f, "MAP = %d(mm)\tMAT = %5.2f(C)\tAET = %5.4f(cm)\tAT = %5.4f(cm)\n\n", Env.ppt, Env.temp, SXW.aet, sum);

	fprintf(f, "Group     \tRelsize\tPR\tResource_cur\n");
	fprintf(f, "-----     \t-------\t-----\t-----\n");
	ForEachGroup(r) {
		sum1 += RGroup[r]->relsize;
		sum2 += RGroup[r]->pr;
		sum3 += _resource_cur[r];
		fprintf(f, "%s\t%.4f\t%.4f\t%.4f\n", RGroup[r]->name, RGroup[r]->relsize, RGroup[r]->pr, _resource_cur[r]);
	}
	fprintf(f, "-----     \t-------\t-----\t-----\n");
	fprintf(f, "%s\t%.4f\t%.4f\t%.4f\n", "Total", sum1, sum2, sum3);

	fprintf(f, "\n------ Production Values Daily Summed Across Types Monthly Averaged -------\n");
	fprintf(f, "Month\tBMass\tPctLive\tLAIlive\tVegCov\tTotAGB\n");
	fprintf(f, "-----\t-----\t-------\t-------\t------\t------\n");

	int doy = 1;
	ForEachMonth(p)
	{
		int days = 31, i;
		double pct_live = 0, lai_live = 0, vegcov = 0, total_agb = 0, biomass = 0;

		if (p == Apr || p == Jun || p == Sep || p == Nov) //all these months have 30 days
			days = 30;
		else if (p == Feb) { //February has either 28 or 29 days
			days = 28;
			if (Is_LeapYear(SW_Model.year))
				days = 29;
		} // all the other months have 31 days

		for (i = doy; i < (doy + days); i++) { //accumulating the monthly values...
			lai_live += (v->tree.lai_live_daily[i])
					+ (v->shrub.lai_live_daily[i])
					+ (v->grass.lai_live_daily[i])
					+ (v->forb.lai_live_daily[i]);
			vegcov += (v->tree.vegcov_daily[i]) + (v->shrub.vegcov_daily[i])
					+ (v->grass.vegcov_daily[i]) + (v->forb.vegcov_daily[i]);
			total_agb += (v->tree.total_agb_daily[i])
					+ (v->shrub.total_agb_daily[i])
					+ (v->grass.total_agb_daily[i])
					+ (v->forb.total_agb_daily[i]);
			pct_live += (v->tree.pct_live_daily[i])
					+ (v->shrub.pct_live_daily[i])
					+ (v->grass.pct_live_daily[i])
					+ (v->forb.pct_live_daily[i]);
			biomass += (v->tree.biomass_daily[i])
					+ (v->shrub.biomass_daily[i])
					+ (v->grass.biomass_daily[i])
					+ (v->forb.biomass_daily[i]);
		}
		doy += days; //updating the doy
		//biomass = (v->tree.biomass[p]) + (v->shrub.biomass[p])
		//		+ (v->grass.biomass[p]) + (v->forb.biomass[p]);		//getting the monthly biomass
		//pct_live = (v->tree.pct_live[p]) + (v->shrub.pct_live[p])
		//		+ (v->grass.pct_live[p]) + (v->forb.pct_live[p]);	//getting the monthly pct_live

		pct_live /= days;
		biomass /= days;
		lai_live /= days; //getting the monthly averages...
		vegcov /= days;
		total_agb /= days;

		fprintf(f, "%4d\t%.2f\t%.3f\t%.3f\t%.2f\t%.2f\n", p + 1, biomass,
				pct_live, lai_live, vegcov, total_agb);
	}

	for (t = 0; t < 4; t++) {
		fprintf(f, "\n------ Active Roots (sum) %s -------\n", vegProdNames[t]);
		fprintf(f, "Layer:");
		ForEachTrPeriod(p)
			fprintf(f, "\t%d", p + 1);
		fprintf(f, "\n");
		for (l = 0; l < SXW.NTrLyrs; l++) {
			fprintf(f, "%d", t + 1);
			ForEachTrPeriod(p)
				fprintf(f, "\t%.4f", _roots_active_sum[Itlp(t, l, p)]);
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
			//ForEachTranspLayer(t) {
			fprintf(f, "%d", t + 1);
			ForEachTrPeriod(p)
				fprintf(f, "\t%.4f", _roots_active_rel[Iglp(r, t, p)]);
			fprintf(f, "\n");
		}
	}

	fprintf(f, "\n------ Transpiration Total Values -------\nPeriod:");
	//ForEachTranspLayer(t)
	for (t = 0; t < SXW.NSoLyrs; t++)
		fprintf(f, "\t\tL%d", t + 1);
	fprintf(f, "\n");

	ForEachTrPeriod(p)
	{
		fprintf(f, "%d : ", p + 1);
		//ForEachTranspLayer(t)
		for (t = 0; t < SXW.NSoLyrs; t++)
			fprintf(f, "\t%.4f", SXW.transpTotal[Ilp(t, p)]);
		fprintf(f, "\n");
	}

	fprintf(f, "Current Soil Water Content:\n");
	ForEachTrPeriod(p)
	{
		fprintf(f, "%d : ", p + 1);
		ForEachSoilLayer(t)
			fprintf(f, "\t%5.4f", SXW.swc[Ilp(t, p)]);
		fprintf(f, "\n");
	}

	fprintf(f, "\n");
	CloseFile(&f);
}


#ifdef DEBUG_MEM
#include "myMemory.h"
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
   NoteMemoryRef(SXW.transp);


   SW_CTL_SetMemoryRefs();
}

#endif

 /***********************************************************/
 //returns the number of transpiration layers correctly for each veg_prod_type
int getNTranspLayers(int veg_prod_type) {
	if(veg_prod_type == 1)
		return SW_Site.n_transp_lyrs_tree;
	else if(veg_prod_type == 2)
		return SW_Site.n_transp_lyrs_shrub;
	else if(veg_prod_type == 3)
		return SW_Site.n_transp_lyrs_grass;
	else if(veg_prod_type == 4)
		return SW_Site.n_transp_lyrs_forb;
	return -1;
}

/***********************************************************/
void free_all_sxw_memory( void ) {
	free_sxw_memory();
	Mem_Free(SXW.f_roots);
	Mem_Free(SXW.f_phen);
	Mem_Free(SXW.f_bvt);
	Mem_Free(SXW.f_prod);
	Mem_Free(SXW.f_watin);
	Mem_Free(SXW.transpTotal);
	Mem_Free(SXW.transpTrees);
	Mem_Free(SXW.transpShrubs);
	Mem_Free(SXW.transpForbs);
	Mem_Free(SXW.transpGrasses);
	if (SXW.debugfile || UseGrid) Mem_Free(SXW.swc);
}

/***********************************************************/
void free_sxw_memory( void ) {
	Mem_Free(_roots_max);
	Mem_Free(_rootsXphen);
	Mem_Free(_roots_active);
	Mem_Free(_roots_active_rel);
	Mem_Free(_roots_active_sum);
	Mem_Free(_phen);
	Mem_Free(_prod_bmass);
	Mem_Free(_prod_pctlive);
}

/***********************************************************/
void load_sxw_memory( RealD* grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive ) {
	//load memory from the grid
	free_sxw_memory();
	_roots_max = Mem_Calloc(SXW.NGrps * SXW.NTrLyrs, sizeof(RealD), "load_sxw_memory()");
	_rootsXphen = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "load_sxw_memory()");
	_roots_active = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "load_sxw_memory()");
	_roots_active_rel = Mem_Calloc(SXW.NGrps * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "load_sxw_memory()");
	_roots_active_sum = Mem_Calloc(4 * SXW.NPds * SXW.NTrLyrs, sizeof(RealD), "load_sxw_memory()");
	_phen = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "load_sxw_memory()");
	_prod_bmass = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "load_sxw_memory()");
	_prod_pctlive = Mem_Calloc(SXW.NGrps * MAX_MONTHS, sizeof(RealD), "load_sxw_memory()");

	memcpy(_roots_max, grid_roots_max, SXW.NGrps * SXW.NTrLyrs * sizeof(RealD));
	memcpy(_rootsXphen, grid_rootsXphen, SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(_roots_active, grid_roots_active, SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(_roots_active_rel, grid_roots_active_rel, SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(_roots_active_sum, grid_roots_active_sum, 4 * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(_phen, grid_phen, SXW.NGrps * MAX_MONTHS * sizeof(RealD));
	memcpy(_prod_bmass, grid_prod_bmass, SXW.NGrps * MAX_MONTHS * sizeof(RealD));
	memcpy(_prod_pctlive, grid_prod_pctlive, SXW.NGrps * MAX_MONTHS * sizeof(RealD));
}

/***********************************************************/
void save_sxw_memory( RealD * grid_roots_max, RealD* grid_rootsXphen, RealD* grid_roots_active, RealD* grid_roots_active_rel, RealD* grid_roots_active_sum, RealD* grid_phen, RealD* grid_prod_bmass, RealD* grid_prod_pctlive ) {
	//save memory to the grid
	memcpy(grid_roots_max, _roots_max, SXW.NGrps * SXW.NTrLyrs * sizeof(RealD));
	memcpy(grid_rootsXphen, _rootsXphen, SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(grid_roots_active, _roots_active, SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(grid_roots_active_rel, _roots_active_rel, SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(grid_roots_active_sum, _roots_active_sum, 4 * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
	memcpy(grid_phen, _phen, SXW.NGrps * MAX_MONTHS * sizeof(RealD));
	memcpy(grid_prod_bmass, _prod_bmass, SXW.NGrps * MAX_MONTHS * sizeof(RealD));
	memcpy(grid_prod_pctlive, _prod_pctlive, SXW.NGrps * MAX_MONTHS * sizeof(RealD));
}
