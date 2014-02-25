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
#include "SW_Control.h"
#include "SW_Model.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_Files.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
SXW_t SXW;

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_SOILWAT SW_Soilwat;

/*************** External Function Declarations ***************/
/***********************************************************/
void rgroup_ResPartIndiv(void);


/*************** Module/Local Variable Declarations ***************/
/***********************************************************/
/* these are initialized and maybe populated here but are used
 * in sxw_steppe.c so they aren't declared static.
 */
/* ----- 3d arrays ------- */
RealD * _roots_phen_lyr_grp, /* relative roots X phen by layer & group */
      * _roots_rel_phen_lyr_grp; /*relative to the total roots_phen_lyr_group */


/* ----- 2D arrays ------- */
/* malloc'ed here for consistency but only used */
/* in sxw_steppe.c and sxw_soilwat.c */

     /* rgroup by layer, ie, group-level values */
RealD * _roots_max,     /* read from root distr. file */
      * _roots_current, /* relsize * max root */
      * _roots_rel,     /* group's proportion  */
      * _roots_phen_totals, /* totals, rel roots X rel phen per group/layer */

     /* rgroup by period */
      * _phen,          /* phenology read from file */
      * _phen_grp_rel,      /* normalized to sum to 1.0 */
      * _transp_grp_totals; /* group's contribution to year's transp */

/* _transp_base is similar to SXW.transp because the function
 * _transp_contribution_by_group works for both.
 */
RealD *_transp_base; /* dyn array of nominal values read from file */


/* simple vectors hold the resource information for each group */
/* curr/equ gives the available/required ratio */
RealF _resource_equ[MAX_RGROUPS],  /* equilibrium resource utilization */
      _resource_cur[MAX_RGROUPS];  /* current resource utilization */

/* addition to meet changes specified at the top of the file */
RealF _Grp_BMass[MAX_RGROUPS];

/* and one 2D vector for the production constants */
RealF _prod_conv[MAX_MONTHS][3];

/* These are only used here so they are static.  */
static char inbuf[FILENAME_MAX];   /* reusable input buffer */
static char _swOutDefName[FILENAME_MAX];
static char *MyFileName;
static char **_files[SXW_NFILES];
static char _debugout[256];
static TimeInt _debugyrs[100], _debugyrs_cnt;


/*************** Local Function Declarations ***************/
/***********************************************************/

static void _read_files(void);
static void _read_times(void);
static void _read_roots_max(void);
static void _read_phen(void);
static void _read_transp(void);
static void _read_prod(void);
static void _read_watin(void);
static void _make_arrays(void);
static void _make_roots_arrays(void);
static void _make_phen_arrays(void);
static void _make_transp_arrays(void);
static void _write_sw_outin(void);
static void _recover_names(void);
static void _read_debugfile(void);
void _print_debuginfo(void);
static void _make_swc_array(void);

/****************** Begin Function Code ********************/
/***********************************************************/
void SXW_Init(void) {
  /* read SOILWAT's input files and initialize some variables */
  /* The shorthand table dimensions are set at the point
   * at which they become defined in _read_files().
   *
   * 2/14/03 - cwb - Added the summation for each group's biomass
   *     as per the notes at the top of the file.
   */
   GrpIndex rg; SppIndex sp;
   /* Sum each group's maximum biomass */
   ForEachGroup(rg) _Grp_BMass[rg] = 0.0;
   ForEachSpecies(sp)
     _Grp_BMass[Species[sp]->res_grp] += Species[sp]->mature_biomass;

   /* end code 2/14/03 */

   _files[0] = &SXW.f_times;
   _files[1] = &SXW.f_roots;
   _files[2] = &SXW.f_phen;
   _files[3] = &SXW.f_transp;
   _files[4] = &SXW.f_prod;
   _files[5] = &SXW.f_watin;

  SXW.NGrps = Globals.grpCount;
  _read_files();
  _read_times();
  _read_watin();
  if (*SXW.debugfile) _read_debugfile();
  _write_sw_outin();

  SW_CTL_init_model(SXW.f_watin);
  SXW.NTrLyrs = SW_Site.n_transp_lyrs;
  if (*SXW.debugfile) SXW.NSoLyrs = SW_Site.n_layers;

  _make_arrays();

  _read_transp();
  _read_roots_max();
  _read_phen();
  _read_prod();


/*  _recover_names(); */

  _sxw_relative_phen();
  _sxw_base_resource_vector();



#ifdef TESTING
  _sxw_test();
  exit(0);
#endif
}


void SXW_InitPlot (void) {
/*======================================================*/
/* Call this from main::Plot_Init() after killing everything
 * so the sxw tables will be reset.
 */
  GrpIndex g;
  RealF sizes[MAX_RGROUPS];

 _sxw_sw_clear_transp();
 _sxw_update_resource();


/* this stuff was taken from Run_Soilwat() but we're now trying
 * to minimize the dynamic effect, so resources are always based
 * on full-sized plants.  So we only need to do this at the
 * beginning of each steppe-model iteration.
 */
  /* 2/28/03: ForEachGroup(g) sizes[g] = RGroup[g]->relsize; */
  /* ForEachGroup(g) sizes[g] = RGroup[g]->relsize; */
  ForEachGroup(g) sizes[g] = 1.0;
  _sxw_update_root_tables(sizes);

  /* compute transpiration values based on current plant sizes */
  _sxw_sw_setup();

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


  SXW.aet = 0.;  /* used to be in sw_setup() but it needs clearing each run */
  _sxw_sw_run();

  /* now compute resource availability for the given plant sizes */
  _sxw_update_resource();

  /* and set environmental variables */
  _sxw_set_environs();


}

void SXW_GetResource( GrpIndex rg, RealF *baseline, RealF *actual) {
/*======================================================*/

  (*baseline) = _resource_equ[rg];
  (*actual)   = _resource_cur[rg];

}

void SXW_PrintDebug(void) {
/*======================================================*/
  TimeInt i;

  for(i=0; i<_debugyrs_cnt; i++) {
    if (SW_Model.year == _debugyrs[i]) {_print_debuginfo(); break;}
  }

}


static void  _read_files(void) {
/*======================================================*/
  /* read list of input files.   */
  FILE *fin;
  int i, nfiles = SXW_NFILES-1;  /* not counting first one */

  SXW.f_files = Parm_name(F_SXW);  /* aliased */
  MyFileName = SXW.f_files;
  fin = OpenFile(MyFileName,"r");

  for(i=0; i < nfiles; i++) {
    if (!GetALine(fin, inbuf)) break;
    *_files[i] = Str_Dup(Str_TrimLeftQ(inbuf));
  }

  if (i < nfiles) {
    LogError(logfp, LOGFATAL, "STEPWAT: %s: Insufficient files found",
            MyFileName);
  }

  CloseFile(&fin);


}


static void  _read_times(void) {
/*======================================================*/
/* transpiration and phenology periods are the same.    */

  FILE *fp;

  MyFileName = SXW.f_times;
  fp = OpenFile(MyFileName,"r");

  if (!GetALine(fp, inbuf)) {
    LogError( logfp, LOGFATAL, "%s: No data found!", MyFileName);
  }

  if (stricmp(inbuf, "week") == 0)
    SXW.NPds = MAX_WEEKS;
  else if (stricmp(inbuf,"month") == 0)
    SXW.NPds = MAX_MONTHS;
  else if (stricmp(inbuf,"day") == 0)
    SXW.NPds = MAX_DAYS;
  else {
    LogError(logfp, LOGFATAL, "%s: Invalid period", MyFileName);
  }

  CloseFile(&fp);

}


static void  _read_roots_max(void) {
/*======================================================*/
  GrpIndex g;
  int cnt=0, lyr;
  char *p;
  FILE *fp;

  MyFileName = SXW.f_roots;
  fp = OpenFile(MyFileName,"r");

  while( GetALine(fp, inbuf) ) {
    p = strtok(inbuf," \t"); /* g'teed to work via GetALine() */

    if ( (g=RGroup_Name2Index(p)) <0 ) {
      LogError(logfp, LOGFATAL, "%s: Invalid group name (%s) found.",
              MyFileName, p);
    }
    cnt++;
    lyr = 0;
    while ( (p=strtok(NULL," \t")) ) {
      _roots_max[Ilg(lyr,g)] = atof(p);
      lyr++;
    }

  }

  if (cnt < Globals.grpCount) {
    LogError(logfp, LOGFATAL,
             "%s: Not enough valid groups found.", MyFileName);
  }

  CloseFile(&fp);
}

static void _read_phen(void) {
/*======================================================*/

  GrpIndex g;
  IntUS cnt=0;
  TimeInt pd;
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
    pd = 0;
    while ((p=strtok(NULL," \t")) ) {
      if (pd == SXW.NPds) {
        LogError(logfp, LOGFATAL,
                 "%s: More than 12 months of data found.", MyFileName);
      }
      _phen[Igp(g,pd)] = atof(p);
      pd++;
    }

  }

  if (cnt < Globals.grpCount) {
    LogError(logfp, LOGFATAL,
             "%s: Not enough valid groups found.", MyFileName);
  }

  CloseFile(&fp);
}

static void _read_prod(void) {
/*======================================================*/
  int x;
  FILE *fp;
  Months mon=Jan;

  MyFileName = SXW.f_prod;
  fp = OpenFile(MyFileName,"r");

  while(GetALine(fp, inbuf) ) {
      x=sscanf( inbuf, "%f %f",
                      &_prod_conv[mon][PC_Bmass],
                      &_prod_conv[mon][PC_Litter]);
      if (x<2) {
         LogError(logfp, LOGFATAL, "%s: invalid record %d.",
                 MyFileName, mon +1);
      }

      if (++mon > Dec)
        break;
   }
   CloseFile(&fp);

   if (mon <= Dec ) {
     LogError(logfp, LOGWARN,
              "%s: No Veg Production values found after month %d",
              MyFileName, mon +1);
   }

}

static void _read_transp(void) {
/*======================================================*/
/* reads a table of baseline transpiration values.
 * checks to make sure the number of periods and the number
 * of layers match those expected by the SOILWAT, so
 * _read_times() must have been called for SXW.NPds and
 * SOILWAT must be initialized for the SW_Site.n_layers.
 *
/*  Input looks like this with no header line:
 * Period  Lyr1     Lyr2     Lyr3      Lyr4    Lyr5
   1,    0.0000,  0.0000,  0.0000,  0.0000,  0.0000
   2,    0.0000,  0.0000,  0.0000,  0.0000,  0.0000
   3,    0.2983,  0.1792,  0.0348,  0.0167,  0.0777

   Separators can be any combination of tab, space, comma
 */
  TimeInt cnt=0, pd;
  LyrIndex lyr;
  int i;
  char *p;
  FILE *fp;

  MyFileName = SXW.f_transp;
  fp = OpenFile(MyFileName,"r");

  while( GetALine(fp, inbuf) ) {
    p = strtok(inbuf,", \t"); /* g'teed to work via GetALine() */
    if (++cnt > SXW.NPds) break;

    i = atoi(p);
    if (i < 1) {
      LogError( logfp, LOGFATAL,
                "%s: Invalid layer number (%d)", MyFileName, i);
    }
    if ((pd = (TimeInt) i-1) == SXW.NPds) break;

    lyr = 0;
    while ((p=strtok(NULL," \t")) ) {
      if (lyr >= SW_Site.n_transp_lyrs) {
        LogError(logfp, LOGFATAL,
                        "%s: Too many layers of transpiration found.\n"
                        "Compare with soil layer definitions for SOILWAT.",
                        MyFileName);
      }
      _transp_base[Ilp(lyr,pd)] = atof(p);
      lyr++;
    }

  }

  if (cnt > SXW.NPds) {
    LogError(logfp, LOGFATAL, "%s: Too many time periods.  Expected %d.",
            MyFileName, SXW.NPds);
  }

  if (cnt < SXW.NPds) {
    LogError(logfp, LOGFATAL,
             "%s: Too few time periods.  Expected %d.",
             MyFileName, SXW.NPds);
  }

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
     if (++lineno == eOutput) {
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

  switch(SXW.NPds) {
    case MAX_WEEKS:  strcpy(pd, "WK"); break;
    case MAX_MONTHS: strcpy(pd, "MO"); break;
    case MAX_DAYS:   strcpy(pd, "DY"); break;
  }
  fp = OpenFile(_swOutDefName,"w");
  fprintf(fp,"TRANSP  SUM  %s  1  end  transp\n", pd);
  fprintf(fp,"PRECIP  SUM  YR  1  end  precip\n");
  fprintf(fp,"TEMP    AVG  YR  1  end  temp\n");
  if (*SXW.debugfile) {
    fprintf(fp,"AET     SUM  YR  1  end  aet\n");
    fprintf(fp,"SWC     FIN  MO  1  end  swc\n");
  }

  CloseFile(&fp);
}


static void _make_arrays(void) {
/*======================================================*/
/* central point to make all dynamically allocated arrays
 * now that the dimensions are known.
 */

  _make_transp_arrays();
  _make_roots_arrays();
  _make_phen_arrays();
  if (*SXW.debugfile)
    _make_swc_array();



}

static void _make_roots_arrays(void) {
/*======================================================*/
  int size;
  char *fstr = "_make_roots_array()";

  size  = SXW.NGrps * SXW.NTrLyrs ;
  _roots_max     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _roots_rel     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _roots_current = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

  size = SXW.NGrps * SXW.NPds * SXW.NTrLyrs;
  _roots_phen_lyr_grp     = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _roots_rel_phen_lyr_grp = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _roots_phen_totals      = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
}

static void _make_phen_arrays(void) {
/*======================================================*/
  int size;
  char *fstr = "_make_phen_array()";

  size = SXW.NGrps * SXW.NPds;
  _phen         = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);
  _phen_grp_rel = (RealD *) Mem_Calloc(size, sizeof(RealD), fstr);

}

static void _make_transp_arrays(void) {
/*======================================================*/
/* SXW.transp holds year-to-year values and is populated in SOILWAT.
 * _transp_base holds the "normal" values and is read from a file.
 * both are indexed by the macro Ilp().
 */
  char *fstr = "_make_transp_array()";
  int size;

   size = SXW.NGrps * SXW.NPds;
   _transp_grp_totals = (RealD *)
                         Mem_Calloc(size, sizeof(RealD), fstr);

  size = SXW.NPds * SXW.NTrLyrs;
  SXW.transp   = (RealD *) Mem_Calloc( size, sizeof(RealD), fstr);
  _transp_base = (RealD *) Mem_Calloc( size, sizeof(RealD), fstr);


}


static void _make_swc_array(void) {
/*======================================================*/
/* SXW.swc holds year-to-year values and is populated in SOILWAT.
 * it is indexed by the macro Ilp(). only used if soilwat option
 * specified with debugfile
 */
  char *fstr = "_make_swc_array()";
  int size = SXW.NPds * SXW.NSoLyrs;
  SXW.swc   = (RealF *) Mem_Calloc( size, sizeof(RealF *), fstr);


}


static void _recover_names(void) {
/*======================================================*/
  int i, last = SXW_NFILES-1;  /* recall we skipped the first file */
  for (i=0; i < last; i++) {
    Mem_Free(*_files[i]);
  }

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
   int cnt =0;
   TimeInt i;

   f = OpenFile(SXW.debugfile, "r");

   /* get name of output file */
   if (!GetALine(f, inbuf)) { CloseFile(&f); return; }
   strcpy(_debugout, inbuf);

   /* get output years */
   while( GetALine(f, inbuf) ) {
     _debugyrs[cnt++] = atoi(strtok( inbuf, " \t")); /* g'teed via getaline() */
     while (NULL != (date = strtok(NULL, " \t"))) {
       _debugyrs[cnt++] = atoi(date);
     }
   }
   _debugyrs_cnt = cnt;

   sprintf(errstr,"Debugging Transpiration turned on.\n%s will contain"
                  "%d years of output:\n", _debugout, _debugyrs_cnt);

   for(i=0; i<_debugyrs_cnt; i++) {
     sprintf(str, "%d\n", _debugyrs[i]);
     strcat(errstr, str);
   }
   strcat(errstr, "Note that data will always be appended,\n");
   strcat(errstr, "so clear file contents before re-use.\n");
   LogError(logfp,LOGNOTE,errstr);


   CloseFile(&f);

   /* now empty the file prior to the run */
   f = OpenFile(_debugout,"w");
   CloseFile(&f);

}


#include "SW_VegProd.h"
extern SW_VEGPROD SW_VegProd;

void _print_debuginfo(void) {
/*======================================================*/
  SW_VEGPROD *v = &SW_VegProd;
  TimeInt p;
  LyrIndex t;
  FILE *f;
  GrpIndex r;
  RealF y[MAX_LAYERS] = {0,0,0,0,0,0,0,0,0,0},
        sum = 0.;

  f = OpenFile(_debugout, "a");

  fprintf(f,"================== %d =============================\n",
          SW_Model.year);
  fprintf(f,"MAP = %d(mm)\tMAT = %5.2f(C)\tAET = %5.4f(cm)\n",
          Env.ppt, Env.temp, SXW.aet);

  fprintf(f,"Group     \tRelsize\tPR  \t\"avail\"\t\"required\"\n");
  ForEachGroup(r)
    fprintf(f,"%s\t%.4f\t%.4f\t%.4f\t%.4f\n",
            RGroup[r]->name, RGroup[r]->relsize, RGroup[r]->pr,
            RGroup[r]->res_avail, RGroup[r]->res_required);

  fprintf(f,"\n------ Production Values -------\n");
  fprintf(f,"Month\tBMass\tPctLive\tLAIlive\tVegCov\tTotAGB\n");
  ForEachMonth(p) {
    fprintf(f,"%4d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
            p+1,v->biomass[p], v->pct_live[p],v->lai_live[p],
            v->vegcov[p], v->total_agb[p]);
  }

  fprintf(f,"\n------ Transpiration Values -------\nLayer:");
  ForEachTreeTranspLayer(t)
    fprintf(f,"\t%d", t+1);
  fprintf(f,"\n");

  fprintf(f,"Base Transpiration:\n");
  ForEachTrPeriod(p) {
    fprintf(f, "%d : ", p+1);
    ForEachTreeTranspLayer(t)
      fprintf(f,"\t%.14f", _transp_base[Ilp(t,p)]);
    fprintf(f,"\n");
  }

  fprintf(f,"Current Transpiration:\n");
  fprintf(f,"Coeffs:");
  ForEachTreeTranspLayer(t) {
    ForEachGroup(r) y[t] += _roots_current[Ilg(t,r)];
    sum += y[t];
  }
  ForEachTreeTranspLayer(t)
    fprintf(f,"\t%.4f", y[t] / sum);
  fprintf(f,"\n");

  ForEachTrPeriod(p) {
    fprintf(f, "%d : ", p+1);
    ForEachTreeTranspLayer(t)
      fprintf(f,"\t%.14f", SXW.transp[Ilp(t,p)]);
    fprintf(f,"\n");
  }

  fprintf(f,"Current Soil Water Content:\n");
  ForEachTrPeriod(p) {
    fprintf(f, "%d : ", p+1);
    ForEachSoilLayer(t)
      fprintf(f,"\t%5.4f", SXW.swc[Ilp(t,p)]);
    fprintf(f,"\n");
  }

  fprintf(f,"\n");
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
     NoteMemoryRef(*_files[i]);
   }
   NoteMemoryRef(_roots_max);
   NoteMemoryRef(_roots_rel);
   NoteMemoryRef(_roots_current);
   NoteMemoryRef(_roots_phen_lyr_grp);
   NoteMemoryRef(_roots_rel_phen_lyr_grp);
   NoteMemoryRef(_roots_phen_totals);
   NoteMemoryRef(_phen);
   NoteMemoryRef(_phen_grp_rel);
   NoteMemoryRef(_transp_grp_totals);
   NoteMemoryRef(SXW.transp);
   NoteMemoryRef(_transp_base);


   SW_CTL_SetMemoryRefs();
}

#endif
