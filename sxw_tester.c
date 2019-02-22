/********************************************************/
/********************************************************/
/*  Source file: sxw_tester.c
 *
/*  Type: testing module
 *
/*  Purpose: Contains support code to test
 *           the sxw module.
 *
/*  Calls:
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *
/*  History:
/*     (25-Oct-2002) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "ST_steppe.h"
#include "SW_Defines.h"
#include "SW_Site.h"
#include "sxw_funcs.h"
#include "sxw_module.h"
#include "sxw_vars.h"


#define LOG fprintf
#define LOG_NL LOG(fp,"\n")

/************ External Variable Definitions  ***************/
/*              see ST_globals.h                       */
/***********************************************************/

extern SW_SITE SW_Site;

/* ----- 3d arrays ------- */
extern
  RealF * _rootsXphen,
        * _roots_active; /*relative to the total roots_phen_lyr_group */

extern
  RealF * _roots_rel,     /* group's proportion  */
        * _roots_phen_totals,

       /* rgroup by period */
        * _phen_grp_rel;      /* normalized to sum to 1.0 */


static void _print_results(void);
static void _read_test_data( void);

/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_test(void) {
/*======================================================*/

  _read_test_data();

  _sxw_update_resource();
  //_sxw_sw_setup();

  _print_results();

}

static void _print_results(void) {
/*======================================================*/

  LyrIndex lyr;
  TimeInt pd;
  GrpIndex g;
  RealF sum;
  FILE *fp;

  fp = OpenFile("testdata.txt","w");

  LOG(fp,"===============================================\n");
  LOG(fp,"Echoing relative sizes and transpiration values.\n\n");

  ForEachGroup(g) {
    LOG(fp,"%s(%4.3f)\t", RGroup[g]->name, getRGroupRelsize(g));
  }
  LOG_NL;

  ForEachTrPeriod(pd) {
    LOG(fp,"%d", pd);
    ForEachTreeTranspLayer(lyr)
      LOG(fp,"\t%5.4f", SXW.transpTotal[Ilp(lyr,pd)]);
    LOG_NL;
  }

  LOG(fp,"\n===============================================\n");
  LOG(fp,"     Relative rooting distributions.\n\n");
  LOG(fp,"Layer");
  ForEachGroup(g) LOG(fp,"\t%s", RGroup[g]->name);
  LOG_NL;
  ForEachTreeTranspLayer(lyr) {
    LOG(fp,"%d", lyr);
    ForEachGroup(g)
      LOG(fp, "\t%5.4f", _roots_rel[Ilg(lyr,g)]);
    LOG_NL;
  }

  LOG(fp,"\n===============================================\n");
  LOG(fp,"     Relative phenology.\n\n");
  LOG(fp,"Group");
  ForEachTrPeriod(pd) LOG(fp,"\t%d", pd);
  LOG_NL;
  ForEachGroup(g) {
    LOG(fp,"%s", RGroup[g]->name);
    ForEachTrPeriod(pd)
      LOG(fp,"\t%5.4f", _phen_grp_rel[Igp(g,pd)]);
    LOG_NL;
  }

  LOG(fp,"\n===============================================\n");
  LOG(fp,"  Check rel vals for rel phen X rel roots.\n\n");
  LOG_NL;
  ForEachGroup(g) {
    LOG(fp,"\n  ------  %s  -------\n", RGroup[g]->name);
    ForEachTrPeriod(pd) LOG(fp,"\t%d", pd);
    LOG_NL;
    ForEachTreeTranspLayer(lyr) {
      LOG(fp,"%d", lyr);
      ForEachTrPeriod(pd)
        LOG(fp,"\t%5.4f", _roots_active[Iglp(g,lyr,pd)]);
    LOG_NL;
    }
  }

  LOG(fp,"\n===============================================\n");
  LOG(fp,"  Totals for rel roots X rel phen.\n\n");
  LOG(fp,"Layer\\Mon");
  ForEachTrPeriod(pd) LOG(fp,"\t%d", pd);
  LOG_NL;
  ForEachTreeTranspLayer(lyr) {
    LOG(fp,"%d", lyr);
    ForEachTrPeriod(pd)
      LOG(fp,"\t%5.4f", _roots_phen_totals[Ilp(lyr,pd)]);
    LOG_NL;
  }

  LOG(fp,"\n===============================================\n");
  LOG(fp,"  Contribution to period's transpiration by group.\n\n");
  LOG(fp,"Group");
  ForEachTrPeriod(pd) LOG(fp,"\t%d", pd);
  LOG_NL;
  ForEachGroup(g) {
    LOG(fp,"%s", RGroup[g]->name);
    ForEachTrPeriod(pd) {
      sum = 0.;
      ForEachTreeTranspLayer(lyr) {
        sum += _roots_active[Iglp(g,lyr,pd)]
            *  SXW.transpTotal[Ilp(lyr,pd)];
      }
      LOG(fp,"\t%5.4f", sum);
    }
    LOG_NL;
  }

  CloseFile(&fp);
}

static void _read_test_data( void) {
/*======================================================*/
/* Read a set of arbitrary (non-baseline) group sizes and
 * transpiration values and populates the RGroup[g].relsizes
 * and the SXW.transp array. This simulates retrieval from
 * soilwat.
 */
  LyrIndex lyr, cnt;
  TimeInt pd;
  char *p;
  FILE *fp;
  char *fname = "sxw_testdata.in";

  fp = OpenFile(fname,"r");

  /* read the relative sizes */
  GetALine(fp, inbuf);
  p = strtok(inbuf," \t");

  /* read the transpiration values */
  cnt=0;
  while( GetALine(fp, inbuf) ) {
    p = strtok(inbuf,", \t"); /* g'teed to work via GetALine() */
    if (++cnt > SXW.NPds) break;

    pd = atoi(p) -1;
    lyr = 0;
    while ((p=strtok(NULL," \t")) ) {
      if (lyr > SW_Site.n_transp_lyrs_tree) {
        LogError(logfp, LOGFATAL, "Too many layers of transpiration found in %s\n"
                        "Compare with soil layer definitions for SOILWAT.",
                        fname);
      }
      SXW.transpTotal[Ilp(lyr,pd)] = atof(p);
      lyr++;
    }

  }

  if (cnt > SXW.NPds) {
    LogError(logfp, LOGFATAL, "Too many time periods found in %s.  Expected %d.",
            fname, SXW.NPds);
  }

  if (cnt < SXW.NPds) {
    LogError(logfp, LOGFATAL, "Not enough time periods found in %s.  Expected %d.",
            fname, SXW.NPds);
  }

  CloseFile(&fp);
}
