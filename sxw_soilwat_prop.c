/********************************************************/
/********************************************************/
/*  Source file: sxw_soilwat.c
/*  Type: module
 *
/*  Purpose: Subroutines to handle interactions with the
 *           SOILWAT model including setting up input/output
 *           files, writing to internal data structures,
 *           and actually running the model.
 *
/*  Called by: sxw.c
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
/*  History:
/*     (14-Apr-2002) -- INITIAL CODING - cwb
 *
 *     28-Feb-02 - cwb - The model runs but plants die
 *         soon after establishment in a way that suggests
 *         chronic stretching of resources.  At this time
 *         I'm setting the input to soilwat to be based on
 *         full-sized plants to provide maximum typical
 *         transpiration and interpreting that as "available"
 *         resource. The affected routines are all the ones
 *         called within sxw_sw_setup():
 *           _update_productivity();
 *           _update_pct_live();
 *         See also sxw.c.
 *     19-Jun-2003 - cwb - Added RealD (double precision)
 *                 types for internal dynamic matrices and
 *                 other affected variables.  See notes in
 *                 sxw.c.
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "Times.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "SW_Defines.h"
#include "sxw.h"
#include "sxw_module.h"
#include "SW_Control.h"
#include "SW_Model.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Files.h"


/*************** Global Variable Declarations ***************/
/***********************************************************/
#include "sxw_vars.h"

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_VEGPROD SW_VegProd;

/*********** Local/Module Variable Declarations ************/
/***********************************************************/
extern RealD *_roots_current,
             *_phen;          /* phenology read from file */
       RealF _prod_conv[MAX_MONTHS][3];

extern RealF _Grp_BMass[];  /* added 2/28/03 */

/*************** Local Function Declarations ***************/
/***********************************************************/
static void _update_transp_coeff(void);
static void _update_productivity(void);
static void _update_pct_live(void);


/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_sw_setup (void) {
/*======================================================*/
  _update_transp_coeff();
  _update_productivity();
  _update_pct_live();


}

void _sxw_sw_run(void) {
/*======================================================*/

  SW_Model.year = SW_Model.startyr + Globals.currYear -1;
  SW_CTL_run_current_year();

}

void _sxw_sw_clear_transp(void) {
/*======================================================*/

  Mem_Set(SXW.transp, 0, SXW.NPds * SXW.NTrLyrs * sizeof(RealD));

}


static void _update_transp_coeff(void) {
/*======================================================*/
/* copy the relative root distribution to soilwat's layers */
/* POTENTIAL BUG:  if _roots_current is all zero but for some
 *   reason the productivity values (esp. %live) are >0,
 *   there can be transpiration but nowhere for it to come
 *   from in the soilwat model.
 */
  SW_LAYER_INFO *y ;
  GrpIndex g;
  LyrIndex t;
  RealF sum=0.;

#if (0)
  /* reset transp coeff sums */
  ForEachTranspRegion(t) SW_Site.sum_transp_coeff[t] = 0.;
#endif

  ForEachTreeTranspLayer(t) { y = SW_Site.lyr[t];
    y->transp_coeff = 0.;
    ForEachGroup(g)
      y->transp_coeff += (RealF)_roots_current[Ilg(t,g)];
    sum += y->transp_coeff;
  }

  /* normalize coefficients to 1.0 */
  ForEachTreeTranspLayer(t) { y = SW_Site.lyr[t];
    y->transp_coeff /= sum;
#if(0)
    SW_Site.sum_transp_coeff[y->my_transp_rgn]
       += y->transp_coeff;
#endif
  }
}

static void _update_productivity(void) {
/*======================================================*/
/* must convert STEPPE's plot-level biomass to SOILWAT's
 * sq.meter - based biomass.
 *
 * 2/28/03 - cwb - as per the notes at the top of this file,
 *     RGroup_GetBiomass is being changed to
 *     SUM(RGroup[g]->mature_biomass)
 */
  GrpIndex g;
  Months m;
  SW_VEGPROD *v = &SW_VegProd;
  RealF max = 0.0, bmass = 0.0;

  ForEachGroup(g) bmass += _Grp_BMass[g] / Globals.plotsize *.5;
  ForEachMonth(m) {
    v->biomass[m] = bmass * _prod_conv[m][PC_Bmass];
    max = fmax(max, v->biomass[m]);
  }

  ForEachMonth(m) v->litter[m] = max * _prod_conv[m][PC_Litter];

  SW_VPD_init();

#if 0
/* This chunk goes in the foreachmonth() loop, if needed */
    ForEachGroup(g) {
/*2/28/03  v->biomass[m] += (RGroup_GetBiomass(g) *_prod_conv[m][PC_Bmass])

      /* v->biomass[m] += RGroup_GetBiomass(g) * RGroup[g]->max_per_sqm
                        * _prod_conv[m][PC_Bmass]; */
    }
#endif
}

static void _update_pct_live(void) {
/*======================================================*/
/* 2/28/03 - cwb - adding changes as per notes at top of file.
 */
  RealF *sums;
  TimeInt p, i;
  Months m;
  /* 3/25/03 - it seems that the earlier method of calculating %live
     isn't very similar to some values found in old prod files,
     so I'm just doing a simple substitution here for the time, until
     I decide which is the best way to go.  This only works for
     Npds = 12.
     */
  RealF pcts[12] = {0., 0., .1, .4, .7, .8, .9, .5, .2, .1, .0, .0};

  sums = (RealF *) Mem_Calloc(SXW.NPds, sizeof(RealF),
                              "_sxw_phen_rel()");

  ForEachTrPeriod(p) sums[p] = pcts[p];

  /* pct_live is always in months, so if the number of
   * transpiration/phenology periods is different
   * (ie greater) than months, we have to sum the values
   */
  ForEachMonth(m) SW_VegProd.pct_live[m] = 0.;
  switch (SXW.NPds) {
    case MAX_DAYS:
         for( i=1; i <= 365; i++ ) {
           m = doy2month(i);
           SW_VegProd.pct_live[m] += sums[i];
         }
         break;

    case MAX_WEEKS:
         for(i=0; i < 52; i++) {  /* 52/4=13 > MAX_MONTHS */
           m = (Months) (i / 4);
           SW_VegProd.pct_live[m] += sums[i];
         }
         break;

    case MAX_MONTHS:
         ForEachMonth(m)
           SW_VegProd.pct_live[m] = sums[m];

         break;

    default:
      LogError(logfp,LOGFATAL,"PGMR: Invalid NPDS in update_pct_live()");
  }


  Mem_Free(sums);

}
