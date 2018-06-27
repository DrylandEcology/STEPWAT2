/********************************************************/
/********************************************************/
/*  Source file: sxw.h
 *  Type: header
 *  Purpose: Contains pertinent declarations, global variables,
 *           etc required to support new functions that
 *           interface STEPPE with SOILWAT.
 *  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model. */
/*  History:
 *     (14-Apr-2002) -- INITIAL CODING - cwb */
/********************************************************/
/********************************************************/

#ifndef SXW_DEF
#define SXW_DEF

/* comment the next line to use STEPPE-provided values
 * of biomass and size when computing SOILWAT parameters
 * or scaling transpiration.  If defined, the macro
 * sets the sizes to the maximum.
 */
/*#define SXW_BYMAXSIZE*/


#include "generic.h"
#include "SW_Times.h"
#include "ST_defines.h"
#include "sw_src/SW_Defines.h"

int getNTranspLayers(int veg_prod_type);
void free_all_sxw_memory( void );

struct stepwat_st {
  // ------ Values from SOILWAT2:
  // Note: the function `SW_OUT_set_SXWrequests` specifies the required
  // output time periods and output aggregation types

  // transpXXX: monthly sum of soilwat's transpiration by soil layer
  // * these are dynamic arrays that are indexed by Ilp()
  RealD *transpTotal, // total transpiration, i.e., sum across vegetation types
        *transpVeg[NVEGTYPES]; // transpiration as contributed by vegetation types
  RealF *swc; // monthly mean SWCbulk for each soil layer

  // fixed monthly array:
  RealF ppt_monthly[MAX_MONTHS];  // monthly sum of soilwat's precipitation
  RealF temp_monthly[MAX_MONTHS];  // monthly mean soilwat's air temperature

  // annual values:
  RealF temp,   // annual mean soilwat's air temperature
        ppt,    // annual sum of soilwat's precipitation
        aet;    // annual sum of soilwat's evapotranspiration


  // ------ Resource values partitioned by STEPWAT2:
  // current years 'resources' partitioned to each STEPWAT resource group:
  RealF transp_SWA[MAX_RGROUPS];

  /* `SWA`-related resource functionality is currently not implemented because
    the former implementation was incorrect (#133 and #138),
    please, see https://github.com/DrylandEcology/STEPWAT2/issues/136 for
    details.

  RealF *sum_dSWA_repartitioned; // required for `_SWA_contribution_by_group`
  */


  // ------ Size variables:
  TimeInt NPds;  /* number of transp periods= maxdays, maxweeks, maxmonths */
  IntUS NTrLyrs, /* # transp. layers taken from SOILWAT */
        NGrps;   /* # plant groups taken from STEPPE */
  IntUS NSoLyrs;  /* number of soil layers defined */

  // ------ These are file names:
  char  *f_files,  /* list of input files for sxw */
        *f_roots,  /* root distributions */
        *f_phen,   /* phenology */
        *f_bvt,    /* biomass vs transpiration 12/29/03 */
        *f_prod,   /* biomass to prod. conv. nos. */
        *f_watin;  /* soilwat's input file */

  // ------ DEBUG stuff:
  char *debugfile; /* added in ST_Main(), read to get debug instructions */
};

#define SXW_NFILES 5

typedef struct stepwat_st SXW_t;

#define ForEachTrPeriod(i) for((i)=0; (i)< SXW.NPds; (i)++)


/* convert 3-d index to actual array index for
   group/layer/phenology 3d table */
#define Iglp(g,l,p) (((g)*SXW.NTrLyrs*SXW.NPds) + ((l)*SXW.NPds) + (p))

/* convert 3-d index to actual array index for
 * veg-prod-type/layer/phenology
 */
#define Itlp(t,l,p) (((t)*SXW.NTrLyrs*SXW.NPds) + ((l)*SXW.NPds) + (p))

/* convert 4-d index to actual array index for
 * veg-prod-type/crit-value/layer/phenology
 */
 //veg_type, new_critical_value, layer, timeperiod
#define Itclp(t,c,l,p) (((t)*SXW.NTrLyrs*SXW.NPds) + ((c)*NVEGTYPES) + ((l)*SXW.NPds) + (p)) // c*4 is because there are 4 critical values

// for use with avg values
// year, layer, timeperiod, avg/std
#define Iylp(y,l,p,pd,x) (((y)*Globals.runModelYears * SXW.NTrLyrs * 4 * SXW.NPds * 2) + ((l)*SXW.NTrLyrs * 4 * SXW.NPds * 2) + ((p) * 4 * SXW.NPds * 2) + ((pd)*4*2) + ((x) * 2))

// for soilwat average and standard deviation
// year, timeperiod, choice (avg or std), timeperiod (dy, wk, mo, yr)
// difference between p and b is p is current period within period (ie. for day it could be 0 to 364) and b is just timeperiod (0 to 3 where 0 is day and 3 is yr)
#define Iypc(y,p,c,b) (((y)*Globals.runModelYears * SXW.NPds * NVEGTYPES) + ((p)*SXW.NPds * NVEGTYPES) + ((c) * NVEGTYPES) + (b))

// veg type, layer, timeperiod
#define Ivlp(v,l,p) (((v)*NVEGTYPES * SXW.NTrLyrs * SXW.NPds) + ((l)*SXW.NTrLyrs * SXW.NPds) + ((p)*SXW.NPds))

/* convert 2d layer by period indices to
  layer/phenology 1D index */
#define Ilp(l,p) ((l)*SXW.NPds + (p))

/* convert 2d group by period indices to
   group/phenology 1D index */
#define Igp(g,p) ((g)*SXW.NPds + (p))

/* convert 2d group by layer indices to
   layer/period 1D index */
#define Ilg(l,g) ((l)*SXW.NGrps + (g))

#endif
