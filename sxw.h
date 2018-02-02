/********************************************************/
/********************************************************/
/*  Source file: sxw.h
/*  Type: header
/*  Purpose: Contains pertinent declarations, global variables,
 *           etc required to support new functions that
 *           interface STEPPE with SOILWAT.
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
/*  History:
/*     (14-Apr-2002) -- INITIAL CODING - cwb
/*
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

int getNTranspLayers(int veg_prod_type);

struct stepwat_st {
  RealD *transpTotal; /* points to dynamic array indexed by Ilp() */
  RealD *transpTrees;
  RealD *transpShrubs;
  RealD *transpForbs;
  RealD *transpGrasses;

  RealD *transpTotal_avg,
        *transpTrees_avg,
        *transpShrubs_avg,
        *transpForbs_avg,
        *transpGrasses_avg;

  RealF  temp,   /* soilwat's MAT */
         ppt;    /* soilwat's MAP */
  TimeInt NPds;  /* number of transp periods= maxdays, maxweeks, maxmonths */
  IntUS NTrLyrs, /* # transp. layers taken from SOILWAT */
        NGrps;   /* # plant groups taken from STEPPE */
  IntUS NSoLyrs;  /* number of soil layers defined */

  /* These are file names */
  char  *f_files,  /* list of input files for sxw */
        *f_roots,  /* root distributions */
        *f_phen,   /* phenology */
        *f_bvt,    /* biomass vs transpiration 12/29/03 */
        *f_prod,   /* biomass to prod. conv. nos. */
        *f_watin;  /* soilwat's input file */

  /* DEBUG stuff */
  char *debugfile; /* added in ST_Main(), read to get debug instructions */
  RealF *swc, /* dynamic array(Ilp) of SWC from SOILWAT */
         aet;
              /* soilwat's evapotranspiration for the year */
  RealD  surfaceTemp;   /* soilwat's surfaceTemp */

  // PPT variables
  int    yearInterval; // keep track of years
  int    curMonth;

  RealF PPT_sum,
        PPT_rain,
        PPT_snow_fall,
        PPT_snow_melt,
        PPT_snow_loss;

  RealF *SWA_grass_avg, // 2D array to store SWA vals ([days of year][number of max layers])
        *SWA_shrub_avg,
        *SWA_tree_avg,
        *SWA_forb_avg;

  RealF *SWA_master, // 4D array to store SWA for all veg_types
        *dSWAbulk, // 4D array to store actual available SWA
        *dSWA_repartitioned; // 4D array to store repartioned SWA values

  RealF *sum_dSWA_repartitioned;

  RealF transp_SWA[MAX_YEARS][11]; // store the sum of SWA and transp for each year and resource. transp_SWA[year][steppe_resource_group]
};

struct soilwat_average{
  RealF *soilinfilt_avg,
        *runoff_total_avg,
        *surface_runoff_avg,
        *surface_runon_avg,
        *runoff_snow_avg,
        *vwcbulk_avg,
        *vwcmatric_avg,
        *swamatric_avg,
        *swabulk_avg,
        *swpmatric_avg,
        *surfacewater_avg,
        *evapsoil_avg,
        *evapsurface_total_avg,
        *evapsurface_tree_avg,
        *evapsurface_shrub_avg,
        *evapsurface_forb_avg,
        *evapsurface_grass_avg,
        *evapsurface_litter_avg,
        *evapsurface_water_avg,

        *interception_total_avg,
        *interception_tree_avg,
        *interception_shrub_avg,
        *interception_forb_avg,
        *interception_grass_avg,
        *interception_litter_avg,

        *lyrdrain_avg,

        *hydred_total_avg,
        *hydred_tree_avg,
        *hydred_shrub_avg,
        *hydred_forb_avg,
        *hydred_grass_avg,

        *pet_avg,
        *wetday_avg,
        *snowpack_water_eqv_avg,
        *snowpack_depth_avg,
        *deepswc_avg,
        *soiltemp_avg,
        *estab_avg,
        *max_temp_avg,
        *min_temp_avg,
        *avg_temp_avg,
        *aet_avg,
        *val_snowmelt_avg,
        *val_snowloss_avg;

    RealD *surfaceTemp_avg;

    // Carbon Variables
    RealD *biomass_grass_avg,
          *biomass_shrub_avg,
          *biomass_tree_avg,
      		*biomass_forb_avg,
      		*biomass_total_avg,

      		*biolive_grass_avg,
      		*biolive_shrub_avg,
      		*biolive_tree_avg,
      		*biolive_forb_avg,
      		*biolive_total_avg,

      		*bio_mult_grass_avg,
      		*bio_mult_shrub_avg,
      		*bio_mult_tree_avg,
      		*bio_mult_forb_avg,

      		*wue_mult_grass_avg,
      		*wue_mult_shrub_avg,
      		*wue_mult_tree_avg,
      		*wue_mult_forb_avg;

    RealF *swc_avg;
};

#define SXW_NFILES 5

typedef struct stepwat_st SXW_t;
typedef struct soilwat_average SXW_avg;

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
#define Iylp(y,l,p,x) (((y)*Globals.runModelYears * SXW.NTrLyrs * SXW.NPds * 2) + ((l)*SXW.NTrLyrs * SXW.NPds * 2) + ((p)*SXW.NPds * 2) + ((x)*2))

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
