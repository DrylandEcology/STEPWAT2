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
  RealD *transpTotal_avg; /* points to dynamic array indexed by Ilp() */
  RealD *transpTrees_avg;
  RealD *transpShrubs_avg;
  RealD *transpForbs_avg;
  RealD *transpGrasses_avg;

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

  /* % Cover from prod.in */
  /*float grass_cover,
        shrub_cover,
        tree_cover,
        forbs_cover;*/
  RealD critSoilWater[4]; // storing values in same order as defined in rgroup.in (0=tree, 1=shrub, 2=grass, 3=forb)

  /* DEBUG stuff */
  char *debugfile; /* added in ST_Main(), read to get debug instructions */
  RealF *swc, /* dynamic array(Ilp) of SWC from SOILWAT */
         aet;     /* soilwat's evapotranspiration for the year */
  RealD  surfaceTemp;   /* soilwat's surfaceTemp */

  // PPT variables
  int    yearInterval; // keep track of years
  int    curMonth;
  RealF  PPTVal[500]; // array to store ppt vals
  RealF  PPT_day[500];
  RealF  PPT_week[500];
  RealF  PPT_month[500];

  // store converted SWA values
  float SWAbulk_grass[366][25], // 2D array to store SWA vals ([days of year][number of max layers])
        SWAbulk_shrub[366][25],
        SWAbulk_tree[366][25],
        SWAbulk_forb[366][25],
        SWAbulk_grass_avg[MAX_YEARS][25], // 2D array to store SWA vals ([days of year][number of max layers])
        SWAbulk_shrub_avg[MAX_YEARS][25],
        SWAbulk_tree_avg[MAX_YEARS][25],
        SWAbulk_forb_avg[MAX_YEARS][25];

  RealF transp_SWA[MAX_YEARS][11]; // store the sum of SWA and transp for each year and resource. transp_SWA[year][steppe_resource_group]


  // 2D array to store 4 critical values per layer
  float SWCbulk[4][8]; // TODO: first value needs to be (number of layers * plant types) - not hardcoded
  float SWCoriginal[500][8]; // storing SWC values here instead of in *swc since that does not have enough storage for more than month timestep

  int curInterval;
  int tempInt;

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
