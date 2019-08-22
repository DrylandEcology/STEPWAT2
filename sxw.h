/**
 * \file sxw.h
 * \brief Contains function declarations, global variables,
 *        etc required to interface STEPPE with SOILWAT.
 * 
 * The structs, functions, variables, and macros defined in this file
 * facilitate the conversion of SOILWAT output into STEPWAT input. 
 * 
 * \author
 *     Kyle Palmquist\n
 *     Daniel Schlaepfer\n
 *     Chandler Haukap\n
 *     Chris Bennett\n
 * 
 * \date 22 August 2019
 * 
 * \ingroup SXW
 */

#ifndef SXW_DEF
#define SXW_DEF

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

// The number of transpiration values retained by transp_data
#define MAX_WINDOW 100

/** 
 * \brief Stores statistics on transpiration over a window defined in inputs.
 * 
 * transp_data stores three arrays: ratios, transp, and sum of squares. These arrays form 
 * a moving window that stores "size" years worth of previous transpiration data.
 * average, ratio_average, and sum_of_sqrs all store summaries of their respective arrays
 * so that you do not have to iterate through the arrays every time to get information. 
 * size and add_here deal directly with manipulating the arrays.
 * 
 * \ingroup SXW
 */
struct transp_data {
  // average of all transp/ppt values currently inside ratios[]. It should be updated every time
  // a value is added to transp[].
  RealF ratio_average;

  // average of all transpiration currently inside transp[]. It should be updated every time a
  // value is added to ratios[].
  RealF average;

  // sum of all (xi - mean)^2 values currently inside SoS_array[]. It should be updated every
  // time a value is added to SoS_array[].
  RealF sum_of_sqrs;
  
  // The size off all arrays(see below) is MAX_WINDOW. How much of that window we use is determined
  // by size. size should be set before using a transp_data struct, and must be between 1 and MAX_WINDOW
  int size;

  // oldest_index has different behavior when the arrays are full versus when the arrays are partially full.
  // Partially full: oldest_index points to the next empty spot in the array. This means that oldest_index
  // should start at 0 and increment every time values are added to the arrays.
  // Full: oldest_index is the array index of the values that have been in the window the longest.
  // every time a value is overwritten oldest_index should be incremented by one. This behavior turns
  // the arrays into First in First out queues.
  int oldest_index;

  // ratios[] stores (transpiration/precipitation) values. It is used to keep track of
  // what value needs to be removed from the moving average.
  RealF ratios[MAX_WINDOW]; // transp/ppt
  
  // transp[] stores transpiration values. It is used to keep track of what value needs to be
  // removed from the moving ratio_average.
  RealF transp[MAX_WINDOW]; // transp

  //SoS_array[] stores the sum of squares values (xi - mean)^2 for the previous (size) iterations.
  RealF SoS_array[MAX_WINDOW]; //(xi - mean)^2
};

typedef struct transp_data transp_t;

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
