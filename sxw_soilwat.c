/**
 * \file sxw_soilwat.c
 * \brief Handles function calls to [SOILWAT2](\ref sw_src)
 * 
 * Functions in this file set up SOILWAT2 to be run inside of STEPWAT2.
 * This means allocating memory and running SOILWAT.
 * 
 *  History:
 *     (14-Apr-2002) -- INITIAL CODING - cwb
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
 *	08/01/2012 - DLM - updated _update_productivity() function
 *          to use the 3 different VegProds now used in soilwat 
 * 
 * \author CWB (initial programming)
 * \date 14 April 2002
 * \author DLM
 * \date 1 August 2012
 * \ingroup SXW_PRIVATE
 */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */


#include "sw_src/include/generic.h"
#include "sw_src/include/filefuncs.h"
#include "sw_src/include/myMemory.h"
#include "sw_src/include/Times.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/include/SW_Defines.h"
#include "sxw.h" // externs `*SXWResources`
#include "sxw_module.h"
#include "sxw_vars.h"
#include "sw_src/include/SW_Control.h"
#include "sw_src/include/SW_Weather.h"
#include "sw_src/include/SW_Model.h"
#include "sw_src/include/SW_Site.h"
#include "sw_src/include/SW_SoilWater.h"
#include "sw_src/include/SW_VegProd.h"
#include "sw_src/include/SW_Files.h"
#include "sw_src/include/SW_Sky.h"



/*************** Local Function Declarations ***************/
/***********************************************************/
static void _update_transp_coeff(void);
static void _update_productivity(RealF size[]);


/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_sw_setup (RealF sizes[]) {
/*======================================================*/
	int doy, k;
	SW_VEGPROD *v = &SoilWatAll.VegProd;

  _update_productivity(sizes);
  _update_transp_coeff();

  for (doy = 1; doy <= MAX_DAYS; doy++) {
    ForEachVegType(k) {
      v->veg[k].litter_daily[doy] = 0.;
      v->veg[k].biomass_daily[doy] = 0.;
      v->veg[k].pct_live_daily[doy] = 0.;
      v->veg[k].lai_conv_daily[doy] = 0.;
    }
  }
}

/** @brief Handle the weather generator to create new daily weather for current year

  This function takes the place of SOILWAT2's `SW_WTH_read()` and
  `SW_WTH_finalize_all_weather()`.

  Note: `Env_Generate()` (via `_sxw_sw_run()`) sets SOILWAT2's "year"
  (`SW_Model.year`) to `SW_Model.startyr + Globals->currYear - 1`.
  SOILWAT2 expects current year's weather values at element
  `SW_Model.year - SW_Weather.startYear` (see `SW_WTH_new_day()`).

  Note: As an alternative to the function here that is called by `Env_Generate()`,
  STEPWAT2 could generate weather for all years at once if
  each grid cell would store a local copy of `SW_Weather`.
*/
void _sxw_generate_weather(void) {
  SW_WEATHER *w = &SoilWatAll.Weather;
  SW_SKY *sky = &SoilWatAll.Sky;

  deallocateAllWeather(w->allHist, w->n_years);
  w->n_years = 1;
  w->startYear = SoilWatAll.Model.startyr + Globals->currYear - 1;
  allocateAllWeather(&w->allHist, w->n_years, &LogInfo);

  if (!w->use_weathergenerator_only) {
    LogError(
      &LogInfo,
      LOGERROR,
      "STEPWAT2 expects 'use_weathergenerator_only'."
    );
  }

  // Make sure monthly flags are set to interpolate monthly values into daily values
  w->use_humidityMonthly = swTRUE;
  w->use_cloudCoverMonthly = swTRUE;
  w->use_windSpeedMonthly = swTRUE;

  readAllWeather(
    w->allHist,
    w->startYear,
    w->n_years,
    swTRUE, // `use_weathergenerator_only`
    w->name_prefix, // not used because `use_weathergenerator_only`
	  w->use_cloudCoverMonthly,
    w->use_humidityMonthly,
    w->use_windSpeedMonthly,
    w->n_input_forcings,
    w->dailyInputIndices,
    w->dailyInputFlags,
    sky->cloudcov,
    sky->windspeed,
    sky->r_humidity,
    SoilWatAll.Model.cum_monthdays,
    SoilWatAll.Model.days_in_month,
    &LogInfo
    );

  finalizeAllWeather(&SoilWatAll.Markov, w, SoilWatAll.Model.cum_monthdays,
                      SoilWatAll.Model.days_in_month, &LogInfo); // run the weather
}


void _sxw_sw_run(void) {
/*======================================================*/
  SW_GEN_OUT *swGenOut = &SoilWatAll.GenOutput;
   LyrIndex lyrno;
   TimeInt month;
   int vegType;

 	SoilWatAll.Model.year = SoilWatAll.Model.startyr + Globals->currYear-1;

   // Copy global values to SOILWAT2's SW_ALL to work with correct values
   SoilWatAll.Model.runModelIterations = SuperGlobals.runModelIterations;
   SoilWatAll.Model.runModelYears = SuperGlobals.runModelYears;

 	SW_CTL_run_current_year(&SoilWatAll, SoilWatOutputPtrs, &LogInfo);

   // Copy the results of SOILWAT2's SXW values
   for(lyrno = 0; lyrno < SXW->NSoLyrs; lyrno++) {
     ForEachTrPeriod(month) {
       SXW->swc[Ilp(lyrno, month)] = swGenOut->swc[lyrno][month];
       SXW->transpTotal[Ilp(lyrno, month)] = swGenOut->transpTotal[lyrno][month];

       SXW->ppt_monthly[month] = swGenOut->ppt_monthly[month];
       SXW->temp_monthly[month] = swGenOut->temp_monthly[month];

       ForEachVegType(vegType) {
         SXW->transpVeg[vegType][Ilp(lyrno, month)] =
                                   swGenOut->transpVeg[vegType][lyrno][month];
       }
     }
   }

   SXW->temp = swGenOut->temp;
   SXW->ppt = swGenOut->ppt;
   SXW->aet = swGenOut->aet;
}

void _sxw_sw_clear_transp(void) {
/*======================================================*/
	int k;

	Mem_Set(SXW->transpTotal, 0, SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	ForEachVegType(k) {
		Mem_Set(SXW->transpVeg[k], 0, SXW->NPds * SXW->NSoLyrs * sizeof(RealD));
	}
}

static void _update_transp_coeff(void) {
    /*======================================================*/
    /* copy the relative root distribution to soilwat's layers */
    /* POTENTIAL BUG:  if _roots_current is all zero but for some
     *   reason the productivity values (esp. %live) are >0,
     *   there can be transpiration but nowhere for it to come
     *   from in the soilwat model.
     */
    GrpIndex g;
    LyrIndex l;
    RealF sum[NVEGTYPES] = {0.};

    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_TREES) {
        SoilWatAll.Site.transp_coeff[SW_TREES][l] = 0.;

        ForEachGroup(g) {
          if (RGroup[g]->veg_prod_type == SW_TREES)
              if (getNTranspLayers(SW_TREES))
                  SoilWatAll.Site.transp_coeff[SW_TREES][l] += 
                      (RealF) SXWResources->_roots_max[Ilg(l, g)] * 
                            RGroup[g]->rgroupFractionOfVegTypeBiomass;
        }
        sum[SW_TREES] += SoilWatAll.Site.transp_coeff[SW_TREES][l];
    }

    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_SHRUB) {
        SoilWatAll.Site.transp_coeff[SW_SHRUB][l] = 0.;

        ForEachGroup(g)
        if (RGroup[g]->veg_prod_type == SW_SHRUB) {
            if (getNTranspLayers(SW_SHRUB))
                SoilWatAll.Site.transp_coeff[SW_SHRUB][l] += 
                    (RealF) SXWResources->_roots_max[Ilg(l, g)] *
                        RGroup[g]->rgroupFractionOfVegTypeBiomass;

            /*printf("* lyr=%d, group=%s(%d), type=%d, tl=%d, rootmax=%f, relsize2=%f, trco=%f\n",
              l, RGroup[g]->name, g, RGroup[g]->veg_prod_type, getNTranspLayers(RGroup[g]->veg_prod_type),
              SXWResources->_roots_max[Ilg(l, g)], getRGroupRelsize(g), y->transp_coeff[SW_SHRUB]);
             */
        }
        sum[SW_SHRUB] += SoilWatAll.Site.transp_coeff[SW_SHRUB][l];
    }

    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_GRASS) {
        SoilWatAll.Site.transp_coeff[SW_GRASS][l] = 0.;

        ForEachGroup(g) {
          if (RGroup[g]->veg_prod_type == SW_GRASS)
              if (getNTranspLayers(SW_GRASS))
                  SoilWatAll.Site.transp_coeff[SW_GRASS][l] +=
                       (RealF) SXWResources->_roots_max[Ilg(l, g)] * 
                            RGroup[g]->rgroupFractionOfVegTypeBiomass;
         }
         sum[SW_GRASS] += SoilWatAll.Site.transp_coeff[SW_GRASS][l];
    }

    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_FORBS) {
        SoilWatAll.Site.transp_coeff[SW_FORBS][l] = 0.;

        ForEachGroup(g) {
          if (RGroup[g]->veg_prod_type == SW_FORBS)
              if (getNTranspLayers(SW_FORBS))
                  SoilWatAll.Site.transp_coeff[SW_FORBS][l] +=
                     (RealF) SXWResources->_roots_max[Ilg(l, g)] *
                        RGroup[g]->rgroupFractionOfVegTypeBiomass;
        }
        sum[SW_FORBS] += SoilWatAll.Site.transp_coeff[SW_FORBS][l];
    }

    /* normalize coefficients to 1.0 If sum is 0, then the transp_coeff is also 0. */
    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_TREES) {
      if (!ZRO(sum[SW_TREES])) SoilWatAll.Site.transp_coeff[SW_TREES][l] /= sum[SW_TREES];
    }
    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_SHRUB) {
      if (!ZRO(sum[SW_SHRUB])) SoilWatAll.Site.transp_coeff[SW_SHRUB][l] /= sum[SW_SHRUB];
    }
    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_GRASS) {
      if (!ZRO(sum[SW_GRASS])) SoilWatAll.Site.transp_coeff[SW_GRASS][l] /= sum[SW_GRASS];
    }
    ForEachTranspLayer(l, SoilWatAll.Site.n_transp_lyrs, SW_FORBS) {
      if (!ZRO(sum[SW_FORBS])) SoilWatAll.Site.transp_coeff[SW_FORBS][l] /= sum[SW_FORBS];
    }

    /*printf("'_update_transp_coeff': ShrubTranspCoef: ");
    ForEachSoilLayer(l) {
      printf("[%d] = %.4f - ", l, SW_Site.lyr[l]->transp_coeff[SW_SHRUB]);
    }
    printf("\n");
     */

}

static void _update_productivity(RealF sizes[]) {
    
    GrpIndex g;
    TimeInt m;
    IntUS k;

    SW_VEGPROD *v = &SoilWatAll.VegProd;
    RealF totbmass = 0.0,
            *bmassg,
    vegTypeBiomass[NVEGTYPES] = {0.};

    bmassg = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups,
                                 sizeof(RealF), "_update_productivity", &LogInfo);


    // totbmass: total biomass in g/m2
    // vegTypeBiomass: biomass for each of the SOILWAT2 vegetation types in g/m2
    ForEachGroup(g) {
        bmassg[g] = sizes[g] / Globals->plotsize; // gram per plot -> gram per m2
        totbmass += bmassg[g];

        // Sum biomass of resource groups per SOILWAT2 vegetation type
        vegTypeBiomass[RGroup[g]->veg_prod_type] += bmassg[g];
    }

    // Calculate the cover fraction of each SOILWAT2 vegetation type (sum = 1):
    // here, approximate cover with the contribution of vegetation type biomass
    // to total biomass
    if (GT(totbmass, 0.)) {
        ForEachVegType(k) {
          v->veg[k].cov.fCover = vegTypeBiomass[k] / totbmass;
        }

        //TODO: figure how to calculate bareground fraction.
        v->bare_cov.fCover = 0;
    } else {

        ForEachVegType(k) {
            v->veg[k].cov.fCover = 0.0;
        }
        v->bare_cov.fCover = 1;
    }

    // Calculate the biomass contribution of a STEPWAT2 resource group relative
    // to its SOILWAT2 vegetation type
    ForEachGroup(g) {
        // Get biomass of SOILWAT2 vegetation type for current resource group
        if (GT(vegTypeBiomass[RGroup[g]->veg_prod_type], 0.))
            RGroup[g]->rgroupFractionOfVegTypeBiomass = bmassg[g] / vegTypeBiomass[RGroup[g]->veg_prod_type];
        else
            RGroup[g]->rgroupFractionOfVegTypeBiomass = 0;
    }

    // Calculate monthly biomass, litter, and pct live:
    // Note: SOILWAT2 expects biomass per vegetation type as if that type
    //       covered 100% of a square meter --> scale biomass by 1 / fCover
    ForEachMonth(m) {

        ForEachVegType(k) {
            v->veg[k].pct_live[m] = 0.;
            v->veg[k].biomass[m] = 0.;
            v->veg[k].litter[m] = 0.;
        }

        if (GT(totbmass, 0.)) {
            ForEachGroup(g) {
              k = RGroup[g]->veg_prod_type;
              v->veg[k].pct_live[m] += SXWResources->_prod_pctlive[Igp(g, m)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;

              v->veg[k].biomass[m] += SXWResources->_prod_bmass[Igp(g, m)] * 
                                      bmassg[g] / v->veg[k].cov.fCover;

              v->veg[k].litter[m] += vegTypeBiomass[k] * SXWResources->_prod_litter[g][m];
            }
        }
    }

    Mem_Free(bmassg);
}
