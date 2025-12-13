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
 *          to use the different VegProds now used in soilwat
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

#include <stdlib.h>

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
	SW_VEGPROD_SIM *v = &SoilWatRun.VegProdSim;

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
  SW_WEATHER_INPUTS *w = &SoilWatRun.WeatherIn;
  SW_WEATHER_HIST **wh = &SoilWatRun.RunIn.weathRunAllHist;
  SW_SKY_INPUTS *sky = &SoilWatRun.RunIn.SkyRunIn;

  deallocateAllWeather(&SoilWatRun.RunIn.weathRunAllHist);
  w->n_years = 1;
  w->startYear = SoilWatRun.ModelIn.startyr + Globals->currYear - 1;
  SW_WTH_allocateAllWeather(wh, w->n_years, &LogInfoSW);

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
    *wh,
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
    w->fixWeatherData,
    sky->cloudcov,
    sky->windspeed,
    sky->r_humidity,
    SoilWatRun.RunIn.ModelRunIn.elevation,
    SoilWatRun.ModelSim.cum_monthdays,
    SoilWatRun.ModelSim.days_in_month,
    &LogInfoSW
    );

  finalizeAllWeather(&SoilWatRun.MarkovIn, w, *wh, SoilWatRun.ModelSim.cum_monthdays,
                      SoilWatRun.ModelSim.days_in_month, NULL, FALSE, &LogInfoSW); // run the weather
}


void _sxw_sw_run(void) {
/*======================================================*/
   SW_OUT_RUN *OutRun = &SoilWatRun.OutRun;
   LyrIndex lyrno;
   TimeInt month;
   int vegType;

 	SoilWatRun.ModelSim.year = SoilWatRun.ModelIn.startyr + Globals->currYear-1;

   // Copy global values to SOILWAT2's SW_ALL to work with correct values
   SoilWatRun.ModelSim.runModelIterations = SuperGlobals.runModelIterations;
   SoilWatRun.ModelIn.runModelYears = SuperGlobals.runModelYears;
 	SW_CTL_run_current_year(&SoilWatRun, &SoilWatDomain.OutDom, &LogInfoSW);


   // Copy the results of SOILWAT2's SXW values
   for(lyrno = 0; lyrno < SXW->NSoLyrs; lyrno++) {
     ForEachTrPeriod(month) {
       SXW->swc[Ilp(lyrno, month)] = OutRun->swc[lyrno][month];
       SXW->transpTotal[Ilp(lyrno, month)] = OutRun->transpTotal[lyrno][month];

       SXW->ppt_monthly[month] = OutRun->ppt_monthly[month];
       SXW->temp_monthly[month] = OutRun->temp_monthly[month];

       ForEachVegType(vegType) {
         SXW->transpVeg[vegType][Ilp(lyrno, month)] =
                                   OutRun->transpVeg[vegType][lyrno][month];
       }
     }
   }

   SXW->temp = OutRun->temp;
   SXW->ppt = OutRun->ppt;
   SXW->aet = OutRun->aet;
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
    unsigned int k;
    RealF sum[NVEGTYPES] = {0.};

    ForEachVegType(k) {
        ForEachTranspLayer(l, SoilWatRun.SiteSim.n_transp_lyrs, k) {
            SoilWatRun.RunIn.SoilRunIn.transp_coeff[k][l] = 0.;

            ForEachGroup(g) {
                if (RGroup[g]->veg_prod_type == k) {
                    if (getNTranspLayers(k)) {
                        SoilWatRun.RunIn.SoilRunIn.transp_coeff[k][l] +=
                            SXWResources->_roots_max[Ilg(l, g)] *
                                RGroup[g]->rgroupFractionOfVegTypeBiomass;
                    }
                }
            }

            sum[k] += SoilWatRun.RunIn.SoilRunIn.transp_coeff[k][l];
        }
    }

    /* normalize coefficients to 1.0 If sum is 0, then the transp_coeff is also 0. */
    ForEachVegType(k) {
        ForEachTranspLayer(l, SoilWatRun.SiteSim.n_transp_lyrs, k) {
            if (!ZRO(sum[k])) {
                SoilWatRun.RunIn.SoilRunIn.transp_coeff[k][l] /= sum[k];
            }
        }
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

    SW_VEGPROD_RUN_INPUTS *v = &SoilWatRun.RunIn.VegProdRunIn;
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

              if (GT(v->veg[k].cov.fCover, 0.)) {
                v->veg[k].pct_live[m] += SXWResources->_prod_pctlive[Igp(g, m)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;

                v->veg[k].biomass[m] += SXWResources->_prod_bmass[Igp(g, m)] *
                                        bmassg[g] / v->veg[k].cov.fCover;

                v->veg[k].litter[m] += vegTypeBiomass[k] * SXWResources->_prod_litter[g][m];
              }
            }
        }
    }

    free(bmassg);
}
