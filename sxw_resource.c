/**
 * \file sxw_resource.c
 * \brief Translates transpiration values from [SOILWAT](\ref sw_src) into 
 *        [Steppe](\ref STEPPE) resource values.
 * 
 * \author CWB (initial programming)
 * \date 21 May 2002
 * \ingroup SXW_PRIVATE
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include "sw_src/rands.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sxw_module.h"
#include "sxw_vars.h"
#include "sw_src/pcg/pcg_basic.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
/* for steppe, see ST_globals.h */

//extern SW_SITE SW_Site;
//extern SW_SOILWAT SW_Soilwat;
//extern SW_VEGPROD SW_VegProd;

extern SXW_resourceType* SXWResources;

/* ------ Running Averages ------ */
extern
  transp_t* transp_window;

extern
  pcg32_random_t resource_rng;

//void _print_debuginfo(void);

/*************** Local Function Declarations ***************/
/***********************************************************/

static void _transp_contribution_by_group(RealF use_by_group[]);

/***********************************************************/
/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_root_phen(void) {
/*======================================================*/
/* should only be called once, after root distr. and
 * phenology tables are read
 */

	LyrIndex y;
	GrpIndex g;
	TimeInt p;

	for (y = 0; y < (Globals->grpCount * SXW->NPds * SXW->NTrLyrs); y++)
		SXWResources->_rootsXphen[y] = 0;

	ForEachGroup(g)
	{
		int nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		for (y = 0; y < nLyrs; y++) {
			ForEachTrPeriod(p) {
				SXWResources->_rootsXphen[Iglp(g, y, p)] = SXWResources->_roots_max[Ilg(y, g)] * SXWResources->_phen[Igp(g, p)];
			}
		}
	}
}

void _sxw_update_resource(void) {
/*======================================================*/
/* Determines resources available to each STEPPE functional group each year.
 * The first step is to get the current biomass for each STEPPE functional group.
 * Second, we re-calculate the relative active roots in each layer in each month
 * using the updated functional group biomass. Third, we divide transpiration to
 * each STEPPE functional group based on the matching of active roots in each
 * soil layer and month with transpiration in each layer and month. Finally, we
 * scale resources available (cm) to resources in terms of grams of biomass */

  RealF *sizes;
  GrpIndex g;

  sizes = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF), "_sxw_update_resource");

	ForEachGroup(g)
	{
		sizes[g] = 0.;
//printf("_sxw_update_resource()RGroup Name= %s, RGroup[g]->regen_ok=%d \n ", RGroup[g]->name, RGroup[g]->regen_ok);
		if (!RGroup[g]->regen_ok)
			continue;
		sizes[g] = RGroup_GetBiomass(g);
	}

    /* Update the active relative roots based on current biomass values */
	_sxw_update_root_tables(sizes);

	/* Assign transpiration (resource availability) to each STEPPE functional group */
	_transp_contribution_by_group(SXWResources->_resource_cur);

        /* Scale transpiration resources by a constant, bvt, to convert resources
         * (cm) to biomass that can be supported by those resources (g/cm) */
	ForEachGroup(g)
	{
//printf("%s: _resource_cur pre-multiplication = %f. Post multiplication = %f.\n",RGroup[g]->name, 
                                //_resource_cur[g], SXWResources->_resource_cur[g] * RGroup[g]->_bvt);
		SXWResources->_resource_cur[g] = SXWResources->_resource_cur[g] * RGroup[g]->_bvt;
	}

    Mem_Free(sizes);
}

void _sxw_update_root_tables( RealF sizes[] ) {
/*======================================================*/
/* Updates the active relative roots array based on sizes, which contains the groups'
 * actual biomass in grams. This array in utilized in partitioning of transpiration
 * (resources) to each STEPPE functional group. */

	GrpIndex g;
	LyrIndex l, nLyrs;
	TimeInt p;
	RealD x, sum_all_veg_types;
	int t;

	/* Set some things to zero where 4 refers to Tree, Shrub, Grass, Forb */
	Mem_Set(SXWResources->_roots_active_sum, 0, NVEGTYPES * SXW->NPds * SXW->NTrLyrs * sizeof(RealD));

        /* Calculate the active roots in each month and soil layer for each STEPPE
         * functional group based on the functional group biomass this year */
	ForEachGroup(g)
	{
		nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		for (l = 0; l < nLyrs; l++) {
			ForEachTrPeriod(p)
			{
				x = SXWResources->_rootsXphen[Iglp(g, l, p)] * sizes[g];
				SXWResources->_roots_active[Iglp(g, l, p)] = x;
				SXWResources->_roots_active_sum[Itlp(RGroup[g]->veg_prod_type, l, p)] += x;
			}
		}
	}

	/* Rescale SXWResources->_roots_active_sum to represent the relative "activity" of a 
         * STEPPE group's roots in a given layer in a given month */
  /* Details: for each soil layer `l` and each month (trperiod) `p`, the sum
     of `_roots_active_rel[Iglp(g, l, p)]` across rgroups `g` must be 1;
     see its use by function _transp_contribution_by_group() */
	ForEachGroup(g)
	{
		nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		
		for (l = 0; l < nLyrs; l++) {
			ForEachTrPeriod(p)
			{
        sum_all_veg_types = 0;
        ForEachVegType(t) {
          sum_all_veg_types += SXWResources->_roots_active_sum[Itlp(t, l, p)];
        }

				SXWResources->_roots_active_rel[Iglp(g, l, p)] = ZRO(sum_all_veg_types) ? 0. 
                                         : SXWResources->_roots_active[Iglp(g, l, p)] / sum_all_veg_types;
			}
		}
	}

}

static void _transp_contribution_by_group(RealF use_by_group[]) {
    /*======================================================*/
    /* use_by_group is the amount of transpiration (cm) assigned to each STEPPE
     * functional group. Must call _update_root_tables() before this.
     * Compute each group's amount of transpiration from SOILWAT2
     * based on its biomass, root distribution, and phenological
     * activity. This represents "normal" resources or transpiration each year.
     * In cases where transpiration is significantly below the mean due to low
     * biomass (e.g. fire years), additional transpiration is added. */

    GrpIndex g;
    TimeInt p;
    int nLyrs, t;
    LyrIndex l;
    RealD *transp, proportion_total_resources, average_proportion;
    RealF sumUsedByGroup = 0., sumTranspTotal = 0., TranspRemaining = 0.;
    RealF transp_ratio, transp_ratio_sd;
    transp_window->added_transp = 0;

    // year 0 is a set up year. No need to calculate transpiration.
    // if there are multiple iterations the last year will run twice;
    // once for data and once for tear down. The second "last year" is
    // equivalent to year 0.
    if(Globals->currYear == 0 || Globals->currYear == transp_window->lastYear) {
      transp_window->average = 0;
      transp_window->ratio_average = 0;
      transp_window->sum_of_sqrs = 0;
      transp_window->oldest_index = 0;
      return; //no point calculating anything since SOILWAT hasn't run
    }

    ForEachGroup(g) //Steppe functional group
    {
        use_by_group[g] = 0.;
        transp = SXW->transpTotal;

        //Loops through each month and calculates amount of transpiration for each STEPPE functional group
        //according to whether that group has active living roots in each soil layer for each month
        /* Details: for each soil layer `l` and each month (trperiod) `p`, the
           sum of `_roots_active_rel[Iglp(g, l, p)]` across rgroups `g` must
           be 1; only if they sum to 1, can they can be used as proper weights
           for `transp[Ilp(l, p)]` to be split up among rgroups */
        ForEachTrPeriod(p) {
            nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
            for (l = 0; l < nLyrs; l++) {
                use_by_group[g] += (RealF) (SXWResources->_roots_active_rel[Iglp(g, l, p)] * transp[Ilp(l, p)]);
            }
        }

        sumUsedByGroup += use_by_group[g];
        //printf(" sumUsedByGroup in transp=%f \n",sumUsedByGroup);
    }

    // Rescale transpiration based on space (min_res_req)
    ForEachGroup(g){
      if(use_by_group[g] != 0){ //Prevents NaN errors
        //Proportion of total traspiration that this group would recieve without accounting for space
        proportion_total_resources = use_by_group[g] / sumUsedByGroup;
        //Average of the two proportions. This equally weights proportion_total_resources and min_res_req.
        average_proportion = (proportion_total_resources + RGroup[g]->min_res_req) / 2;
        //Recalculate this group's resources.
        use_by_group[g] = sumUsedByGroup * average_proportion;
      }
    }

    //sumUsedByGroup needs to be synced because of potention floating point errors
    sumUsedByGroup = 0;
    ForEachGroup(g){
      sumUsedByGroup += use_by_group[g];
    }

    //Very small amounts of transpiration remain and not perfectly partitioned to functional groups.
    //This check makes sure any remaining transpiration is divided proportionately among groups.
    ForEachTrPeriod(p) {
        for (t = 0; t < SXW->NSoLyrs; t++)
            sumTranspTotal += SXW->transpTotal[Ilp(t, p)];
    }

    TranspRemaining = sumTranspTotal - sumUsedByGroup;

    /* ------------- Begin testing to see if additional transpiration is necessary ------------- */

    transp_ratio = sumTranspTotal / SXW->ppt;

    // Determines if the current year transpiration/ppt is greater than 1 standard deviations away
    // from the mean. If TRUE, add additional transpiration.
    if(transp_window->size >= Globals->currYear) //we need to do a running average
    {
        // add transpiration to the window
        transp_window->transp[transp_window->oldest_index] = sumTranspTotal;
        //update the average
        transp_window->average = get_running_mean(Globals->currYear,transp_window->average, sumTranspTotal);
        //add the ratio value to the window
        transp_window->ratios[transp_window->oldest_index] = transp_ratio;
        //save the last mean. we will need it to calculate the sum of squares
        RealF last_ratio = transp_window->ratio_average;
        //calculate the running mean
        transp_window->ratio_average = get_running_mean(Globals->currYear,transp_window->ratio_average,transp_ratio);
        //calculate the running sum of squares
        RealF ssqr = get_running_sqr(last_ratio, transp_window->ratio_average, transp_ratio);
        //add the calculated sum of squares to the running total
        transp_window->sum_of_sqrs += ssqr;
        //add the calculated sum of squares to the array
        transp_window->SoS_array[transp_window->oldest_index] = ssqr;
        //calculate the standard deviation
        transp_ratio_sd = final_running_sd(Globals->currYear, transp_window->sum_of_sqrs);

    } else { //we need to do a moving window
        //add the new value, subtract the old value from the average;
        transp_window->average += (sumTranspTotal - transp_window->transp[transp_window->oldest_index])/transp_window->size;
        //add the new value, subtract the old value from the ratio average;
        transp_window->ratio_average += (transp_ratio - transp_window->ratios[transp_window->oldest_index])/transp_window->size;
        //put the new transpiration in the window
        transp_window->transp[transp_window->oldest_index] = sumTranspTotal;
        //put the new ratio in the window
        transp_window->ratios[transp_window->oldest_index] = transp_ratio;
        // calculate the new sum of squares value
        RealF ssqr = (transp_ratio - transp_window->ratio_average) * (transp_ratio - transp_window->ratio_average);
        // add the new sum of squares, subtract the old.
        transp_window->sum_of_sqrs += ssqr - transp_window->SoS_array[transp_window->oldest_index];
        // update the sum of squares window.
        transp_window->SoS_array[transp_window->oldest_index] =  ssqr;
        //calculate the standard deviation
        transp_ratio_sd = final_running_sd(transp_window->size, transp_window->sum_of_sqrs);
    }

    //printf("Year %d: ratio = %f, mean = %f, sos = %f sd = %f\n",Globals->currYear,
              // transp_ratio,transp_window->ratio_average, transp_window->sum_of_sqrs, transp_ratio_sd);

    // If this year's transpiration is notably low (1 sd below the mean), add additional transpired water
    if (transp_ratio < (transp_window->ratio_average - 1 * transp_ratio_sd)) {
            RealF min = transp_window->ratio_average - transp_ratio_sd;
            RealF max = transp_window->ratio_average + transp_ratio_sd;

            // This transpiration will be added
            transp_window->added_transp = (1 - transp_ratio / RandUniFloatRange(min, max, &resource_rng)) * transp_window->average;
            if(transp_window->added_transp < 0){
                LogError(logfp, LOGNOTE, "sxw_resource: Added transpiration less than 0.\n");
            }
            //printf("Year %d:\tTranspiration to add: %f\n",Globals->currYear,add_transp);

            /* -------------------- Recalculate the window including added_transp in the current year -------------------- */
            RealF added_transp_ratio = transp_window->added_transp / SXW->ppt;
            //add added_transp to the average. This is technically part of the current year, so no need to subtract anything.
            transp_window->average += transp_window->added_transp/transp_window->size;
            //add added_transp ratio to the ratio average. This is technically part of the current year, so no need to subtract anything.
            transp_window->ratio_average += (added_transp_ratio)/transp_window->size;
            //put the new transpiration in the window. Note: oldest_index has not been incremented, so it points to what was just added
            transp_window->transp[transp_window->oldest_index] += transp_window->added_transp;
            //put the new ratio in the window. Note: oldest_index has not been incremented, so it points to what was just added
            transp_window->ratios[transp_window->oldest_index] += added_transp_ratio;
            // calculate the new sum of squares value
            RealF totalTranspRatio = (sumTranspTotal + transp_window->added_transp)/SXW->ppt;
            RealF ssqr = (totalTranspRatio - transp_window->ratio_average) * (totalTranspRatio - transp_window->ratio_average);
            // Subtract the sum of squares calculated above, which was stored in the array. Replace it with what was just calculated.
            transp_window->sum_of_sqrs += ssqr - transp_window->SoS_array[transp_window->oldest_index];
            // replace the sum of squares with what we just calculated
            transp_window->SoS_array[transp_window->oldest_index] =  ssqr;

            /* Adds the additional transpiration to the remaining transpiration
             * so it can be distributed proportionally to the functional groups. */
            TranspRemaining += transp_window->added_transp;
    }

    //oldest_index++ accounting for the array bounds
    transp_window->oldest_index = (transp_window->oldest_index + 1) % transp_window->size;

    /* ------------ End testing to see if additional transpiration is necessary ---------- */

    // If there is transpiration this year add the remaining (and added) transpiration to each functional group.
    // Else the sum of transpiration is 0 and for each group there is also zero transpiration.
    if (!ZRO(sumUsedByGroup)) {

        ForEachGroup(g) {
            use_by_group[g] += (use_by_group[g] / sumUsedByGroup) * TranspRemaining;
            //printf("for groupName= %s, after sum use_by_group[g]= %f \n",RGroup[g]->name,use_by_group[g]);
        }
    }

    //remember the last year. When setting up for a new iteration the same year will appear twice, and we want to skip it the second time
    transp_window->lastYear = Globals->currYear; 
}
