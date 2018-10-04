/********************************************************/
/********************************************************/
/*  Source file: sxw_resource.c
 *  Type: module
 *  Purpose: Compute resource vector for STEPPE based on
 *           transpiration values from SOILWAT.
 *  Dependency:  sxw.c
 *  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model. */
/*  History:
 *     (21-May-2002) -- INITIAL CODING - cwb
 *     19-Jun-2003 - cwb - Added RealD (double precision)
 *                 types for internal dynamic matrices and
 *                 other affected variables.  See notes in
 *                 sxw.c. */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include "generic.h"
#include "rands.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "SW_Defines.h"
#include "sxw.h"
#include "sxw_module.h"
#include "sxw_vars.h"
#include "SW_Control.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Files.h"
#include "SW_Times.h"
#include "sw_src/pcg/pcg_basic.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
/* for steppe, see ST_globals.h */

//extern SW_SITE SW_Site;
//extern SW_SOILWAT SW_Soilwat;
//extern SW_VEGPROD SW_VegProd;


/*************** Local Variable Declarations ***************/
/***********************************************************/
/* malloc'ed and maybe read in sxw.c but used here */
/* ----- 3d arrays ------- */
extern
  RealD * _rootsXphen, /* relative roots X phen by layer & group */
        * _roots_active, /*relative to the total roots_phen_lyr_group */
        * _roots_active_rel;


/* ----- 2D arrays ------- */

extern
       /* rgroup by layer */
  RealD * _roots_max,     /* root distribution with depth for STEPPE functional groups, read from input */
        * _roots_active_sum, /* active roots in each month and soil layer for STEPPE functional groups in the current year */

       /* rgroup by period */
        * _phen;          /* phenologic activity for each month for STEPPE functional groups, read from input */

extern
  RealF _resource_pr[MAX_RGROUPS],  /* resource convertable to pr */
        _resource_cur[MAX_RGROUPS]; /* resources currently availible by group*/
        RealF added_transp; /* transpiration added for the current year */

extern
  RealF _bvt;

/* ------ Running Averages ------ */
extern
  transp_t transp_window;

extern 
  pcg32_random_t resource_rng;

// _transp_contribution_by_group() has different behavior for production years and setup years.
// this variable keeps track of what year last year was, to make sure that the year actually
// incremented. If it did NOT, then this is a setup year.
int lastYear;

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

	for (y = 0; y < (Globals.grpCount * SXW.NPds * SXW.NTrLyrs); y++)
		_rootsXphen[y] = 0;

	ForEachGroup(g)
	{
		int nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		for (y = 0; y < nLyrs; y++) {
			ForEachTrPeriod(p) {
				_rootsXphen[Iglp(g, y, p)] = _roots_max[Ilg(y, g)] * _phen[Igp(g, p)];
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
  
  RealF sizes[MAX_RGROUPS] = {0.};
  GrpIndex g;

  #ifdef SXW_BYMAXSIZE
    int i;
    SppIndex sp;
    ForEachGroup(g) {
      sizes[g] = 0.;
      if (RGroup[g]->regen_ok) {
        ForEachGroupSpp(sp, g, i) {
          sizes[g] += Species[sp]->mature_biomass;
        }
      }
    }
  #else
	ForEachGroup(g)
	{
		//RGroup[g]->veg_prod_type
		sizes[g] = 0.;
//printf("_sxw_update_resource()RGroup Name= %s, RGroup[g]->regen_ok=%d \n ", RGroup[g]->name, RGroup[g]->regen_ok);
		if (!RGroup[g]->regen_ok)
			continue;
		sizes[g] = RGroup_GetBiomass(g);
	}
  #endif
        
        /* Update the active relative roots based on current biomass values */
	_sxw_update_root_tables(sizes);
        
	/* Assign transpiration (resource availability) to each STEPPE functional group */
	_transp_contribution_by_group(_resource_cur);
        
        /* Scale transpiration resources by a constant, bvt, to convert resources 
         * (cm) to biomass that can be supported by those resources (g/cm) */
	ForEachGroup(g)
	{
//printf("for groupName= %smresource_cur prior to multiplication: %f\n",RGroup[g]->name, _resource_cur[g]);
		_resource_cur[g] = _resource_cur[g] * _bvt;
//printf("for groupName= %s, resource_cur post multiplication: %f\n\n",Rgroup[g]->name, _resource_cur[g]);
	}
/* _print_debuginfo(); */
}

void _sxw_update_root_tables( RealF sizes[] ) {
/*======================================================*/
/* Updates the active relative roots array based on sizes, which contains the groups'
 * actual biomass in grams. This array in utilized in partitioning of transpiration
 * (resources) to each STEPPE functional group. */

	GrpIndex g;
	LyrIndex l;
	TimeInt p;
	RealD x;
	int t,nLyrs;

	/* Set some things to zero where 4 refers to Tree, Shrub, Grass, Forb */
	Mem_Set(_roots_active_sum, 0, 4 * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));
        
        /* Calculate the active roots in each month and soil layer for each STEPPE
         * functional group based on the functional group biomass this year */
	ForEachGroup(g)
	{
		t = RGroup[g]->veg_prod_type-1;
		nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		for (l = 0; l < nLyrs; l++) {
			ForEachTrPeriod(p)
			{
				x = _rootsXphen[Iglp(g, l, p)] * sizes[g];
				_roots_active[Iglp(g, l, p)] = x;
				_roots_active_sum[Itlp(t, l, p)] += x;
			}
		}
	}

	/* Rescale _roots_active_sum to represent the relative "activity" of a 
         * STEPPE group's roots in a given layer in a given month */
	ForEachGroup(g)
	{
		int t = RGroup[g]->veg_prod_type-1;
		int nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		for (l = 0; l < nLyrs; l++) {
			ForEachTrPeriod(p)
			{
				_roots_active_rel[Iglp(g, l, p)] =
				ZRO(_roots_active_sum[Itlp(t,l,p)]) ?
						0. :
						_roots_active[Iglp(g, l, p)]
								/ _roots_active_sum[Itlp(t,l, p)];
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
    LyrIndex l;
    int t;
    RealD *transp;
    RealF sumUsedByGroup = 0., sumTranspTotal = 0., TranspRemaining = 0.;
    RealF transp_ratio;
    added_transp = 0;
    RealF transp_ratio_sd;

    // year 0 is a set up year. No need to calculate transpiration.
    // if there are multiple iterations the last year will run twice;
    // once for data and once for tear down. The second "last year" is
    // equivalent to year 0.
    if(Globals.currYear == 0 || Globals.currYear == lastYear) {
      transp_window.average = 0;
      transp_window.ratio_average = 0;
      transp_window.sum_of_sqrs = 0;
      if(Globals.transp_window > MAX_WINDOW){
                LogError(logfp, LOGNOTE, "sxw_resource: Transp_window specified in inputs is greater than maximum window.\nInput: %d\nMaximum: %d\nSetting window to %d.\n",Globals.transp_window, MAX_WINDOW,MAX_WINDOW);
                Globals.transp_window = MAX_WINDOW;
      }
      transp_window.size = Globals.transp_window;
      transp_window.oldest_index = 0;
      return; //no point calculating anything since SOILWAT hasn't run
    }

    ForEachGroup(g) //Steppe functional group
    {
        use_by_group[g] = 0.;
        t = RGroup[g]->veg_prod_type - 1;

        switch (t) {
            case 0://Tree
                transp = SXW.transpVeg[SW_TREES];
                break;
            case 1://Shrub
                transp = SXW.transpVeg[SW_SHRUB];
                break;
            case 2://Grass
                transp = SXW.transpVeg[SW_GRASS];
                break;
            case 3://Forb
                transp = SXW.transpVeg[SW_FORBS];
                break;
            default:
                transp = SXW.transpTotal;
                break;
        }

        //Loops through each month and calculates amount of transpiration for each STEPPE functional group
        //according to whether that group has active living roots in each soil layer for each month

        ForEachTrPeriod(p) {
            int nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
            for (l = 0; l < nLyrs; l++) {
                use_by_group[g] += (RealF) (_roots_active_rel[Iglp(g, l, p)] * transp[Ilp(l, p)]);
            }
        }
        //printf("for groupName= %s, use_by_group[g] in transp= %f \n",RGroup[g]->name,use_by_group[g] );

        sumUsedByGroup += use_by_group[g];
        //printf(" sumUsedByGroup in transp=%f \n",sumUsedByGroup);
    }

    //Very small amounts of transpiration remain and not perfectly partitioned to functional groups.
    //This check makes sure any remaining transpiration is divided proportionately among groups.

    ForEachTrPeriod(p) {
        for (t = 0; t < SXW.NSoLyrs; t++)
            sumTranspTotal += SXW.transpTotal[Ilp(t, p)];
    }

    TranspRemaining = sumTranspTotal - sumUsedByGroup;
    //printf(" sumTranspTotal=%f, sumUsedByGroup=%f  TranspRemaining=%f"\n", sumTranspTotal, sumUsedByGroup, TranspRemaining);

    /* ------------- Begin testing to see if additional transpiration is necessary ------------- */

    transp_ratio = sumTranspTotal / SXW.ppt;

    // Determines if the current year transpiration/ppt is greater than 1 standard deviations away
    // from the mean. If TRUE, add additional transpiration.
    if(transp_window.size >= Globals.currYear) //we need to do a running average
    {
        // add transpiration to the window
        transp_window.transp[transp_window.oldest_index] = sumTranspTotal;
        //update the average
        transp_window.average = get_running_mean(Globals.currYear,transp_window.average, sumTranspTotal);
        //add the ratio value to the window
        transp_window.ratios[transp_window.oldest_index] = transp_ratio;
        //save the last mean. we will need it to calculate the sum of squares
        RealF last_ratio = transp_window.ratio_average;
        //calculate the running mean
        transp_window.ratio_average = get_running_mean(Globals.currYear,transp_window.ratio_average,transp_ratio);
        //calculate the running sum of squares
        RealF ssqr = get_running_sqr(last_ratio, transp_window.ratio_average, transp_ratio);
        //add the calculated sum of squares to the running total
        transp_window.sum_of_sqrs += ssqr;
        //add the calculated sum of squares to the array
        transp_window.SoS_array[transp_window.oldest_index] = ssqr;
        //calculate the standard deviation
        transp_ratio_sd = final_running_sd(Globals.currYear, transp_window.sum_of_sqrs);

    } else { //we need to do a moving window
        //add the new value, subtract the old value from the average;
        transp_window.average += (sumTranspTotal - transp_window.transp[transp_window.oldest_index])/transp_window.size;
        //add the new value, subtract the old value from the ratio average;
        transp_window.ratio_average += (transp_ratio - transp_window.ratios[transp_window.oldest_index])/transp_window.size;
        //put the new transpiration in the window
        transp_window.transp[transp_window.oldest_index] = sumTranspTotal;
        //put the new ratio in the window
        transp_window.ratios[transp_window.oldest_index] = transp_ratio;
        // calculate the new sum of squares value
        RealF ssqr = (transp_ratio - transp_window.ratio_average) * (transp_ratio - transp_window.ratio_average);
        // add the new sum of squares, subtract the old.
        transp_window.sum_of_sqrs += ssqr - transp_window.SoS_array[transp_window.oldest_index];
        // update the sum of squares window.
        transp_window.SoS_array[transp_window.oldest_index] =  ssqr;
        //calculate the standard deviation
        transp_ratio_sd = final_running_sd(transp_window.size, transp_window.sum_of_sqrs);
    }

    //printf("Year %d: ratio = %f, mean = %f, sos = %f sd = %f\n",Globals.currYear,
              // transp_ratio,transp_window.ratio_average, transp_window.sum_of_sqrs, transp_ratio_sd);

    // If this year's transpiration is notably low (1 sd below the mean), add additional transpired water
    if (transp_ratio < (transp_window.ratio_average - 1 * transp_ratio_sd)) {
        //printf("Year %d below 1 sd: ratio = %f, average = %f, sd = %f\n", Globals.currYear,transp_ratio,
                                                          //transp_window.ratio_average, transp_ratio_sd);
            RealF min = transp_window.ratio_average - transp_ratio_sd;
            RealF max = transp_window.ratio_average + transp_ratio_sd;

            // This transpiration will be added 
            added_transp = (1 - transp_ratio / RandUniFloatRange(min, max, &resource_rng)) * transp_window.average;
            if(added_transp < 0){
                LogError(logfp, LOGNOTE, "sxw_resource: Added transpiration less than 0.\n");
            }
            //printf("Year %d:\tTranspiration to add: %f\n",Globals.currYear,add_transp);
            //printf("TranspRemaining: %f\tTranspRemaining+add_transp: %f\n",TranspRemaining,add_transp+TranspRemaining);

            /* -------------------- Recalculate the window including added_transp in the current year -------------------- */
            RealF added_transp_ratio = added_transp / SXW.ppt;
            //add added_transp to the average. This is technically part of the current year, so no need to subtract anything.
            transp_window.average += added_transp/transp_window.size;
            //add added_transp ratio to the ratio average. This is technically part of the current year, so no need to subtract anything.
            transp_window.ratio_average += (added_transp_ratio)/transp_window.size;
            //put the new transpiration in the window. Note: oldest_index has not been incremented, so it points to what was just added
            transp_window.transp[transp_window.oldest_index] += added_transp;
            //put the new ratio in the window. Note: oldest_index has not been incremented, so it points to what was just added
            transp_window.ratios[transp_window.oldest_index] += added_transp_ratio;
            // calculate the new sum of squares value
            RealF totalTranspRatio = (sumTranspTotal + added_transp)/SXW.ppt;
            RealF ssqr = (totalTranspRatio - transp_window.ratio_average) * (totalTranspRatio - transp_window.ratio_average);
            // Subtract the sum of squares calculated above, which was stored in the array. Replace it with what was just calculated.
            transp_window.sum_of_sqrs += ssqr - transp_window.SoS_array[transp_window.oldest_index];
            // replace the sum of squares with what we just calculated
            transp_window.SoS_array[transp_window.oldest_index] =  ssqr;

            //printf("Year %d: ratio = %f, mean = %f, sos = %f\n",Globals.currYear,
               //transp_ratio+added_transp_ratio,transp_window.ratio_average, transp_window.sum_of_sqrs);
            
            /* Adds the additional transpiration to the remaining transpiration 
             * so it can be distributed proportionally to the functional groups. */
            TranspRemaining += added_transp;
    }

    //oldest_index++ accounting for the array bounds
    transp_window.oldest_index = (transp_window.oldest_index + 1) % transp_window.size;

    /* ------------ End testing to see if additional transpiration is necessary ---------- */

    // If there is transpiration this year add the remaining (and added) transpiration to each functional group.
    // Else the sum of transpiration is 0 and for each group there is also zero transpiration.
    if (!ZRO(sumUsedByGroup)) {

        ForEachGroup(g) {
            use_by_group[g] += (use_by_group[g] / sumUsedByGroup) * TranspRemaining;
            //printf("for groupName= %s, after sum use_by_group[g]= %f \n",RGroup[g]->name,use_by_group[g]);
        }

        /*printf("'_transp_contribution_by_group': Group = %s, SXW.transp_SWA[g] = %f \n",
          RGroup[g]->name, SXW.transp_SWA[g]);
         */
    }
    
    //remember the last year. When setting up for a new iteration the same year will appear twice, and we want to skip it the second time
    lastYear = Globals.currYear; 
}
