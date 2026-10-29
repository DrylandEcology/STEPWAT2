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

#include <stdlib.h>

#include "sw_src/include/generic.h"
#include "sw_src/include/rands.h"
#include "sw_src/include/filefuncs.h"
#include "sw_src/include/myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/include/SW_Defines.h"
#include "sxw.h" // externs `*SXWResources`, `transp_window`, `resource_rng`
#include "sxw_module.h"
#include "sxw_vars.h"
#include "sw_src/include/SW_Control.h"
#include "sw_src/include/SW_Site.h"
#include "sw_src/include/SW_SoilWater.h"
#include "sw_src/include/SW_VegProd.h"
#include "sw_src/include/SW_Files.h"
#include "sw_src/include/SW_Times.h"

static RealD *propuptake = NULL;

static unsigned int prop_nIters;
static unsigned int prop_nYears;
static unsigned int prop_nGroups;
static unsigned int prop_nLayers;
static unsigned int prop_nMonths;
static unsigned int nextYear = 1;

#define PROP_IDX(i, y, g, l, p) \
    ((((((size_t)(i) * prop_nYears + (y)) \
        * prop_nGroups + (g)) \
        * prop_nLayers + (l)) \
        * prop_nMonths + (p)))


/*************** Global Variable Declarations ***************/
/***********************************************************/
/* for steppe, see ST_globals.h */

//void _print_debuginfo(void);

/*************** Local Function Declarations ***************/
/***********************************************************/

static void _transp_contribution_by_group(RealF use_by_group[],RealD proportion_by_layer_month[]);

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

  RealF *sizes_live;
  RealD *proportion_by_layer_month;
  GrpIndex g;
  TimeInt p;

  sizes_live = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF),
                                   "_sxw_update_resource", &LogInfo);

	ForEachGroup(g)
	{
		sizes_live[g] = 0.;
//printf("_sxw_update_resource()RGroup Name= %s, RGroup[g]->regen_ok=%d \n ", RGroup[g]->name, RGroup[g]->regen_ok);
		if (!RGroup[g]->regen_ok)
			continue;
		sizes_live[g] = RGroup_GetBiomass(g) * RGroup[g]->live_biomass;
	}

    /* Update the active relative roots based on current biomass values */
	_sxw_update_root_tables(sizes_live);

    size_t nEntries =
    (size_t)Globals->grpCount *
    SXW->NTrLyrs *
    SXW->NPds;

 proportion_by_layer_month = (RealD *)Mem_Calloc(
    nEntries,
    sizeof(RealD),
    "_sxw_update_resource",
    &LogInfo
);

	/* Assign transpiration (resource availability) to each STEPPE functional group */
	_transp_contribution_by_group(SXWResources->_resource_cur,proportion_by_layer_month);
   
   if(Globals->currIter >0 && (Globals->currYear == nextYear)){ 
    nextYear++;
    if(nextYear > SuperGlobals.runModelYears) nextYear =1;
    ForEachGroup(g) {
    ForEachTrPeriod(p) {
        int nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);

        for (LyrIndex l = 0; l < nLyrs; l++) {

            size_t src = Iglp(g, l, p);
            size_t dst = PROP_IDX(
                        Globals->currIter - 1,
                        Globals->currYear - 1,
                        g, l, p
            );

    propuptake[dst] = proportion_by_layer_month[src]; 
        }
    }
    }
    }
        free(proportion_by_layer_month);
        
        /* Scale transpiration resources by a constant, bvt, to convert resources
         * (cm) to biomass that can be supported by those resources (g/cm) */
	ForEachGroup(g)
	{
//printf("%s: _resource_cur pre-multiplication = %f. Post multiplication = %f.\n",RGroup[g]->name, 
                                //_resource_cur[g], SXWResources->_resource_cur[g] * RGroup[g]->_bvt);
		SXWResources->_resource_cur[g] = SXWResources->_resource_cur[g] * RGroup[g]->_bvt;
	}

    free(sizes_live);
}

void _sxw_update_root_tables( RealF sizes_live[] ) {
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
				x = SXWResources->_rootsXphen[Iglp(g, l, p)] * sizes_live[g];
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

static void _transp_contribution_by_group(RealF use_by_group[],RealD proportion_by_layer_month[]) {
    /*======================================================*/
    /* use_by_group is the amount of transpiration (cm) assigned to each STEPPE
     * functional group. Must call _update_root_tables() before this.
     * Compute each group's amount of transpiration from SOILWAT2
     * based on its biomass, root distribution, and phenological
     * activity. This represents "normal" resources or transpiration each year. */

    GrpIndex g;
    TimeInt p;
    int nLyrs, t;
    LyrIndex l;
    RealD *transp, proportion_total_resources = 0.;
    RealF sumUsedByGroup = 0., sumTranspTotal = 0., TranspRemaining = 0.;
    // Initialize monthly transpiration proportions for all groups and soil layers.
    size_t nEntries =
        (size_t) Globals->grpCount *
        SXW->NTrLyrs *
        SXW->NPds;

    memset(
        proportion_by_layer_month,
        0,
        nEntries * sizeof(*proportion_by_layer_month)
    );

    ForEachGroup(g) //Steppe functional group
    {
        use_by_group[g] = 0.;
        transp = SXW->transpTotal;

        //Loops through each month and calculates amount of transpiration for each STEPPE functional group
        //according to whether that group has active living roots in each soil layer for each month
        /* Details: for each soil layer `l` and each month (trperiod) `p`, the
           sum of `_roots_active_rel[Iglp(g, l, p)]` across rgroups `g` must
           be 1; only if they sum to 1, can they be used as proper weights
           for `transp[Ilp(l, p)]` to be split up among rgroups */
        ForEachTrPeriod(p) {
            nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
            for (l = 0; l < nLyrs; l++) {
                RealD totalTransp = transp[Ilp(l, p)];
                RealD groupTransp =
                    SXWResources->_roots_active_rel[Iglp(g, l, p)]
                    * totalTransp;

                // Accumulate transpiration for the original annual group allocation.
                use_by_group[g] += (RealF) groupTransp;

                // Calculate the group's proportion of transpiration for this layer and month.
                if (!ZRO(totalTransp)) {
                    proportion_by_layer_month[Iglp(g, l, p)] =
                        groupTransp / totalTransp;
                }
            }
        }

        sumUsedByGroup += use_by_group[g];
        //printf(" sumUsedByGroup in transp=%f \n",sumUsedByGroup);
    }

    // Normalize monthly transpiration proportions within each iteration.
    // For each month-layer combination, proportions across all resource groups
    // should sum to 1 if at least one group has active relative roots.
    // Groups with zero proportions remain zero.
    // If all groups have zero proportions, no resources are allocated.
    ForEachTrPeriod(p) {
        for (l = 0; l < SXW->NTrLyrs; l++) {
            RealD sumProportions = 0.;

            ForEachGroup(g) {
                sumProportions +=
                    proportion_by_layer_month[Iglp(g, l, p)];
            }

            // Redistribute proportions only among groups with nonzero values.
            if (!ZRO(sumProportions)) {
                ForEachGroup(g) {
                    proportion_by_layer_month[Iglp(g, l, p)] /=
                        sumProportions;
                }
            }
        }
    }

    //Very small amounts of transpiration remain and not perfectly partitioned to functional groups.
    //This check makes sure any remaining transpiration is divided proportionately among groups.
    ForEachTrPeriod(p) {
        for (t = 0; t < SXW->NSoLyrs; t++)
            sumTranspTotal += SXW->transpTotal[Ilp(t, p)];
    }

    TranspRemaining = sumTranspTotal - sumUsedByGroup;

    // If there is transpiration this year add the remaining transpiration to each functional group.
    // Else the sum of transpiration is 0 and for each group there is also zero transpiration.
    if (!ZRO(sumUsedByGroup)) {

        ForEachGroup(g) {
            use_by_group[g] += (use_by_group[g] / sumUsedByGroup) * TranspRemaining;
            //printf("for groupName= %s, after sum use_by_group[g]= %f \n",RGroup[g]->name,use_by_group[g]);
        }
    }

    ForEachGroup(g) {
        if (!ZRO(sumUsedByGroup) || use_by_group[g] != 0) {
            // Proportion of total transpiration that this group receives
            proportion_total_resources = use_by_group[g] / sumUsedByGroup;
        }
    }
}

void init_propuptake(unsigned int nIters, unsigned int nYears){

    prop_nIters = nIters;
    prop_nYears = nYears;
    prop_nGroups = Globals->grpCount;
    prop_nLayers = SXW->NTrLyrs;
    prop_nMonths = SXW->NPds;

    size_t total =
        (size_t)nIters * nYears *
        prop_nGroups * prop_nLayers * prop_nMonths;

    propuptake = calloc(total, sizeof(*propuptake));
}

void write_propuptake_csv(void)
{
    FILE *fp = fopen(Parm_name(F_PROPUPTAKE), "w");
    if (!fp) return;

    fprintf(fp,
        "Year,Month,Layer,Rgroup,propuptake_mean,propuptake_std\n");

    for (unsigned int y = 0; y < prop_nYears; y++) {
        for (unsigned int p = 0; p < prop_nMonths; p++) {
            for (unsigned int l = 0; l < prop_nLayers; l++) {
                for (unsigned int g = 0; g < prop_nGroups; g++) {

                    double mean = 0.;
                    double sos = 0.;

                    for (unsigned int i = 0; i < prop_nIters; i++) {
                        double value =
                            propuptake[PROP_IDX(i, y, g, l, p)];

                        double oldMean = mean;

                        mean = get_running_mean(i + 1, mean, value);

                        sos += get_running_sqr(
                            oldMean, mean, value
                        );
                    }

                    double std =
                        final_running_sd(prop_nIters, sos);

                    fprintf(fp, "%u,%u,%u,%s,%.10f,%.10f\n",
                            y + 1, p + 1, l + 1,
                            RGroup[g]->name,
                            mean, std);
                }
            }
        }
    }

    fclose(fp);
    free(propuptake);
}
