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
  RealD * _roots_max,     /* read from root distr. file */
        * _roots_active_sum,

       /* rgroup by period */
        * _phen;          /* phenology read from file */

extern
  RealF _resource_pr[MAX_RGROUPS],  /* resource convertable to pr */
        _resource_cur[MAX_RGROUPS];

extern
  RealF _bvt;

//void _print_debuginfo(void);

/*************** Local Function Declarations ***************/
/***********************************************************/

static void _transp_contribution_by_group(RealF use_by_group[]);

//static void _SWA_contribution_by_group(RealF use_by_group[]);


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

	_sxw_update_root_tables(sizes);
	// add transpiration to `SXW.transp_SWA`:
	_transp_contribution_by_group(_resource_cur);
	// add SWA to `SXW.transp_SWA`:
  //_SWA_contribution_by_group(SXW.sum_dSWA_repartitioned);

	ForEachGroup(g)
	{
    _resource_cur[g] = SXW.transp_SWA[g];
    //printf("for groupName= %smresource_cur prior to multiplication: %f\n",RGroup[g]->name, _resource_cur[g]);
		_resource_cur[g] = _resource_cur[g] * _bvt;
    //printf("for groupName= %s, resource_cur post multiplication: %f\n\n",Rgroup[g]->name, _resource_cur[g]);
	}
/* _print_debuginfo(); */
}

void _sxw_update_root_tables( RealF sizes[] ) {
/*======================================================*/
/* sizes is a simple array that contains the groups'
 * actual biomass in grams in group order.
 *
 */

	GrpIndex g;
	LyrIndex l;
	TimeInt p;
	RealD x;
	int t,nLyrs;

	/* set some things to zero 4-Tree,Shrub,Grass,Forb*/
	Mem_Set(_roots_active_sum, 0, 4 * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));

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

	/* normalize the previous roots X phen table */
	/* the relative "activity" of a group's roots in a
	 * given layer in a given month is obtained by dividing
	 * the cross product by the totals from above */
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
	/* use_by_group is the amount of transpiration (cm) assigned
	 * to each STEPPE functional group.
	 * Must call _update_root_tables() before this.
	 * Compute each group's amount of transpiration from SOILWAT2
	 * based on its biomass, root distribution, and phenological
	 * activity. */

	GrpIndex g;
	TimeInt p;
	LyrIndex l;
	int t;
	RealD *transp;
	RealF sumUsedByGroup = 0., sumTranspTotal = 0., TranspRemaining = 0.;

	ForEachGroup(g) //Steppe functional group
	{
    use_by_group[g] = 0.;
		t = RGroup[g]->veg_prod_type-1;

		switch(t) {
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

        //Loops through each month and calculates amount of transpiration for each steppe functional group
        //according to whether that group has active living roots in each soil layer for each month
            ForEachTrPeriod(p)
            {
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
	ForEachTrPeriod(p)
	{
		for (t = 0; t < SXW.NSoLyrs; t++)
			sumTranspTotal += SXW.transpTotal[Ilp(t, p)];
	}

    TranspRemaining = sumTranspTotal - sumUsedByGroup;
    //printf(" sumTranspTotal=%f, sumUsedByGroup=%f  TranspRemaining=%f \n",sumTranspTotal,sumUsedByGroup,TranspRemaining);

    ForEachGroup(g)
    {
        if (!ZRO(use_by_group[g])) {
            use_by_group[g] += (use_by_group[g]/sumUsedByGroup) * TranspRemaining;
            //printf("for groupName= %s, after sum use_by_group[g]= %f \n",RGroup[g]->name,use_by_group[g] );

            SXW.transp_SWA[g] = use_by_group[g];
            //printf("for groupName= %s, SXW.transp_SWA[g] in transp= %f \n",RGroup[g]->name,SXW.transp_SWA[g]);
        }

        /*printf("'_transp_contribution_by_group': Group = %s, SXW.transp_SWA[g] = %f \n",
          RGroup[g]->name, SXW.transp_SWA[g]);
        */

    }
}

/*======================================================*/
/* use_by_group is the amount of plant available soil water (SWA,cm)
 * assigned to each STEPPE functional group.
 * Must call _update_root_tables() before this.
 * Compute each group's amount of SWA from SOILWAT2
 * based on its biomass, root distribution, and phenological
 * activity. */
/*
static void _SWA_contribution_by_group(RealF use_by_group[]) {

	GrpIndex g;
	TimeInt p;
	LyrIndex l;
	int t;
	RealF sumUsedByGroup = 0.;

	ForEachGroup(g) //Steppe functional group
	{
		use_by_group[g] = 0.;
		t = RGroup[g]->veg_prod_type-1;
        //printf("g, t || %d, %d\n", g, t);

		ForEachTrPeriod(p)
		{
			for (l = 0; l < SXW.NSoLyrs; l++) {
        //printf("%d,%d,%d\n", t,l,p);
				use_by_group[g] += (RealF) (_roots_active_rel[Iglp(g, l, p)] * SXW.sum_dSWA_repartitioned[Ivlp(t,l,p)]);
                //printf("for groupName= %s, layerIndex: %d, month: %d, in swa loop use_by_group[g]= %f \n",RGroup[g]->name,l,p,use_by_group[g]);
			}
		}
		//printf("for groupName= %s, use_by_group[g] in swa= %f \n",RGroup[g]->name,use_by_group[g]);

        sumUsedByGroup += use_by_group[g];
        //printf(" sumUsedByGroup in swa=%f \n",sumUsedByGroup);

        SXW.transp_SWA[g] += use_by_group[g];
        //printf("SXW.transp_SWA[%d]: %f\n", g, SXW.transp_SWA[g]);
	}
}
*/
