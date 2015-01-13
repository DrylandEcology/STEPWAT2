/********************************************************/
/********************************************************/
/*  Source file: sxw_resource.c
 *
/*  Type: module
 *
/*  Purpose: Compute resource vector for STEPPE based on
 *           transpiration values from SOILWAT.
 *
/*  Dependency:  sxw.c
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *
/*  History:
/*     (21-May-2002) -- INITIAL CODING - cwb
 *     19-Jun-2003 - cwb - Added RealD (double precision)
 *                 types for internal dynamic matrices and
 *                 other affected variables.  See notes in
 *                 sxw.c.
/*
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
#include "SW_Model.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Files.h"
#include "SW_Times.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
/* for steppe, see ST_globals.h */

//extern SW_SITE SW_Site;
//extern SW_MODEL SW_Model;
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
        _resource_cur[MAX_RGROUPS];  /* current resource utilization */

extern
  RealF _bvt;

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

  RealF sizes[MAX_RGROUPS] = {0.};
  GrpIndex g;
  SppIndex sp;
  int i;

  #ifdef SXW_BYMAXSIZE
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
		if (!RGroup[g]->regen_ok)
			continue;
		if (RGroup[g]->max_age == 1) {
			ForEachGroupSpp(sp,g,i)
				sizes[g] += Species[sp]->mature_biomass * .75;
		} else {
			sizes[g] = RGroup_GetBiomass(g);
		}
	}
  #endif

	_sxw_update_root_tables(sizes);
	_transp_contribution_by_group(_resource_cur);

	ForEachGroup(g)
		_resource_pr[g] =
				ZRO(sizes[g]) ? 0.0 : _resource_cur[g] * _bvt / sizes[g];
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
				_roots_active_sum[Itlp(t,l, p)] += x;
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
	/*
	 * use_by_group is the vector to be used in the resource
	 *        availability calculation, ie, the output.

	 * must call _update_root_tables() before this.
	 *
	 */

	/* compute each group's contribution to the
	 * transpiration values retrieved from SOILWAT based
	 * on its relative size, its root distribution, and
	 * its phenology (activity).
	 */

	GrpIndex g;
	TimeInt p;
	LyrIndex l;
	int t;
	RealD *transp;
	RealF sumUsedByGroup = 0., sumTranspTotal = 0.;

	ForEachGroup(g)
	{
		use_by_group[g] = 0.; /* clear */
		t = RGroup[g]->veg_prod_type-1;
		switch(t) {
		case 0://Tree
			transp = SXW.transpTrees;
			break;
		case 1://Shrub
			transp = SXW.transpShrubs;
			break;
		case 2://Grass
			transp = SXW.transpGrasses;
			break;
		case 3://Forb
			transp = SXW.transpForbs;
			break;
		default:
			transp = SXW.transpTotal;
			break;
		}
		ForEachTrPeriod(p)
		{
			int nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
			for (l = 0; l < nLyrs; l++) {
				use_by_group[g] += (RealF) (_roots_active_rel[Iglp(g, l, p)] * transp[Ilp(l, p)]);
			}
		}
		sumUsedByGroup += use_by_group[g];
	}
	//Extra Transp is proportionately given to the Rgroups.
	ForEachTrPeriod(p)
	{
		for (t = 0; t < SXW.NSoLyrs; t++)
			sumTranspTotal += SXW.transpTotal[Ilp(t, p)];
	}

	ForEachGroup(g)
	{
		use_by_group[g] = (use_by_group[g]/sumUsedByGroup) * sumTranspTotal;
	}
}



