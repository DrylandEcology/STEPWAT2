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

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_SOILWAT SW_Soilwat;
extern SW_VEGPROD SW_VegProd;


/*************** Local Variable Declarations ***************/
/***********************************************************/
/* malloc'ed and maybe read in sxw.c but used here */
/* ----- 3d arrays ------- */
extern
  RealD * _roots_phen_lyr_grp, /* relative roots X phen by layer & group */
        * _roots_rel_phen_lyr_grp, /*relative to the total roots_phen_lyr_group */
      * _transp_base_use_by_group,
      * _transp_curr_use_by_group;


/* ----- 2D arrays ------- */

       /* rgroup by layer */
extern
  RealD * _roots_max,     /* read from root distr. file */
        * _roots_current, /* relsize * max root */
        * _roots_rel,     /* group's proportion  */
        * _roots_phen_totals, /* totals, rel roots X rel phen per group/layer */

       /* rgroup by period */
        * _phen,          /* phenology read from file */
        * _phen_grp_rel,      /* normalized to 1.0 per group */
        * _transp_grp_totals; /* group's contribution to year's transp */

extern RealD *_transp_base; /* nominal values read from file */

extern
  RealF _resource_equ[MAX_RGROUPS],  /* equilibrium resource utilization */
        _resource_cur[MAX_RGROUPS];  /* current resource utilization */

void _print_debuginfo(void);

/*************** Local Function Declarations ***************/
/***********************************************************/

static void _transp_contribution_by_group(
               RealD *transp,
               RealF use_by_group[]);


/***********************************************************/
/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_relative_phen(void) {
/*======================================================*/
/* should only be called once, after phenology table is read */

  RealD *sums;
  GrpIndex g;
  TimeInt p;

  sums = (RealD *) Mem_Calloc(SXW.NPds, sizeof(RealD),
                              "_sxw_phen_rel()");

  /* get totals by period */
  ForEachGroup(g) {
    ForEachTrPeriod(p) {
      sums[p] += _phen[Igp(g,p)];
    }
  }

  /* compute relative phenological activity per group per period */
  ForEachTrPeriod(p) {
    ForEachGroup(g)
      _phen_grp_rel[Igp(g,p)] = !ZRO(sums[p])
                              ? _phen[Igp(g,p)] / sums[p]
                              : 0.0;
  }


  Mem_Free(sums);
}

void _sxw_base_resource_vector(void) {
/*======================================================*/
/* computes the equilibrium resource demand against which
 * current resource use is compared for PR value.
 */
/* call once after equilibrium transpiration is read. */

  RealF sizes[MAX_RGROUPS];
  GrpIndex g;

  ForEachGroup(g) sizes[g] = 1.0;

  _sxw_update_root_tables(sizes);

  _transp_contribution_by_group( _transp_base,
                                 _resource_equ);
/* _print_debuginfo(); */
}

void _sxw_update_resource(void) {
/*======================================================*/
/* call after each run of SOILWAT to compute resource vector */
/* computing the groups' transpirations gives the amount
 * of resource available.  The equilibrium transpiration
 * (previously computed) is the amount of resource required
 * at full stocking or normal conditions. However, the
 * requirements must be scaled to the size of the group so
 * for example a small usage for a small group is not considered
 * inadequate resources compared to "typical".
 *
 * 2/28/03 - cwb - added assignment to res_avail and res_req
 *         because it's needed in rgroup_ResPartIndiv(). check
 *         there for more info.
 * 3/25/2003 - cwb - scaling requirements by 1/density accounts
 *         for rapidity with which a plant uses the resources of
 *         the plot.
 * 15-May-03 (cwb) - improved resource partitioning algorithm
 *         in rgroup_PartResource() obviates the need here.
 */

/*
  GrpIndex rg;
  GroupType *g;
*/
IntUS l,p;
RealD sum=0.;

    ForEachTreeTranspLayer(l) {
      ForEachTrPeriod(p) sum += SXW.transp[Ilp(l,p)];
    }


  _transp_contribution_by_group( SXW.transp,
                                 _resource_cur);
/*
required and available now get computed by _rgroup_PartResource()

  ForEachGroup(rg) {  g = RGroup[rg];
    g->res_avail    = _resource_cur[rg];
    g->res_required = _resource_equ[rg]
                      * (g->relsize / g->max_density);
    g->pr = ZRO(g->res_avail)
                  ? 1.
                  : g->res_required / g->res_avail;
  }
*/

}

void _sxw_update_root_tables( RealF relsizes[] ) {
/*======================================================*/
/* relsizes is a simple array that contains the groups'
 * relative size values in the actual group order
 *
 * 13-May-03 - adding code to factor in groups' resource space.
 *             this is intended as a temporary measure to better
 *             account for competition, but I'll be surprised
 *             if it doesn't stay in.
 */

  GrpIndex g;
  LyrIndex l;
  TimeInt p;
  RealD x;
  RealD sums[MAX_LAYERS];

  /* set some things to zero */
  ForEachTreeTranspLayer(l) sums[l] = 0.;
  Mem_Set( _roots_phen_totals, 0,
           SXW.NGrps * SXW.NPds * SXW.NTrLyrs * sizeof(RealD));

 /* actual root sizes */
  ForEachGroup(g) {
    ForEachTreeTranspLayer(l) {
      x = (RealF) (relsizes[g] * RGroup[g]->min_res_req)
          * _roots_max[Ilg(l,g)];

      _roots_current[Ilg(l,g)] = x;
      sums[l] += x;
    }
  }

  /* a group's roots in a layer as a proportion of the total
   * roots in the layer */
  ForEachGroup(g) {
    ForEachTreeTranspLayer(l)
      _roots_rel[Ilg(l,g)] = ZRO(sums[l]) ? 0.
                           : _roots_current[Ilg(l,g)] / sums[l];
  }

  /* update the 3d tables */
  /* multiply proportional roots X normalized phenology
   * to get cross product.  Total across groups for later
   * normalization */
  ForEachGroup(g) {
    ForEachTreeTranspLayer(l) {
      ForEachTrPeriod(p) {
        x = _roots_rel[Ilg(l,g)] * _phen_grp_rel[Igp(g,p)];
        _roots_phen_lyr_grp[Iglp(g,l,p)] = x;
        _roots_phen_totals[Ilp(l,p)] += x;
      }
    }
  }

  /* normalize the previous roots X phen table */
  /* the relative "activity" of a group's roots in a
   * given layer in a given month is obtained by dividing
   * the cross product by the totals from above */
  ForEachGroup(g) {
    ForEachTreeTranspLayer(l) {
      ForEachTrPeriod(p)
        _roots_rel_phen_lyr_grp[Iglp(g,l,p)] =
           ZRO(_roots_phen_totals[Ilp(l,p)])
           ? 0.
           : _roots_phen_lyr_grp[Iglp(g,l,p)]
             / _roots_phen_totals[Ilp(l,p)];
    }
  }

}


static void _transp_contribution_by_group(
               RealD *transp,
               RealF use_by_group[]) {
/*======================================================*/
/* transp is usually SXW.transp but the first time is
 *        _transp_base,
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

  ForEachGroup(g) {
    use_by_group[g] = 0.;  /* clear */
    ForEachTrPeriod(p) {
      ForEachTreeTranspLayer(l) {
        use_by_group[g] += (RealF)
          (_roots_rel_phen_lyr_grp[Iglp(g,l,p)]
          * transp[Ilp(l,p)]);
      }
    }
  }

}



