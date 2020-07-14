/********************************************************/
/********************************************************/
/*  Source file: sxw_soilwat.c
 *  Type: module
 *  Purpose: Subroutines to handle interactions with the
 *           SOILWAT model including setting up input/output
 *           files, writing to internal data structures,
 *           and actually running the model.
 *  Called by: sxw.c
 *  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model. */
/*  History:
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
 *          to use the 3 different VegProds now used in soilwat */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include "sw_src/generic.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sw_src/Times.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/SW_Defines.h"
#include "sxw.h"
#include "sxw_module.h"
#include "sw_src/SW_Control.h"
#include "sw_src/SW_Model.h"
#include "sw_src/SW_Site.h"
#include "sw_src/SW_SoilWater.h"
#include "sw_src/SW_VegProd.h"
#include "sw_src/SW_Files.h"


/*************** Global Variable Declarations ***************/
/***********************************************************/
#include "sxw_vars.h"

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_VEGPROD SW_VegProd;
extern SXW_resourceType* SXWResources;

/*************** Local Function Declarations ***************/
/***********************************************************/
static void _update_transp_coeff(void);
static void _update_productivity(RealF size[]);


/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_sw_setup (RealF sizes[]) {
/*======================================================*/
	int doy, k;
	SW_VEGPROD *v = &SW_VegProd;

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

void _sxw_sw_run(void) {
/*======================================================*/
	SW_Model.year = SW_Model.startyr + Globals->currYear-1;
	SW_CTL_run_current_year();
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
    SW_LAYER_INFO *y;
    GrpIndex g;
    LyrIndex l;
    RealF sum[NVEGTYPES] = {0.};

    ForEachTreeTranspLayer(l) {
        y = SW_Site.lyr[l];
        y->transp_coeff[SW_TREES] = 0.;
        ForEachGroup(g) {
          if (RGroup[g]->veg_prod_type == SW_TREES)
              if (getNTranspLayers(SW_TREES))
                  y->transp_coeff[SW_TREES] += (RealF) SXWResources->_roots_max[Ilg(l, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
        }
        sum[SW_TREES] += y->transp_coeff[SW_TREES];
    }

    ForEachShrubTranspLayer(l) {
        y = SW_Site.lyr[l];
        y->transp_coeff[SW_SHRUB] = 0.;
        ForEachGroup(g)
        if (RGroup[g]->veg_prod_type == SW_SHRUB) {
            if (getNTranspLayers(SW_SHRUB))
                y->transp_coeff[SW_SHRUB] += (RealF) SXWResources->_roots_max[Ilg(l, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;

            /*printf("* lyr=%d, group=%s(%d), type=%d, tl=%d, rootmax=%f, relsize2=%f, trco=%f\n",
              l, RGroup[g]->name, g, RGroup[g]->veg_prod_type, getNTranspLayers(RGroup[g]->veg_prod_type),
              SXWResources->_roots_max[Ilg(l, g)], getRGroupRelsize(g), y->transp_coeff[SW_SHRUB]);
             */
        }
        sum[SW_SHRUB] += y->transp_coeff[SW_SHRUB];
    }

    ForEachGrassTranspLayer(l) {
        y = SW_Site.lyr[l];
        y->transp_coeff[SW_GRASS] = 0.;
        ForEachGroup(g) {
          if (RGroup[g]->veg_prod_type == SW_GRASS)
              if (getNTranspLayers(SW_GRASS))
                  y->transp_coeff[SW_GRASS] += (RealF) SXWResources->_roots_max[Ilg(l, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
        }
        sum[SW_GRASS] += y->transp_coeff[SW_GRASS];
    }

    ForEachForbTranspLayer(l) {
        y = SW_Site.lyr[l];
        y->transp_coeff[SW_FORBS] = 0.;
        ForEachGroup(g) {
          if (RGroup[g]->veg_prod_type == SW_FORBS)
              if (getNTranspLayers(SW_FORBS))
                  y->transp_coeff[SW_FORBS] += (RealF) SXWResources->_roots_max[Ilg(l, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
        }
        sum[SW_FORBS] += y->transp_coeff[SW_FORBS];
    }

    /* normalize coefficients to 1.0 If sum is 0, then the transp_coeff is also 0. */
    ForEachTreeTranspLayer(l) {
      if (!ZRO(sum[SW_TREES])) SW_Site.lyr[l]->transp_coeff[SW_TREES] /= sum[SW_TREES];
    }
    ForEachShrubTranspLayer(l) {
      if (!ZRO(sum[SW_SHRUB])) SW_Site.lyr[l]->transp_coeff[SW_SHRUB] /= sum[SW_SHRUB];
    }
    ForEachGrassTranspLayer(l) {
      if (!ZRO(sum[SW_GRASS])) SW_Site.lyr[l]->transp_coeff[SW_GRASS] /= sum[SW_GRASS];
    }
    ForEachForbTranspLayer(l) {
      if (!ZRO(sum[SW_FORBS])) SW_Site.lyr[l]->transp_coeff[SW_FORBS] /= sum[SW_FORBS];
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

    SW_VEGPROD *v = &SW_VegProd;
    RealF totbmass = 0.0,
            *bmassg,
    vegTypeBiomass[NVEGTYPES] = {0.};

    bmassg = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF), "_update_productivity");


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
