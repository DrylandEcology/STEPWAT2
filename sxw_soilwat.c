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
#include "generic.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "Times.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "SW_Defines.h"
#include "sxw.h"
#include "sxw_module.h"
#include "SW_Control.h"
#include "SW_Model.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Files.h"


/*************** Global Variable Declarations ***************/
/***********************************************************/
#include "sxw_vars.h"

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_VEGPROD SW_VegProd;

/*********** Local/Module Variable Declarations ************/
/***********************************************************/
extern RealD *_roots_max;
extern RealD _prod_litter[MAX_MONTHS];
extern RealD * _prod_bmass;
extern RealD * _prod_pctlive;

#ifdef SXW_BYMAXSIZE
extern RealF _Grp_BMass[];  /* added 2/28/03 */
#endif

/*************** Local Function Declarations ***************/
/***********************************************************/
static void _update_transp_coeff(RealF relsize[]);
static void _update_productivity(void);


/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_sw_setup (RealF sizes[]) {
/*======================================================*/
	int doy, k;
	SW_VEGPROD *v = &SW_VegProd;

        _update_productivity();
	_update_transp_coeff(sizes);

#ifndef SXW_BYMAXSIZE
  for (doy = 1; doy <= MAX_DAYS; doy++) {
    ForEachVegType(k) {
      v->veg[k].litter_daily[doy] = 0.;
      v->veg[k].biomass_daily[doy] = 0.;
      v->veg[k].pct_live_daily[doy] = 0.;
      v->veg[k].lai_conv_daily[doy] = 0.;
    }
  }

	SW_VPD_init();
#endif

}

void _sxw_sw_run(void) {
/*======================================================*/
	SW_Model.year = SW_Model.startyr + Globals.currYear-1;
	SW_CTL_run_current_year();
}

void _sxw_sw_clear_transp(void) {
/*======================================================*/
	int k;

	Mem_Set(SXW.transpTotal, 0, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
	ForEachVegType(k) {
		Mem_Set(SXW.transpVeg[k], 0, SXW.NPds * SXW.NSoLyrs * sizeof(RealD));
	}
}

static void _update_transp_coeff(RealF relsize[]) {
    /*======================================================*/
    /* copy the relative root distribution to soilwat's layers */
    /* POTENTIAL BUG:  if _roots_current is all zero but for some
     *   reason the productivity values (esp. %live) are >0,
     *   there can be transpiration but nowhere for it to come
     *   from in the soilwat model.
     */
    SW_LAYER_INFO *y;
    GrpIndex g;
    LyrIndex t;
    RealF sum[NVEGTYPES] = {0.};

    ForEachTreeTranspLayer(t) {
        y = SW_Site.lyr[t];
        y->transp_coeff[SW_TREES] = 0.;
        ForEachGroup(g)
        if (RGroup[g]->veg_prod_type == 1)
            if (getNTranspLayers(RGroup[g]->veg_prod_type))
                y->transp_coeff[SW_TREES] += (RealF) _roots_max[Ilg(t, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
        sum[SW_TREES] += y->transp_coeff[SW_TREES];
    }

    ForEachShrubTranspLayer(t) {
        y = SW_Site.lyr[t];
        y->transp_coeff[SW_SHRUB] = 0.;
        ForEachGroup(g)
        if (RGroup[g]->veg_prod_type == 2) {
            if (getNTranspLayers(RGroup[g]->veg_prod_type))
                y->transp_coeff[SW_SHRUB] += (RealF) _roots_max[Ilg(t, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;

            /*printf("* lyr=%d, group=%s(%d), type=%d, tl=%d, rootmax=%f, relsize1=%f, relsize2=%f, trco=%f\n",
              t, RGroup[g]->name, g, RGroup[g]->veg_prod_type, getNTranspLayers(RGroup[g]->veg_prod_type),
              _roots_max[Ilg(t, g)], relsize[g], RGroup[g]->relsize, y->transp_coeff[SW_SHRUB]);
             */
        }
        sum[SW_SHRUB] += y->transp_coeff[SW_SHRUB];
    }

    ForEachGrassTranspLayer(t) {
        y = SW_Site.lyr[t];
        y->transp_coeff[SW_GRASS] = 0.;
        ForEachGroup(g)
        if (RGroup[g]->veg_prod_type == 3)
            if (getNTranspLayers(RGroup[g]->veg_prod_type))
                y->transp_coeff[SW_GRASS] += (RealF) _roots_max[Ilg(t, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
        sum[SW_GRASS] += y->transp_coeff[SW_GRASS];
    }

    ForEachForbTranspLayer(t) {
        y = SW_Site.lyr[t];
        y->transp_coeff[SW_FORBS] = 0.;
        ForEachGroup(g)
        if (RGroup[g]->veg_prod_type == 4)
            if (getNTranspLayers(RGroup[g]->veg_prod_type))
                y->transp_coeff[SW_FORBS] += (RealF) _roots_max[Ilg(t, g)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
        sum[SW_FORBS] += y->transp_coeff[SW_FORBS];
    }

    /* normalize coefficients to 1.0 If sum is 0, then the transp_coeff is also 0. */
    ForEachTreeTranspLayer(t)
    if (!ZRO(sum[SW_TREES])) SW_Site.lyr[t]->transp_coeff[SW_TREES] /= sum[SW_TREES];
    ForEachShrubTranspLayer(t)
    if (!ZRO(sum[SW_SHRUB])) SW_Site.lyr[t]->transp_coeff[SW_SHRUB] /= sum[SW_SHRUB];
    ForEachGrassTranspLayer(t)
    if (!ZRO(sum[SW_GRASS])) SW_Site.lyr[t]->transp_coeff[SW_GRASS] /= sum[SW_GRASS];
    ForEachForbTranspLayer(t)
    if (!ZRO(sum[SW_FORBS])) SW_Site.lyr[t]->transp_coeff[SW_FORBS] /= sum[SW_FORBS];

    /*printf("'_update_transp_coeff': ShrubTranspCoef: ");
    ForEachSoilLayer(t) {
      printf("[%d] = %.4f - ", t, SW_Site.lyr[t]->transp_coeff[SW_SHRUB]);
    }
    printf("\n");
     */

}

static void _update_productivity(void) {
    /*======================================================*/
    /* must convert STEPPE's plot-level biomass to SOILWAT's
     * sq.meter - based biomass.
     *
     * 2/28/03 - cwb - as per the notes at the top of this file,
     *     RGroup_GetBiomass is being changed to
     *     SUM(RGroup[g]->mature_biomass)
     *
     * 12-Apr-2004 - cwb - update biomass, %Live, and litter with
     *     actual biomass if SXW_BYMAXSIZE defined or max biomass
     *     otherwise.
     */
    GrpIndex g;
    TimeInt m;
    IntUS k;

    SW_VEGPROD *v = &SW_VegProd;
    RealF totbmass = 0.0,
            bmassg[MAX_RGROUPS] = {0.},
    vegTypeBiomass[4] = {0.};

#ifdef SXW_BYMAXSIZE
#define Biomass(g)  _Grp_BMass[g]
#else
#define Biomass(g)  RGroup_GetBiomass(g)
#endif

    /* get total biomass for the plot in sq.m */
    ForEachGroup(g) {
        RGroup[g]->rgroupFractionOfVegTypeBiomass = 0.0;
        bmassg[g] = Biomass(g) / Globals.plotsize;
        totbmass += bmassg[g];
        if (1 == RGroup[g]->veg_prod_type) { //tree
            vegTypeBiomass[0] += bmassg[g];
        } else if (2 == RGroup[g]->veg_prod_type) { //shrub
            vegTypeBiomass[1] += bmassg[g];
        } else if (3 == RGroup[g]->veg_prod_type) { //grass
            vegTypeBiomass[2] += bmassg[g];
        } else if (4 == RGroup[g]->veg_prod_type) {
            vegTypeBiomass[3] += bmassg[g];
        }
    }

    ForEachGroup(g) {
        if (GT(vegTypeBiomass[RGroup[g]->veg_prod_type - 1], 0.))
            RGroup[g]->rgroupFractionOfVegTypeBiomass = bmassg[g] / vegTypeBiomass[RGroup[g]->veg_prod_type - 1];
        else
            RGroup[g]->rgroupFractionOfVegTypeBiomass = 0;
    }

    /* compute monthly biomass, litter, and pct live per month */
    ForEachMonth(m) {

        ForEachVegType(k) {
            v->veg[k].pct_live[m] = 0.;
            v->veg[k].biomass[m] = 0.;
            v->veg[k].litter[m] = 0.;
        }

        if (GT(totbmass, 0.)) {

            ForEachGroup(g) {
                if (1 == RGroup[g]->veg_prod_type) { //tree
                    v->veg[SW_TREES].pct_live[m] += _prod_pctlive[Igp(g, m)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
                    v->veg[SW_TREES].biomass[m] += _prod_bmass[Igp(g, m)] * bmassg[g];
                } else if (2 == RGroup[g]->veg_prod_type) { //shrub
                    v->veg[SW_SHRUB].pct_live[m] += _prod_pctlive[Igp(g, m)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
                    v->veg[SW_SHRUB].biomass[m] += _prod_bmass[Igp(g, m)] * bmassg[g];
                } else if (3 == RGroup[g]->veg_prod_type) { //grass
                    v->veg[SW_GRASS].pct_live[m] += _prod_pctlive[Igp(g, m)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
                    v->veg[SW_GRASS].biomass[m] += _prod_bmass[Igp(g, m)] * bmassg[g];
                } else if (4 == RGroup[g]->veg_prod_type) { //forb
                    v->veg[SW_FORBS].pct_live[m] += _prod_pctlive[Igp(g, m)] * RGroup[g]->rgroupFractionOfVegTypeBiomass;
                    v->veg[SW_FORBS].biomass[m] += _prod_bmass[Igp(g, m)] * bmassg[g];
                }
            }

            v->veg[SW_TREES].litter[m] = (vegTypeBiomass[0] * _prod_litter[m]);
            v->veg[SW_SHRUB].litter[m] = (vegTypeBiomass[1] * _prod_litter[m]);
            v->veg[SW_GRASS].litter[m] = (vegTypeBiomass[2] * _prod_litter[m]);
            v->veg[SW_FORBS].litter[m] = (vegTypeBiomass[3] * _prod_litter[m]);
        }
    }

    if (GT(totbmass, 0.)) {
        //if (ZRO(biomass))
        //	biomass = 1;
        v->veg[SW_TREES].cov.fCover = (vegTypeBiomass[0] / totbmass);
        v->veg[SW_SHRUB].cov.fCover = (vegTypeBiomass[1] / totbmass);
        v->veg[SW_GRASS].cov.fCover = (vegTypeBiomass[2] / totbmass);
        v->veg[SW_FORBS].cov.fCover = (vegTypeBiomass[3] / totbmass);
        //TODO: figure how to calculate bareground fraction.
        v->bare_cov.fCover = 0;
    } else {

        ForEachVegType(k) {
            v->veg[k].cov.fCover = 0.0;
        }
        v->bare_cov.fCover = 1;
    }
#undef Biomass
}
