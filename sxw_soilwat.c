/********************************************************/
/********************************************************/
/*  Source file: sxw_soilwat.c
/*  Type: module
 *
/*  Purpose: Subroutines to handle interactions with the
 *           SOILWAT model including setting up input/output
 *           files, writing to internal data structures,
 *           and actually running the model.
 *
/*  Called by: sxw.c
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
/*  History:
/*     (14-Apr-2002) -- INITIAL CODING - cwb
 *
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
 *		08/01/2012 - DLM - updated _update_productivity() function to use the 3 different VegProds now used in soilwat...
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
extern RealD *_roots_max,
             *_phen;          /* phenology read from file */
       RealF _prod_conv[MAX_MONTHS][3];

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


  _update_transp_coeff(sizes);
  _update_productivity();

#ifndef SXW_BYMAXSIZE
  SW_VPD_init();
#endif

}

void _sxw_sw_run(void) {
/*======================================================*/

  SW_Model.year = SW_Model.startyr + Globals.currYear -1;
  SW_CTL_run_current_year();

}

void _sxw_sw_clear_transp(void) {
/*======================================================*/

  Mem_Set(SXW.transp, 0, SXW.NPds * SXW.NTrLyrs * sizeof(RealD));

}


static void _update_transp_coeff(RealF relsize[]) {
/*======================================================*/
/* copy the relative root distribution to soilwat's layers */
/* POTENTIAL BUG:  if _roots_current is all zero but for some
 *   reason the productivity values (esp. %live) are >0,
 *   there can be transpiration but nowhere for it to come
 *   from in the soilwat model.
 */
  SW_LAYER_INFO *y ;
  GrpIndex g;
  LyrIndex t;
  RealF sum1=0., sum2=0., sum3=0.;

  ForEachTreeTranspLayer(t) { 
  	y = SW_Site.lyr[t];
    y->transp_coeff_tree = 0.;
    ForEachGroup(g)
    	if(getNTranspLayers(RGroup[g]->veg_prod_type == 1))
      		y->transp_coeff_tree += (RealF)_roots_max[Ilg(t,g)] * relsize[g];
    sum1 += y->transp_coeff_tree;
  }
  
  ForEachShrubTranspLayer(t) { 
  	y = SW_Site.lyr[t];
    y->transp_coeff_shrub = 0.;
    ForEachGroup(g)
    	if(getNTranspLayers(RGroup[g]->veg_prod_type == 2))
      		y->transp_coeff_shrub += (RealF)_roots_max[Ilg(t,g)] * relsize[g];
    sum2 += y->transp_coeff_shrub;
  }
  
  ForEachGrassTranspLayer(t) { 
  	y = SW_Site.lyr[t];
    y->transp_coeff_grass = 0.;
    ForEachGroup(g)
    	if(getNTranspLayers(RGroup[g]->veg_prod_type == 3))
      		y->transp_coeff_grass += (RealF)_roots_max[Ilg(t,g)] * relsize[g];
    sum3 += y->transp_coeff_grass;
  }

  /* normalize coefficients to 1.0 */
  ForEachTreeTranspLayer(t)
    SW_Site.lyr[t]->transp_coeff_tree /= sum1;
  ForEachShrubTranspLayer(t)
    SW_Site.lyr[t]->transp_coeff_shrub /= sum2;
  ForEachGrassTranspLayer(t)
    SW_Site.lyr[t]->transp_coeff_grass /= sum3;


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
  Months m;

  SW_VEGPROD *v = &SW_VegProd;
  RealF totbmass = 0.0,
        bmassg[MAX_RGROUPS] = {0.},
        cumprop=0, cumprop1=0, cumprop2=0, cumprop3=0,
        biomass=0, biomass1=0, biomass2=0, biomass3=0,
        props[MAX_MONTHS] = {0.},
        props1[MAX_MONTHS] = {0.},
        props2[MAX_MONTHS] = {0.},
        props3[MAX_MONTHS] = {0.};

#ifdef SXW_BYMAXSIZE
  #define Biomass(g)  _Grp_BMass[g]
#else
  #define Biomass(g)  RGroup_GetBiomass(g)
#endif

  /* get total biomass for the plot in sq.m */
  ForEachGroup(g) {
    bmassg[g] = Biomass(g) / Globals.plotsize;
    totbmass += bmassg[g];
  }

  /* compute monthly biomass, litter, and pct live per month */
  ForEachMonth(m) {
    if ( GT(totbmass, 0.) ) {
      ForEachGroup(g) {
        props[m] += _phen[Igp(g,m)] * bmassg[g];
        
        if(1 == RGroup[g]->veg_prod_type)				//tree
        	props1[m] += _phen[Igp(g,m)] * bmassg[g];
        else if(2 == RGroup[g]->veg_prod_type)			//shrub
        	props2[m] += _phen[Igp(g,m)] * bmassg[g];
        else if(3 == RGroup[g]->veg_prod_type)			//grass
        	props3[m] += _phen[Igp(g,m)] * bmassg[g];
      }
      props[m] /= totbmass;
      props1[m] /= totbmass;
      props2[m] /= totbmass;
      props3[m] /= totbmass;
      
      v->tree.pct_live[m] = pow(props1[m], 0.2);
      v->shrub.pct_live[m] = pow(props2[m], 0.2);
      v->grass.pct_live[m] = pow(props3[m], 0.2);

      cumprop += props[m];
      cumprop1 += props1[m];
      cumprop2 += props2[m];
      cumprop3 += props3[m];
      
      v->tree.biomass[m] = (totbmass + v->tree.pct_live[m] * cumprop1 * totbmass);
      v->tree.litter[m]  = (v->tree.biomass[m] * _prod_conv[m][PC_Litter]);
      v->shrub.biomass[m] = (totbmass + v->shrub.pct_live[m] * cumprop2 * totbmass);
      v->shrub.litter[m]  = (v->shrub.biomass[m] * _prod_conv[m][PC_Litter]);
      v->grass.biomass[m] = (totbmass + v->grass.pct_live[m] * cumprop3 * totbmass);
      v->grass.litter[m]  = (v->grass.biomass[m] * _prod_conv[m][PC_Litter]);
      
      biomass1 += v->tree.biomass[m];
      biomass2 += v->shrub.biomass[m];
      biomass3 += v->grass.biomass[m];
      
    } else {
      v->tree.pct_live[m] = 0.;
      v->tree.biomass[m]  = 0.;
      v->tree.litter[m]   = 0.;
      v->shrub.pct_live[m] = 0.;
      v->shrub.biomass[m]  = 0.;
      v->shrub.litter[m]   = 0.;
      v->grass.pct_live[m] = 0.;
      v->grass.biomass[m]  = 0.;
      v->grass.litter[m]   = 0.;
    }

  }
	
  biomass = biomass1 + biomass2 + biomass3;
  if(biomass == 0)
  	biomass = 1;
  v->fractionTree = (biomass1 / biomass);
  v->fractionShrub = (biomass2 / biomass);
  v->fractionGrass = (biomass3 / biomass);	
#undef Biomass
}

