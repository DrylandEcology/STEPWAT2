/********************************************************/
/********************************************************/
/*	Source file: SW_Flow_lib.cif(toDebugSoilTemp)
	Type: module
	Application: SOILWAT - soilwater dynamics simulator
	Purpose: Water flow subroutines that can be used
	        as a more or less independent library of
	        soil water flow routines.  These routines
	        are designed to operate independently of
	        the soilwater model's data structures.
	
	See Also: SoWat_flow.c  SoWat_flow_subs.c
	
	History:
	(4/10/2000) -- INITIAL CODING - cwb
	10/19/2010	(drs) added function hydraulic_redistribution()
	11/16/2010	(drs) added function forest_intercepted_water()
						renamed evap_litter_stdcrop() -> evap_litter_veg()
	01/06/2011	(drs) drainout is calculated in infiltrate_water_high() and infiltrate_water_low(), but was not added up -> fixed it by changing (*drainout) = drainlw to (*drainout) += drainlw in function infiltrate_water_low()
	01/06/2011	(drs) drainlw was not set to 0. in the for-loop, i.e. a layer with swc below pwp inherited incorrectly drainlw from layer above instead of drainlw=0
	02/22/2011	(drs) evap_litter_veg() set aet implicitely to 0. and then added litter evap; changed such that litter evap is now added to whatever value aet has previously
	07/21/2011	(drs) infiltrate_water_high() & infiltrate_water_low();
						- added variables: saturated water content - swcsat, impermeability, water standing above soil - standingWater
						- percolation is adjusted by (1-impermeability)
						- if a lower soil layer becomes saturated then water is pushed back up even above the soil surface into standingWater
	07/22/2011	(drs) included evaporation from standingWater into evap_litter_veg_surfaceWater() previously called evap_litter_veg()
	07/22/2011	(drs) included evaporation from standingWater into reduce_rates_by_surfaceEvaporation() previously called reduce_rates_by_intercepted()
	08/22/2011	(drs) renamed bs_ev_tr_loss() to EsT_partitioning()
	08/22/2011	(drs) added EsT_partitioning_forest() to partition bare-soil evaporation (Es) and transpiration (T) in forests
	09/08/2011	(drs) renamed stcrp_intercepted_water() to grass_intercepted_water();
						added shrub_intercepted_water() as current copy from grass_intercepted_water();
						added double scale to each xx_intercepted_water() function to scale for snowdepth effects and vegetation type fraction
						added double a,b,c,d to each xx_intercepted_water() function for for parameters in intercepted rain = (a + b*veg) + (c+d*veg) * ppt
	09/09/2011	(drs) added xx_EsT_partitioning() for each vegetation type (grass, shrub, tree); added double lai_param as parameter in exp-equation
	09/09/2011	(drs) replaced evap_litter_veg_surfaceWater() with evap_fromSurface() to be called for each surface water pool seperately
	09/09/2011	(drs) replaced SHD_xx constanst with input parameters in pot_transp()
	09/09/2011	(drs) added double Es_param_limit to pot_soil_evap() to scale and limit bare-soil evaporation rate (previously hidden constant in code)
	09/09/2011	(drs) renamed reduce_rates_by_surfaceEvaporation() to reduce_rates_by_unmetEvapDemand()
	09/11/2011	(drs) added double scale to hydraulic_redistribution() function to scale for vegetation type fraction
	09/13/2011	(drs)	added double scale, a,b,c,d to each litter_intercepted_water() function for parameters in litter intercepted rain = ((a + b*blitter) + (c+d*blitter) * ppt) * scale
	09/21/2011	(drs)	reduce_rates_by_unmetEvapDemand() is obsolete, complete E and T scaling in SW_Flow.c
	01/28/2012	(drs)	unsaturated flow drains now to swcmin instead of swcwp: changed infiltrate_water_low() accordingly
	02/04/2012	(drs)	in 'remove_from_soil': control that swc_avail >= 0
	03/23/2012	(drs)	excluded hydraulic redistribution from top soil layer (assuming that this layer is <= 5 cm deep)
	03/26/2012	(drs)	infiltrate_water_low(): allow unsaturated percolation in first layer
	05/24/2012  (DLM) started coding soil_temperature function
	05/25/2012  (DLM)  added all this fun crazy linear regression stuff to soil_temperature function
	05/29/2012  (DLM) still working on soil_temperature function, linear regression stuff should work now
	05/30/2012  (DLM) got rid of nasty segmentation fault error in soil_temperature function, also tested math seems correct after checking by hand.  added the ability to change the value of deltaX
	05/30/2012  (DLM) added soil_temp_error variable... it simply keeps track of whether or not an error has been reported in the soil_temperature function.  0 for no, 1 for yes.
	05/30/2012  (DLM) added # of lyrs check & maxdepth check at the beginning of soil_temperature function to make sure code doesn't blow up... if there isn't enough lyrs (ie < 2) or the maxdepth is too little (ie < deltaX * 2), the function quits out and reports an error to the user
	05/31/2012  (DLM) added theMaxDepth variable to soil_temperature function to allow the changing of the maxdepth of the equation, also soil_temperature() function now stores most regression data in a structure to reduce redundant regression calculations & speeds things up
	05/31/2012  (DLM) added soil_temperature_init() function for use in the soil_temperature function to initialize the values for use in the regressions... it is not apart of the header file, because it's not meant to be an external function
	06/01/2012  (DLM) edited soil_temperature function(), changed deltaT variable from hours to seconds, also changed some of the regression calculations so that swc, fc, & wp regressions are scaled properly... results are actually starting to look usable!
	06/13/2012  (DLM) soil_temperature function no longer extrapolates values for regression layers that are out of the bounds of the soil layers... instead they are now set to the last soil layers values.  extrapolating code is still in the function and can be commented out and used if wishing to go back to extrapolating the values...
*/
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "generic.h"
#include "SW_Defines.h"
#include "SW_Flow_lib.h"
#include "SW_Flow_subs.h"
#include "Times.h"

/* *************************************************** */
/*                Module-Level Variables               */
/* --------------------------------------------------- */

unsigned int soil_temp_error = 0;  // simply keeps track of whether or not an error has been reported in the soil_temperature function.  0 for no, 1 for yes.
unsigned int soil_temp_init = 0;   // simply keeps track of whether or not the regression values for the soil_temperature function have been initialized.  0 for no, 1 for yes.

ST_RGR_VALUES stValues; // keeps track of the regression values, for use in soil_temperature function

/* *************************************************** */
/* *************************************************** */
/*              Local Function Definitions             */
/* --------------------------------------------------- */


void grass_intercepted_water(double *pptleft, double *wintgrass,
								double ppt,      double vegcov,	double scale,
								double a, double b, double c, double d)  {
/**********************************************************************
PURPOSE: Calculate the water intercepted by grasses.

HISTORY:
4/30/92  (SLC)
7/1/92   (SLC) Reset pptleft to 0 if less than 0 (due to round off)
1/19/93  (SLC) Check if vegcov is zero (in case there was no biomass),
then no standing crop interception is possible.
15-Oct-03 (cwb) replaced Parton's original equations with new ones
developed by John Bradford based on Corbet and Crouse 1968.
Replaced the following code:
par1 = LE(vegcov, 8.5) ?  0.9 + 0.04 * vegcov
: 1.24 + (vegcov-8.5) *.35;
par2 = LE(vegcov, 3.0) ? vegcov * .33333333
: 1. + (vegcov-3.)*0.182;
*wintstcr = par1 * .026 * ppt + 0.094 * par2;

21-Oct-03 (cwb) added MAX_WINTLIT line

INPUTS:
ppt     - precip. for the day
vegcov  - vegetation cover for the day (based on monthly biomass
values, see the routine "initprod")

OUTPUT:
pptleft -  precip. left after interception by standing crop.
wintstcr - amount of water intercepted by standing crop.
**********************************************************************/
	double  intcpt, slope;
	
	if ( GT(vegcov, 0.)  &&  GT(ppt, 0.) ) {
		intcpt = b * vegcov + a;
		slope  = d * vegcov + c;
		
		*wintgrass = (intcpt + slope * ppt) * scale;
		
		*wintgrass = fmin(*wintgrass, ppt);
		*wintgrass = fmin(*wintgrass, MAX_WINTSTCR);
		*pptleft = fmax( ppt - *wintgrass, 0.0);
	} else { /*  no precip., so obviously nothing is intercepted by standing crop. */
		*pptleft=ppt;
		*wintgrass=0.0;
	}
}

void shrub_intercepted_water(double *pptleft, double *wintshrub,
								double ppt,      double vegcov,	double scale,
	                            double a, double b, double c, double d)  {
/**********************************************************************
PURPOSE: Calculate the water intercepted by shrubs
**********************************************************************/
	double  intcpt, slope;
	
	if ( GT(vegcov, 0.)  &&  GT(ppt, 0.) ) {
		intcpt = b * vegcov + a;
		slope  = d * vegcov + c;
		
		*wintshrub = (intcpt + slope * ppt) * scale;
		
		*wintshrub = fmin(*wintshrub, ppt);
		*wintshrub = fmin(*wintshrub, MAX_WINTSTCR);
		*pptleft = fmax( ppt - *wintshrub, 0.0);
	} else { /*  no precip., so obviously nothing is intercepted by standing crop. */
		*pptleft=ppt;
		*wintshrub=0.0;
	}
}


void tree_intercepted_water(	double *pptleft, double *wintfor,
								double ppt,      double LAI,	double scale,
								double a, double b, double c, double d){
/**********************************************************************
PURPOSE: Calculate water intercepted by forests

HISTORY:
11/16/2010	(drs)

INPUTS:
ppt     - precip. for the day in cm
LAI	- forest LAI in cm/cm
scale - scale interception with fraction of tree vegetation component or with snowdepth-scaler

OUTPUT:
pptleft -  precip. left after interception by forest in cm.
wintfor - amount of water intercepted by forest in cm.
**********************************************************************/

	double  intcpt, slope;
	
	if ( GT(LAI, 0.)  &&  GT(ppt, 0.) ) {
		intcpt = b * LAI + a;
		slope  = d * LAI + c;
	
		*wintfor = (intcpt + slope * ppt) * scale;
	
		*wintfor = fmin(*wintfor, ppt);
		*wintfor = fmin(*wintfor, MAX_WINTFOR);
		*pptleft = fmax( ppt - *wintfor, 0.0);
	} else {	/*  no precip., so obviously nothing is intercepted by forest. */
		*pptleft=ppt;
		*wintfor=0.0;
	}
}



void litter_intercepted_water(	double *pptleft, double *wintlit,
								double blitter,	double scale,
								double a, double b, double c, double d) {
/**********************************************************************
PURPOSE: Calculate water intercepted by litter

HISTORY:
4/30/92  (SLC)
7/1/92   (SLC) Reset pptleft to 0 if less than 0 (due to round off)
6-Oct-03 (cwb) wintlit = 0 if no litter.
15-Oct-03 (cwb) replaced Parton's original equations with new ones
developed by John Bradford based on Corbet and Crouse, 1968.
Replaced the following code:
par1 = exp((-1. + .45 * log10(blitter+1.)) * log(10.));
*wintlit = (.015 * (*pptleft) + .0635) * exp(par1);

21-Oct-03 (cwb) added MAX_WINTLIT line

INPUTS:
blitter - biomass of litter for the day

OUTPUTS:
pptleft -  precip. left after interception by litter.
wintlit  - amount of water intercepted by litter .
**********************************************************************/

	double intcpt, slope;
	
	if ( ZRO(blitter) ) {
		*wintlit = 0.0;	
	} else if ( GT(*pptleft,  0.0)  ) {
		intcpt = b * blitter + a;
		slope  = d * blitter + c;
		
		*wintlit = (intcpt + slope * (*pptleft)) * scale;
		
		*wintlit = fmin(*pptleft,*wintlit);
		*wintlit = fmin(*wintlit, MAX_WINTLIT);
		*pptleft -= *wintlit;
		*pptleft = fmax(*pptleft, 0.0);
		
	} else {
		*pptleft = 0.0;
		*wintlit = 0.0;
	}
}

void infiltrate_water_high( double swc[],     double drain[],
double *drainout, double pptleft,
unsigned int nlyrs,     double swcfc[], double swcsat[], double impermeability[], double *standingWater ) {
/**********************************************************************
PURPOSE: Infilitrate water into soil layers under high water
conditions.

HISTORY:
4/30/92  (SLC)
1/14/02 - (cwb) fixed off by one error in loop.
10/20/03 - (cwb) added drainout variable to return drainage
out of lowest layer

INPUTS:
swc - soilwater content before drainage.
swcfc    - soilwater content at field capacity.
pptleft  - precip. available to the soil.
nlyrs - number of layers to drain from

OUTPUTS:
drain  - drainage from layers
swc_local - soilwater content after water has been drained
**********************************************************************/

unsigned int  i;
int j;
double d[nlyrs];
double push;

	swc[0] += pptleft;
	(*standingWater) = 0.;
	
	for (i=0; i < nlyrs; i++) {
		/* calculate potential saturated percolation */
		d[i]  = fmax(0., (1. - impermeability[i]) * (swc[i] - swcfc[i]) );
		drain[i] = d[i];
		
		if(i < nlyrs-1){ /* percolate up to next-to-last layer */
			swc[i+1] += d[i];
			swc[i]   -= d[i];
		} else { /* percolate last layer */
			(*drainout) = d[i];
			swc[i] -= (*drainout);		
		}
	}


/* adjust (i.e., push water upwards) if water content of a layer is now above saturated water content */
	for (j=nlyrs; j >= 0; j--) {
		if( GT(swc[j], swcsat[j]) ) {
			push = swc[j] - swcsat[j];
			swc[j] -= push;
			if(j > 0) {
				drain[j-1] -= push;
				swc[j-1] += push;
			} else {
				(*standingWater) = push;
			}
		}
	}
}


double petfunc(unsigned int doy,        double avgtemp,   double rlat,
double reflec,   double humid,     double windsp,
double cloudcov, double transcoeff ) {
/***********************************************************************
PURPOSE: Calculate the potential evapotranspiration rate using
pennmans equation (1948)

HISTORY:
4/30/92  (SLC)

INPUTS: (GLOBAL inputs from soilw.inc)
Time.Model:
doy            - current day
month           - current month.
sky_parms:
r_humidity(month)   - average relative humidity for the month. (%)
windspeed(month)   - average wind speed for the month. (mph)
cloudcov(month)   - average cloud cover for the month. (%)
transmission(month) - transmission coefficient for the month
weather:
temp_avg         - average temperature for the day
site_parm:
reflection          - reflectivity
latitude       - latitude of the site (in radians)

LOCAL VARIABLES:
solrad - solar radiation (ly/day)
declin - declination (radians)
ahou   - ?
shwave - short wave solar radiation
kelvin - kelvin degrees
arads  - ?
clrsky - relative amount of clear sky
fhumid -
ftemp  -
par1,par2 - parameters in computation of pet.

***********************************************************************/

double   declin,par1,par2,ahou,solrad,
		shwave,kelvin,arads,clrsky,ftemp,
		vapor, result;

vapor = svapor(avgtemp);

/*   calculate the short wave solar radiation */
/*   on a clear day using a equation presented by Sellers(1965) */
declin = .401426 *sin(6.283185 *(doy -77.) /365.);
par2 = -tan(rlat) * tan(declin);
par1 = sqrt(1. - (par2*par2));
ahou = fmax(atan2(par1,par2), 0.0);
solrad = 917. * transcoeff
	* (ahou * sin(rlat) *sin(declin)
	+ cos(rlat) *cos(declin) *sin(ahou));
shwave = solrad *.0168 / transcoeff;

/* calculate the PET */
kelvin = avgtemp + 273.;
arads = vapor *3010.21 / (kelvin*kelvin) ;
clrsky = 1. - cloudcov/100.;
humid *= vapor/100.;
ftemp = kelvin *.01;
ftemp = ftemp*ftemp*ftemp*ftemp *.201;
par1 = .35 *(vapor -humid)
	* (1. +.0098 *(windsp*24.));
par2 = shwave * (1.- reflec) * (.18 +.55 * clrsky)
	- ftemp * (.56-.092 * sqrt(humid)) * (.10+.90*clrsky);
result = ((arads*par2 + .27*par1)/(arads+.27)) /10.;

return fmax(result, 0.01);

}

double svapor(double temp) {
/*********************************************************************
PURPOSE: calculate the saturation vapor pressure of water
the clausius-clapeyron equation(hess,1959) is used
HISTORY:
4/30/92  (SLC)

INPUTS:
atemp - average temperature for the day

OUTPUT:
svapor - saturation vapor pressure (mm of hg)

*********************************************************************/
double   par1,par2;

par1 = 1. / (temp+273.);
par2 = log(6.11) + 5418.38 * (.00366-par1);

return (exp(par2)*.75);
}

void transp_weighted_avg( double *swp_avg,   unsigned int   n_tr_rgns,
unsigned int   n_layers,   unsigned int   tr_regions[],
double tr_coeff[], double swc[]) {
/**********************************************************************

PURPOSE: Compute weighted average of soilwater potential to be
used for transpiration calculations.

HISTORY:
Original:
4/30/92  SLC
6/9/93   (SLC) check that (sum_tr_co <> 0) before dividing swp by this
number
4/10/2000 CWB -- began recoding in C
9/21/01   cwb -- adjusted method for determining transpiration
regions to reflect the new design.  removed
tr_reg_min and max, added n_layers and tr_regions[].
4-Mar-02  cwb -- moved this function after ppt enters soil.  Originally,
it was the first function called, so evapotransp was
based on swp_avg prior to wetting.  Also, set return
value as a pointer argument to be more consistent with
the rest of the code.
1-Oct-03  cwb -- Removed sum_tr_coeff[] requirement as it might as well
be calculated here and save the confusion of
having to keep up with it in the rest of the code.

INPUTS:
n_tr_lyrs  - number of layer regions used in weighted average
(typically 3, to represent shallow, mid, & deep depths)
to compute transpiration rate.
n_layers - number of soil layers
tr_regions - list of n_tr_lyrs elements of transp regions each
soil layer belongs to.
tr_coeff  - transpiration coefficients per layer.
sum_tr_coeff - sum of transpiration coefficients per layer structure.
swc -- current swc per layer

1-Oct-03 - local sumco replaces previous sum_tr_coeff[]

OUTPUT:
swp_avg - weighted average of soilwater potential and transpiration
coefficients
**********************************************************************/
	unsigned int r, i;
	double swp, sumco;
	
	*swp_avg = 0;
	for( r=1; r<= n_tr_rgns; r++) {
		swp = sumco = 0.0;
		
		for(i=0; i<n_layers; i++) {
			if (tr_regions[i] == r) {
				swp += tr_coeff[i] * swpotentl(swc[i], i);
				sumco += tr_coeff[i];
			}
		}
		
		swp /= GT(sumco, 0.) ? sumco : 1.;
		
		/* use smallest weighted average of regions */
		(*swp_avg) = (r==1) ? swp : fmin( swp, (*swp_avg));
	
	}
}

void grass_EsT_partitioning (double *fbse, double *fbst, double blivelai, double lai_param) {
/**********************************************************************
PURPOSE: Calculate fraction of water loss from bare soil
evaporation and transpiration

HISTORY:
4/30/92  (SLC)
24-Oct-03 (cwb) changed exp(-blivelai*bsepar1) + bsepar2;
to exp(-blivelai);

INPUTS:
blivelai - live biomass leaf area index

OUTPUTS:
fbse - fraction of water loss from bare soil evaporation.
fbst - "                           " transpiration.

**********************************************************************/
/* CWB- 4/00 Not sure what's the purpose of bsepar2, unless it's a
* fudge-factor to be played with.
*/

	double bsemax  = 0.995;
	
	*fbse = exp(- lai_param * blivelai);
	
	*fbse = fmin(*fbse, bsemax);
	*fbst = 1.-(*fbse);
}

void shrub_EsT_partitioning (double *fbse, double *fbst, double blivelai, double lai_param) {
/**********************************************************************
PURPOSE: Calculate fraction of water loss from bare soil
evaporation and transpiration
**********************************************************************/

	double bsemax  = 0.995;
	
	*fbse = exp(-lai_param * blivelai);
	
	*fbse = fmin(*fbse, bsemax);
	*fbst = 1.-(*fbse);
}

void tree_EsT_partitioning (double *fbse, double *fbst, double blivelai, double lai_param) {
/**********************************************************************
PURPOSE: Calculate fraction of water loss from bare soil
evaporation and transpiration

08/22/2011	(drs)	According to a regression based on a review by Daikoku, K., S. Hattori, A. Deguchi, Y. Aoki, M. Miyashita, K. Matsumoto, J. Akiyama, S. Iida, T. Toba, Y. Fujita, and T. Ohta. 2008. Influence of evaporation from the forest floor on evapotranspiration from the dry canopy. Hydrological Processes 22:4083-4096.
**********************************************************************/

	double bsemax  = 0.995;
	
	*fbse = exp(-lai_param * blivelai);
	
	*fbse = fmin(*fbse, bsemax);
	*fbst = 1.-(*fbse);
}


void pot_soil_evap( double *bserate, unsigned int nelyrs, double ecoeff[],
					double totagb, double fbse, double petday,
					double shift, double shape, double inflec,
					double range, double width[],  double swc[],
					double Es_param_limit)  {
/**********************************************************************
PURPOSE: Calculate potential bare soil evaporation rate.
See 2.11 in ELM doc.

HISTORY:
4/30/92  (SLC)
8/27/92  (SLC) Put in a check so that bserate cannot become
negative.  If total aboveground biomass (i.e.
litter+bimoass) is > 999., bserate=0.
6 Mar 02 (cwb) renamed watrate's parameters (see also SW_Site.h)
shift,  shift the x-value of the inflection point
shape,  slope of the line at the inflection point
inflec, y-value of the inflection point
range;  max y-val - min y-val at the limits
1-Oct-03 - cwb - removed the sumecoeff variable as it should
always be 1.0.  Also removed the line
avswp = sumswp / sumecoeff;

INPUTS:
nelyrs    - number of layers to consider in evaporation
sumecoeff - sum of evaporation coefficients
ecoeff    - array of evaporation coefficients
totagb    - sum of abovegraound biomass and litter
fbse      - fraction of water loss from bare soil evaporation
petday       - potential evapotranspiration rate
width     - array containing width of each layer.
swc  - array of soil water content per layer.

LOCAL:
avswp     - average soil water potential over all layers
evpar1    - input parameter to watrate.

OUTPUTS:
bserate   - bare soil evaporation loss rate. (cm/day)

FUNCTION CALLS:
watrate   - calculate evaporation rate.
swpotentl - compute soilwater potential
**********************************************************************/

	double  x,
	avswp = 0.0,
	sumwidth = 0.0;
	unsigned int    i;
	
	/* get the weighted average of swp in the evap layers */
	for (i=0; i < nelyrs; i++) {
		x = width[i] * ecoeff[i];
		sumwidth += x;
		avswp    += x * swpotentl( swc[i], i);
	}
	
	avswp /= sumwidth;
	
	
	/*  8/27/92 (SLC) if totagb > Es_param_limit, assume soil surface is
	* completely covered with litter and that bare soil
	* evaporation is inhibited.
	*/
	
	if ( GE(totagb,  Es_param_limit) ) {
		*bserate = 0.;
	} else {
		*bserate = petday * watrate(avswp, petday, shift, shape, inflec, range) * (1.-(totagb / Es_param_limit)) * fbse;
	}

}

void pot_transp(double *bstrate, double swpavg, double biolive,
				double biodead,  double fbst,   double petday,
				double swp_shift, double swp_shape, double swp_inflec, double swp_range,
				double shade_scale, double shade_deadmax, double shade_xinflex, double shade_slope, double shade_yinflex, double shade_range) {
/**********************************************************************
PURPOSE: Calculate potential transpiration rate.
See 2.11 in ELM doc.

HISTORY:
4/30/92  (SLC)
7/1/92   (SLC) Fixed bug.  Equation for bstrate not right.
9/1/92   (SLC) Allow transpiration, even if biodead is zero.  This
was in the original model, we will compute shade if biodead
is greater than deadmax.  Otherwise, shadeaf = 1.0
6 Mar 02 (cwb) renamed watrate's parameters (see also SW_Site.h)
shift,  shift the x-value of the inflection point
shape,  slope of the line at the inflection point
inflec, y-value of the inflection point
range;  max y-val - min y-val at the limits

INPUTS:
swpavg    - weighted average of soil water potential (from
function "transp_weighted_avg")
biolive   - biomass of live
biodead   - biomass of dead
fbst      - fraction of water loss from transpiration
petday       - potential evapotranspiration

LOCAL VARIABLES:
shadeaf - shade affect on transpiration rate
scale1  - scale for shade affect
trpar1  - input paramter to watrate
deadmax - maximum biomass of dead, before shade has any affect.

OUTPUTS:
bstrate   - transpiration loss rate. (cm/day)

FUNCTION CALLS:
watrate - compute transpiration rate.
tanfunc - tangent function
**********************************************************************/

	double  par1, par2, shadeaf;
	
	if ( LE(biolive, 0.) ) {
		*bstrate =0.;
		
		} else {
			if ( GE(biodead,  shade_deadmax) ) {
				par1 = tanfunc(biolive, shade_xinflex, shade_yinflex, shade_range, shade_slope);
				par2 = tanfunc(biodead, shade_xinflex, shade_yinflex, shade_range, shade_slope);
				shadeaf = (par1/par2) *(1.0-shade_scale) + shade_scale;
				shadeaf = fmin(shadeaf, 1.0);
			} else {
				shadeaf=1.0;
		}

	
		*bstrate = watrate(swpavg, petday, swp_shift, swp_shape, swp_inflec, swp_range) * shadeaf *  petday * fbst;
	}
}

double watrate(double swp, double petday, double shift, double shape,
double inflec, double range) {
/**********************************************************************
PURPOSE: Calculate the evaporation (or transpiration) rate, as
a function of potential evapotranspiration and soil
water potential. The ratio of evaporation (transpiration)
rate to PET is inversely proportional to soil water
poential (see Fig2.5a,b, pp.39, "Abiotic Section of ELM")

HISTORY:
4/30/92  (SLC)
6 Mar 02 (cwb) - Changed arguments from a/b/c to shift,shape,
inflec, range because I finally found the source
for tanfunc.

INPUTS:
swp - soil water potential (-bars)
petday - potential evapotranspiration rate for the day.
a   - equation parameter (relative to transpiration or evapor. rate)
b   - equation parameter (relative to transpiration or evapor. rate)
(usually b=.06 for evaporation, and b=.07 for transpiration)

OUTPUT:
watrate - rate of evaporation (or transpiration) from the
soil.

**********************************************************************/

	double  par1,par2, result;
	
	if ( LT(petday, .2) )
		par1 = 3.0;
	else if ( LT(petday, .4) )
		par1 = (.4 -petday) * -10. +5.;
	else if ( LT(petday, .6) )
		par1 = (.6 -petday) * -15. +8.;
	else
		par1 = 8.;
	
	par2 = shift - swp;
	
	result = tanfunc(par2, par1, inflec, range, shape);
	
	return ( fmin( fmax( result, 0.0), 1.0) );

}

// void evap_fromSurface( double *evap_pool, double *evap_demand, double *aet) {
// /**********************************************************************
// PURPOSE: Evaporate water from surface water pool, i.e., intercepted water (tree, shrub, grass, litter) or standingWater
// 	call seperately for each pool
// 
// INPUTS:
// evap_pool	- pool of surface water to evaporate from
// evap_demand	- yet unmet evaporative demand of atmosphere
// 
// OUTPUTS:
// evap_pool	- pool of surface water minus evaporated water
// evap_demand	- evaporative demand of atmosphere  minus evaporated water
// aet			- aet + evaporated water
// **********************************************************************/
// 
// 	double evap;
// 	double a=1., b=1.;
// 	
// 	evap = *evap_demand * a * (1. - exp(- b * (*evap_pool)));
// 	evap = fmax(0., fmin(*evap_demand, fmin(*evap_pool, evap )));
// 
// 	*evap_pool -= evap;
// 	*aet += evap;
// 	*evap_demand -= evap;	
// 
// 
// // 	if( GT(*evap_pool, *evap_demand) ) { /* demand is smaller than available water -> entire demand is evaporated */
// // 		*evap_pool -= *evap_demand;
// // 		*aet += *evap_demand;
// // 		*evap_demand = 0.;
// // 	} else { /* demand is larger than available water -> all available is evaporated */
// // 		*evap_demand -= *evap_pool;
// // 		*aet += *evap_pool;
// // 		*evap_pool = 0.;		
// // 	}
// }


void evap_fromSurface( double *water_pool, double *evap_rate, double *aet) {
/**********************************************************************
PURPOSE: Evaporate water from surface water pool, i.e., intercepted water (tree, shrub, grass, litter) or standingWater
	call seperately for each pool

INPUTS:
water_pool	- pool of surface water to evaporate from
evap_rate	- potential evaporation from this pool

OUTPUTS:
water_pool	- pool of surface water minus evaporated water
evap_rate	- actual evaporation from this pool
aet			- aet + evaporated water
**********************************************************************/

	if( GT(*water_pool, *evap_rate) ) { /* potential rate is smaller than available water -> entire potential is evaporated */
		*water_pool -= *evap_rate;
		*aet += *evap_rate;
	} else { /* potential rate is larger than available water -> entire pool is evaporated */
		*evap_rate = *water_pool;
		*aet += *water_pool;
		*water_pool = 0.;		
	}
}



void remove_from_soil( double swc[],   double qty[], double *aet,
unsigned int nlyrs,
double coeff[], double rate,  double swcmin[])  {
/**********************************************************************
PURPOSE: Remove water from the soil.  This replaces earlier versions'
call to separate functions for evaporation and transpiration
which did exactly the same thing but simply looked
different.  I prefer using one function over two to avoid
possible errors in the same transaction.

See Eqns. 2.12 - 2.18 in "Abiotic Section of ELM".

HISTORY: 10 Jan 2002 - cwb - replaced two previous functions with
this one.
12 Jan 2002 - cwb - added aet arg.
4 Dec 2002  - cwb - Adding STEPWAT code uncovered possible
div/0 error. If no transp coeffs, return.
INPUTS:
nlyrs  - number of layers considered in water removal
coeff - coefficients of removal for removal layers, either
evap_coeff[] or transp_coeff[].
rate - removal rate, either soil_evap_rate or soil_transp_rate.
swcmin - lower limit on soilwater content (per layer).

OUTPUTS:
swc  - soil water content adjusted after evaporation
qty - removal quantity from each layer, evap or transp.
aet - 

FUNCTION CALLS:
swpotentl - compute soilwater potential of the layer.
**********************************************************************/

	unsigned int  i;
	double	swpfrac[MAX_LAYERS],
			sumswp = 0.0,
			swc_avail, q;


	for (i=0; i < nlyrs; i++) {
		swpfrac[i] = coeff[i] / swpotentl(swc[i],i);
		sumswp += swpfrac[i];
	}

	if ( ZRO(sumswp) ) return;

	for (i=0; i < nlyrs; i++) {
		q = (swpfrac[i]/sumswp) * rate;
		swc_avail = fmax(0., swc[i] - swcmin[i]);
		qty[i]  = fmin( q, swc_avail);
		swc[i] -= qty[i];
		*aet   += qty[i];
	}
}

void infiltrate_water_low( double swc[],     double drain[],
double *drainout, unsigned int   nlyrs,
double sdrainpar, double sdraindpth,
double swcfc[],
double width[],   double swcmin[], double swcsat[],
double impermeability[], double *standingWater  ) {
/**********************************************************************
PURPOSE:Calculate soilwater drainage for low soil water conditions
See Equation 2.9 in ELM doc.

HISTORY:
4/30/92  (SLC)
7/2/92   (fixed bug.  equation for drainlw needed fixing)
8/13/92 (SLC) Changed call to function which checks lower bound
on soilwater content.  REplaced call to "chkzero" with
the function "getdiff".
- lower bound is used in place of zero as lower bound
here.  Old code used 0.cm water as a lower bound in
low water drainage.
9/22/01 - (cwb) replaced tr_reg_max[] with transp_rgn[]
see INPUTS
1/14/02 - (cwb) fixed off by one error in loop.
6-Oct-03  (cwb) removed condition disallowing gravitational
drainage from transpiration region 1.

INPUTS:
drain - drainage from each layer
ndeeplyr - bottom layer to stop at for drainage
transp_rgn - array of transpiration regions of each layer
sdrainpar - slow drainage parameter
swcfc   - soilwater content at field capacity
swcwp   - soilwater content at wilting point
swc  - soil water content adjusted by drainage
drain - drainage from each soil water layer

OUTPUTS:
swc  - soil water content adjusted by drainage
drain - drainage from each soil water layer
drainout - added low drainout (to already calculated high drainout)
**********************************************************************/

unsigned int i;
int j;
double drainlw = 0.0, swc_avail, drainpot, d[nlyrs], push;


	for (i=0; i < nlyrs; i++) {
		/* calculate potential unsaturated percolation */
		if (LE(swc[i], swcmin[i]) ) { /* in original code was !GT(swc[i], swcwp[i]) equivalent to LE(swc[i], swcwp[i]), but then water is drained to swcmin nevertheless - maybe should be LE(swc[i], swcmin[i]) */
			d[i] = 0.;
		} else {
			swc_avail = fmax(0., swc[i] - swcmin[i]);
			drainpot = GT(swc[i], swcfc[i])
				? sdrainpar
				: sdrainpar * exp( (swc[i] - swcfc[i]) * sdraindpth / width[i]);
			d[i] = (1. - impermeability[i]) * fmin(swc_avail, drainpot);
		}
		drain[i] += d[i];
		
		if(i < nlyrs-1){ /* percolate up to next-to-last layer */
			swc[i+1] += d[i];
			swc[i]   -= d[i];
		} else { /* percolate last layer */
			drainlw = fmax( d[i], 0.0);
			(*drainout) += drainlw;
			swc[i] -= drainlw;
		}
		
	}

/* adjust (i.e., push water upwards) if water content of a layer is now above saturated water content */
	for (j=nlyrs; j >= 0; j--) {
		if( GT(swc[j], swcsat[j]) ) {
			push = swc[j] - swcsat[j];
			swc[j] -= push;
			if(j > 0) {
				drain[j-1] -= push;
				swc[j-1] += push;
			} else {
				(*standingWater) += push;
			}
		}
	}

}


void hydraulic_redistribution(	double swc[], double swcwp[], double lyrRootCo[], double hydred [],
								unsigned int   nlyrs, double maxCondroot,
								double swp50, double shapeCond,
								double scale)
{
/**********************************************************************
PURPOSE:Calculate hydraulic redistribution according to Ryel, Ryel R, Caldwell, Caldwell M, Yoder, Yoder C, Or, Or D, Leffler, Leffler A. 2002. Hydraulic redistribution in a stand of Artemisia tridentata: evaluation of benefits to transpiration assessed with a simulation model. Oecologia 130: 173-184.

HISTORY:
	10/19/2010 (drs)
	11/13/2010 (drs) limited water extraction for hydred to swp above wilting point
	03/23/2012 (drs) excluded hydraulic redistribution from top soil layer (assuming that this layer is <= 5 cm deep)
	
INPUTS:
swc  - soil water content
lyrRootCo - fraction of active roots in layer i
nlyrs  - number of soil layers
maxCondroot - maximum radial soil-root conductance of the entire active root system for water (cm/-bar/day)
swp50 - soil water potential (-bar) where conductance is reduced by 50%
shapeCond - shaping parameter for the empirical relationship from van Genuchten to model relative soil-root conductance for water
scale - fraction of vegetation type to scale hydred

OUTPUTS:
swc  - soil water content adjusted by hydraulic redistribution
hydred - hydraulic redistribtion for each soil water layer (cm/day/layer)

**********************************************************************/

	unsigned int i, j;
	double swp[nlyrs], swpwp[nlyrs], relCondroot[nlyrs], hydredmat[nlyrs][nlyrs], Rx, swa, hydred_sum;

	hydred[0] = 0.;	/* no hydred in top layer */

	for( i=0; i < nlyrs; i++) {
		swp[i] = swpotentl(swc[i],i);
		relCondroot[i] = fmin( 1., fmax(0., 1./(1. + pow(swp[i]/swp50, shapeCond) ) ) );
		swpwp[i] = swpotentl(swcwp[i],i);
		
		hydredmat[0][i] = hydredmat[i][0] = 0.;	/* no hydred in top layer */
	}	

	for( i=1; i< nlyrs; i++ ) {
		hydred[i] = hydredmat[i][i] = 0.;	/* no hydred within any layer */
		
		for( j=i+1; j< nlyrs; j++ ) {
		
			if( LT(swp[i], swpwp[i]) || LT(swp[j], swpwp[j]) ){ /* hydred occurs only if source layer's swp is above wilting point */
			
				if( GT(swc[i], swc[j]) ) {Rx = lyrRootCo[i];} else {Rx = lyrRootCo[j];}
			
				hydredmat[i][j] = maxCondroot * 10./24. * (swp[j] - swp[i]) * fmax(relCondroot[i], relCondroot[j])
								* (lyrRootCo[i] * lyrRootCo[j] / (1. -  Rx)); /* assuming a 10-hour night */
				hydredmat[j][i] = -hydredmat[i][j];
			} else {
				hydredmat[i][j] = hydredmat[j][i] = 0.;
			}
		}
	}

	for( i=0; i < nlyrs; i++) { /* total hydred from layer i cannot extract more than its swa */
		hydred_sum = 0.;
		for( j=0; j< nlyrs; j++ ) {
			hydred_sum += hydredmat[i][j];		
		}
		
		swa = fmax( 0., swc[i] - swcwp[i] );
		if( LT(hydred_sum, 0.) && GT( -hydred_sum, swa) ) {
			for( j=0; j< nlyrs; j++ ) {
				hydredmat[i][j] *= (swa / -hydred_sum);
				hydredmat[j][i] *= (swa / -hydred_sum);
			}
		}
	}
	
	for( i=0; i < nlyrs; i++) {
		for( j=0; j< nlyrs; j++ ) {
			hydred[i] += hydredmat[i][j] * scale;
		}
		swc[i] += hydred[i];
	}



}

/**********************************************************************
PURPOSE: Initialize soil temperature regression values, only needs to be called once (ie the first time the soil_temperature function is called).  this is not included in the header file since it is NOT an external function

HISTORY: 
	05/31/2012 (DLM) initial coding

INPUTS: they are all defined in the soil_temperature function, so go look there to see their meaning as it would be redundant to explain them all here as well.

OUTPUT: none, but places the regression values in stValues struct for use in the soil_temperature function later
**********************************************************************/
void soil_temperature_init( double bDensity[], double width[],
							double oldsTemp[], unsigned int nlyrs, double fc[],
							double wp[], double deltaX, double theMaxDepth,
							double meanAirTemp, unsigned int nRgr) 
{	
	// local vars
	unsigned int i, j, k, x1 = 1, x2 = 1, equal = 1;
	double acc = 0.0;
	// pointers
	unsigned int *x1P = &x1, *x2P = &x2, *equalP = &equal;
	ST_RGR_VALUES *st = &stValues; // just for convenience, so I don't have to type as much
	
	soil_temp_init = 1; // make this value 1 to make sure that this function isn't called more than once... (b/c it doesn't need to be)

	for( i=0; i < nlyrs; i++) {
		acc += width[i];
		st->depths[i] = acc;
	}
	
	for( i=0; i < nRgr + 1; i++)
		*(st->depthsR + i) = ( (deltaX * i) + deltaX); 	
		
	// if there's less then 2 lyrs of soil, or the max layer depth < 30 cm the function quits (& prints out an error message) so it doesn't blow up later...
	if( (nlyrs < 2) || LT(st->depths[nlyrs - 1], deltaX + deltaX) ) { 
		if(!soil_temp_error) { // if the error hasn't been reported yet... print an error to the stderr and one to the logfile
			// fprintf(stderr, "\nSOIL_TEMP FUNCTION ERROR: (there needs to be >= 2 soil layers, with a maximum combined depth of >= %5.2f cm)... soil temperature will NOT be calculated\n", (deltaX + deltaX));
			fprintf(logfp, "\nSOIL_TEMP FUNCTION ERROR: (there needs to be >= 2 soil layers, with a maximum combined depth of >= %5.2f cm)... soil temperature will NOT be calculated\n", (deltaX + deltaX));
			soil_temp_error = 1;
		}
		return; // exits the function
	}
	
	acc = deltaX;
	k = j = 0;
	// linear regression time complexity of this should be something like O(k * nlyrs).  might be able to do it faster... but would also be a lot more code & would require a rewrite (shouldn't matter anyways, because this function is only called once)
	while(LE(acc, st->depths[nlyrs - 1])) {
		st_getBounds(x1P, x2P, equalP, nlyrs, acc, st->depths);
		
		i = -1;
		if(x1 == i) { // sets the values to the first layer of soils values, since theres nothing else that can be done... fc * wp must be scaled appropriately
			st->fcR[k] = (fc[0] / width[0]) * deltaX; // all the division and multiplication is to make sure that the regressions are scaled appropriately, since the widths of the layers & the layers of the regressions, may not be the same
			st->wpR[k] = (wp[0] / width[0]) * deltaX;
			// st->oldsTempR[k] = regression(0.0, st->depths[0], T1, oldsTemp[0], acc); // regression using the temp at the top of the soil, commented out b/c it's giving worse results
			st->oldsTempR[k] = oldsTemp[0]; // no scaling is necessary with temperature & bulk density
			st->bDensityR[k] = bDensity[0];
		} else if( (x1 == x2) || (equal != 0) ) { // sets the values to the layers values, since x1 and x2 are the same, no regression is necessary
			st->fcR[k] = (fc[x1] / width[x1]) * deltaX;
			st->wpR[k] = (wp[x1] / width[x1]) * deltaX;
			st->oldsTempR[k] = oldsTemp[x1];
			st->bDensityR[k] = bDensity[x1];
		} else { // double regression( double x1, double x2, double y1, double y2, double deltaX ), located in generic.c
			st->fcR[k] = regression(st->depths[x1], st->depths[x2], (fc[x1] / width[x1]) * deltaX, (fc[x2] / width[x2]) * deltaX, acc);
			st->wpR[k] = regression(st->depths[x1], st->depths[x2], (wp[x1] / width[x1]) * deltaX, (wp[x2] / width[x2]) * deltaX, acc);
			st->oldsTempR[k] = regression(st->depths[x1], st->depths[x2], oldsTemp[x1], oldsTemp[x2], acc);
			st->bDensityR[k] = regression(st->depths[x1], st->depths[x2], bDensity[x1], bDensity[x2], acc);
		}
		
		if(equal != 0)
			x2 = x1;
		st->x1BoundsR[k] = x1;
		st->x2BoundsR[k] = x2;
		
		k++;
		acc += deltaX;
	}
	
	// this next chunk is commented out... the code is kept here though if if we want to change back to extrapolating the rest of the regression values
	
	// to fill the rest of the regression values, simply use the last two values since there is no more actual data to go off of
	// if k is < 2 this code will blow up... but that should never happen since the function quits if there isn't enough soil data earlier
	/*for( i=k; i < nRgr; i++) {
		st->wpR[i] = regression(st->depthsR[i - 2], st->depthsR[i - 1], st->wpR[i - 2], st->wpR[i - 1], st->depthsR[i]);
		st->fcR[i] = regression(st->depthsR[i - 2], st->depthsR[i - 1], st->fcR[i - 2], st->fcR[i - 1], st->depthsR[i]);
		st->bDensityR[i] = regression(st->depthsR[i - 2], st->depthsR[i - 1], st->bDensityR[i - 2], st->bDensityR[i - 1], st->depthsR[i]);
	}
	
	// getting the average for a regression...
	for( i=0; i < nlyrs; i++) {
		wpAverage += wp[i] / width[i];
		fcAverage += fc[i] / width[i];
		bDensityAverage += bDensity[i];
	}
	wpAverage = deltaX * (wpAverage / (nlyrs + 0.0));
	fcAverage = deltaX * (fcAverage / (nlyrs + 0.0));
	bDensityAverage = (bDensityAverage / (nlyrs + 0.0));
	
	// if the values are too small, we reset them to the average for the regression... it's a safeguard
	for( i=k; i < nRgr; i++ ) {
		if(LT(st->wpR[i], 1))
			st->wpR[i] = wpAverage;
		if(LT(st->fcR[i], 2))
			st->fcR[i] = fcAverage;
		st->bDensityR[i] = bDensityAverage; // just set the bulk density to the average, it seems to be a better approximation...
	}*/
	
	// just use the last soil layers values...
	for( i=k; i < nRgr; i++ ) {
		st->wpR[i] = deltaX * (wp[nlyrs - 1] / width[nlyrs - 1]);
		st->fcR[i] = deltaX * (fc[nlyrs - 1] / width[nlyrs - 1]);
		st->bDensityR[i] = bDensity[nlyrs - 1];
	}
	
	if(k < nRgr) //was k < 11
		st->oldsTempR[k] = regression(st->depths[nlyrs - 1], theMaxDepth, oldsTemp[nlyrs - 1], meanAirTemp, st->depthsR[k]); // to give a slightly better temp approximation
	for( i=k + 1; i < nRgr; i++)
		st->oldsTempR[i] = regression(st->depthsR[i - 1], theMaxDepth, st->oldsTempR[i - 1], meanAirTemp, st->depthsR[i]); // we do temperature differently, since we already have the temperature for the last layer of soil
	
	st->oldsTempR[nRgr] = meanAirTemp; // the soil temp at the last layer of the regression is equal to the meanAirTemp, this is constant so it's the same for yesterdays temp & todays temp

	// getting all the xBounds values for later use in the soil_temperature function...
	for( i=0; i < nlyrs; i++ ) {
		st_getBounds(x1P, x2P, equalP, nRgr + 1, st->depths[i], st->depthsR);
		if(equal != 0)
			x2 = x1;
		st->x1Bounds[i] = x1;
		st->x2Bounds[i] = x2;
	}
}

/**********************************************************************
PURPOSE: Calculate soil temperature for each layer as described in Parton 1978, ch. 2.2.2 Temperature-profile Submodel, regression values are gotten from a mixture of interpolation & extrapolation 

*NOTE* There will be some degree of error because the original equation is written for soil layers of 15 cm.  if soil layers aren't all 15 cm then linear regressions are used to estimate the values (error should be relatively small though).
*NOTE* Function might not work correctly if the maxDepth of the soil is > 180 cm, since Parton's equation goes only to 180 cm
*NOTE* Function will run if maxLyrDepth > maxDepth of the equation, but the results might be slightly off...

HISTORY: 
	05/24/2012 (DLM) initial coding, still need to add to header file, handle if the layer height is > 15 cm properly, & test
	05/25/2012 (DLM) added all this 'fun' crazy linear regression stuff
	05/29/2012 (DLM) still working on the function, linear regression stuff should work now.  needs testing
	05/30/2012 (DLM) got rid of nasty segmentation fault error, also tested math seems correct after checking by hand.  added the ability to change the value of deltaX
	05/30/2012 (DLM) added # of lyrs check & maxdepth check at the beginning to make sure code doesn't blow up... if there isn't enough lyrs (ie < 2) or the maxdepth is too little (ie < deltaX * 2), the function quits out and reports an error to the user
	05/31/2012 (DLM) added theMaxDepth variable to allow the changing of the maxdepth of the equation, also now stores most regression data in a structure to reduce redundant regression calculations & speeds things up
	06/01/2012 (DLM) changed deltaT variable from hours to seconds, also changed some of the regression calculations so that swc, fc, & wp regressions are scaled properly... results are actually starting to look usable!
	06/13/2012 (DLM) no longer extrapolating values for regression layers that are out of the bounds of the soil layers... instead they are now set to the last soil layers values.  extrapolating code is still in the function and can be commented out and used if wishing to go back to extrapolating the values...

INPUTS:
	airTemp - the average daily air temperature in celsius
	pet - the potential evapotranspiration rate
	aet - the actual evapotranspiration rate
	biomass - the standing-crop biomass 
	swc - soil water content
	bDensity - bulk density of the soil layers
	width - width of layers
	oldsTemp - soil layer temperatures from the previous day in celsius
	nlyrs - number of soil layers, must be greater than 1 or the function won't work right
	fc - field capacity for each soil layer
	wp - wilting point for each soil layer
	bmLimiter - biomass limiter constant (300 g/m^2)
	t1Params - constants for the avg temp at the top of soil equation (15, -4, 600) there is 3 of them
	csParams - constants for the soil thermal conductivity equation (0.00070, 0.00030) there is 2 of them
	shParam - constant for the specific heat capacity equation (0.18)
	snowpack - the amount of snow on the ground
	meanAirTemp - the avg air temperature for the month in celsius
	deltaX - the distance between profile points (default is 15 from Parton's equation, wouldn't recommend changing the value from that).  180 must be evenly divisible by this number.
	theMaxDepth - the lower bound of the equation (default is 180 from Parton's equation, wouldn't recommend changing the value from that).
	nRgr - the number of regressions (1 extra value is needed for the sTempR and oldsTempR for the last layer

OUTPUT:
	sTemp - soil layer temperatures in celsius
**********************************************************************/

void soil_temperature(	double airTemp, double pet, double aet, double biomass,  
						double swc[], double bDensity[], double width[],
						double oldsTemp[], double sTemp[], unsigned int nlyrs, 
						double fc[], double wp[], double bmLimiter, 
						double t1Param1, double t1Param2, double t1Param3, 
						double csParam1, double csParam2, double shParam,
						double snowpack, double meanAirTemp, double deltaX,
						double theMaxDepth, unsigned int nRgr )
{
	unsigned int i, j, k, x1 = 1, x2 = 1;
	double T1, cs, sh, sm, pe, deltaT, part1, part2, acc, swcR[nRgr], sTempR[nRgr + 1], maxLyrDepth;
	int toDebugSoilTemp = 0;
	/* local variables explained: 
	
	   	toDebugSoilTemp - 1 to print out debug messages & then exit the program after completing the function, 0 to not.  default is 0.
	   	T1 - the average daily temperature at the top of the soil in celsius
	   	sm - volumetric soil-water content
	   	pe - ratio of the difference between volumetric soil-water content & soil-water content
	    	at the wilting point to the difference between soil water content at field capacity &
	    	soil-water content at wilting point.
	   	cs - soil thermal conductivity
	   	sh - specific heat capacity  
	   	deltaT - time step (24 hr)
	   	depths[nlyrs] - the depths of each layer of soil, calculated in the function
	   	swcR[], sTempR[] - anything with a R at the end of the variable name stands for the regression of that array
	*/
	
	ST_RGR_VALUES *st = &stValues; // just for convenience, so I don't have to type as much
	
	deltaT = 86400.0; // the # of seconds in a day... (24 hrs * 60 mins/hr * 60 sec/min = 86400 seconds)
	
	// calculating T1, the average daily air temperature at the top of the soil
	if( LE(biomass, bmLimiter)) { // bmLimiter = 300
		T1 = airTemp + (t1Param1 * pet * (1. - ((aet/pet) * (1. - (biomass/bmLimiter))) ) ); // t1Param1 = 15; math is correct
		if(toDebugSoilTemp == 1) printf("\nT1 = %5.4f + (%5.4f * %5.4f * (1 - ((%5.4f / %5.4f) * (1 - (%5.4f / %5.4f))) ) )", airTemp, t1Param1, pet, aet, pet, biomass, bmLimiter);
	} else {
		T1 = airTemp + ((t1Param2 * (biomass - bmLimiter)) / t1Param3); // t1Param2 = -4, t1Param3 = 600; math is correct
		if(toDebugSoilTemp == 1) printf("\nT1 = %5.4f + ((%5.4f * (%5.4f - %5.4f)) / %5.4f)", airTemp, t1Param2, biomass, bmLimiter, t1Param3);
	}
	
	if(toDebugSoilTemp == 1) printf("\nAirTemp : %5.4f pet : %5.4f aet : %5.4f biomass : %5.4f bmLimiter : %5.4f", airTemp, pet, aet, biomass, bmLimiter);
	
	if(GT(snowpack, 0.0)) {// if there is snow on the ground, then T1 is simply set to -2
		T1 = -2.0;
		if(toDebugSoilTemp == 1) printf("\nThere is snow on the ground, T1 set to -2\n");
	}
	
	if(!soil_temp_init)
		soil_temperature_init(bDensity, width, oldsTemp, nlyrs, fc, wp, deltaX, theMaxDepth, meanAirTemp, nRgr);	
		
	if(soil_temp_error) // if there is an error found in the soil_temperature_init function, return so that the function doesn't blow up later
		return;

	maxLyrDepth = st->depths[nlyrs - 1];
	
	k = 0; // k keeps track of which layer of the regression we are on...
	acc = deltaX;
	if(toDebugSoilTemp == 1) printf("\nT1 : %5.4f meanAirTemp : %5.4f \nnlyrs : %d maxLyrDepth : %5.4f \n \n", T1, meanAirTemp, nlyrs, maxLyrDepth);
	// linear regression
	while(LE(acc, maxLyrDepth)) {
	
		x1 = st->x1BoundsR[k];
		x2 = st->x2BoundsR[k];
		i = -1;
		if(toDebugSoilTemp == 1) {
			if(x1 != i) { // makes sure that we're not sending st->depths[-1] to printf, b/c as it turns out that causes a nasty segmentation fault error...
				printf("k %d %d - %d depthLow %5.4f acc %5.4f depthHigh %5.4f\n", k, x1, x2, st->depths[x1], acc, st->depths[x2]);
			} else {
				printf("k %d %d - %d depthLow %5.4f acc %5.4f depthHigh %5.4f\n", k, x1, x2, 0.0, acc, st->depths[x2]);
			}
		}
		
		if(x1 == i) { // sets the values to the first layer of soils values, since theres nothing else that can be done...
			swcR[k] = (swc[0] / width[0]) * deltaX; // division & multiplication is to make sure that the values are scaled appropriately, since the width of the layer & the width of a layer of the regression may not be the same
		} else if(x1 == x2) { // sets the values to the layers values, since x1 and x2 are the same, no regression is necessary (scaling still is necessary, however)
			swcR[k] = (swc[x1] / width[x1]) * deltaX;
		} else { // double regression( double x1, double x2, double y1, double y2, double deltaX ), located in generic.c
			part1 = (x2 != i) ? ((swc[x2] / width[x2]) * deltaX) : ((swc[nlyrs - 1] / width[nlyrs - 1]) * deltaX);
			swcR[k] = regression(st->depths[x1], (x2 != i) ? st->depths[x2] : st->depths[nlyrs - 1], (swc[x1] / width[x1]) * deltaX, part1, acc);
		}
		
		k++;
		acc += deltaX;
	}
	
	// uncomment out this next part if wanting to change back to extrapolating the rest of the regression values...
	 	
	// to fill the rest of the regression values, simply use the last two values since there is no more actual data to go off of
	/*for( i=k; i < nRgr; i++)
		swcR[i] = regression(st->depthsR[i - 2], st->depthsR[i - 1], swcR[i - 2], swcR[i - 1], st->depthsR[i]);
		
	
	// getting the average for a regression...
	for( i=0; i < k; i++) {
		swcAverage += swcR[i];
	}
	swcAverage = swcAverage / (k + 0.0);
		
	// resets the regression values if they're too small... it's a safeguard...
	for( i=0; i < nRgr; i++ )
		if(LT(swcR[i], 1))
			swcR[i] = swcAverage;*/
			
	//just use the last layers values...		
	for( i=k; i < nRgr; i++ ) 
		swcR[i] = deltaX * (swc[nlyrs - 1] / width[nlyrs - 1]);
	
	if(toDebugSoilTemp == 1) printf("\nregression values: \n");
	if(toDebugSoilTemp == 1) for( i=0; i < nRgr; i++)
					 printf("k %d swcR %5.4f fcR %5.4f wpR %5.4f oldsTempR %5.4f bDensityR %5.4f \n", i, swcR[i], st->fcR[i], st->wpR[i], st->oldsTempR[i], st->bDensityR[i]);
	if(toDebugSoilTemp == 1) printf("\nlayer values: \n");
	if(toDebugSoilTemp == 1) for( i=0; i < nlyrs; i++)
					 printf("i %d width %5.4f depth %5.4f swc %5.4f fc %5.4f wp %5.4f oldsTemp %5.4f bDensity %5.4f \n", i, width[i], st->depths[i], swc[i], fc[i], wp[i], oldsTemp[i], bDensity[i]);
	
	// FINALLY done with the regressions!!! this is where we calculate the temperature for each soil layer of the regression
	if(toDebugSoilTemp == 1) printf("\n");
	for( i=0; i < nRgr; i++) { // goes to nRgr, because the soil temp of the last regression layer (nRgr) is the meanAirTemp
	
		// first we must calculate cs & sh (& subsequently sm & pe), for use later
		sm = swcR[i];
		pe = (sm - st->wpR[i]) / (st->fcR[i] - st->wpR[i]);
		cs = csParam1 + (pe * csParam2); // csParam1 = 0.0007, csParam2 = 0.0003
		sh = sm + (shParam * (1. - sm)); // shParam = 0.18
		
		if(toDebugSoilTemp == 1) printf("k %d cs %5.4f sh %5.4f\n", i, cs, sh); 
		
		// breaking the equation down into parts to make it easier for me to process
		part1 = cs / (sh * st->bDensityR[i]); 
		
		if(i > 0) { // handles all layers except the first soil layer
			part2 = (sTempR[i - 1] - (2 * st->oldsTempR[i]) + st->oldsTempR[i + 1]) / squared(deltaX);
		} else { // handles the first soil layer, since it needs the temp of the top of the soil
			part2 = (T1 - (2 * st->oldsTempR[0]) + st->oldsTempR[1]) / squared(deltaX);
		}
		
		sTempR[i] = ((part1 * part2) * deltaT) + st->oldsTempR[i];
	}
	sTempR[nRgr] = meanAirTemp; // again... the last layer of the regression is set to the constant meanAirTemp
	
	// MORE REGRESSIONS! to change sTempR into sTemp for outputting correctly
	if(toDebugSoilTemp == 1) printf("\n");
	j = -1; // have to do this to avoid signed / unsigned comparison warning.  it's a pain really but have to do since the program was all written using unsigned ints (which can be no negative value except for -1) for some reason.  keep in mind that -1 is the largest possible value that an unsigned int can be according to c (ie -1 > 100 is a true statement if they're both unsigned ints b/c -1 is converted to UINT_MAX), very confusing.
	for( i=0; i < nlyrs; i++ ) {
		x1 = st->x1Bounds[i];
		x2 = st->x2Bounds[i];
		
		if(toDebugSoilTemp == 1) printf("i %d %d - %d depthLow %5.4f acc %5.4f depthHigh %5.4f\n", i, x1, x2, ( (x1 + 1) * deltaX), st->depths[i], (x2 != j) ? st->depthsR[x2] : -1.0);
		
		if(x1 == j) { // makes a regression with the temp at the top of the soil & the first soil layer...
			//sTemp[i] = regression(0.0, deltaX, T1, sTempR[0], st->depths[i]); // commented out, b/c it was actually giving a worse approximation
			sTemp[i] = regression(deltaX, deltaX + deltaX, sTempR[0], sTempR[1], st->depths[i]);
		} else if(x1 == x2) { // sets the values to the layers values, since x1 and x2 are the same, no regression is necessary
			sTemp[i] = sTempR[x1];
		} else { 
			if(x2 != j)
				sTemp[i] = regression( st->depthsR[x1], st->depthsR[x2], sTempR[x1], sTempR[x2], st->depths[i]);
			else
				sTemp[i] = regression( st->depthsR[x1], (x2 != j) ? st->depthsR[x2] : st->depthsR[nRgr - 1], sTempR[x1], (x2 != j) ? sTempR[x2] : sTempR[nRgr - 1], st->depths[i]);
		}
	}
	
	if(toDebugSoilTemp == 1) printf("\nregression temp values: \n");
	if(toDebugSoilTemp == 1) for( i=0; i < nRgr + 1; i++)
		 			printf("k %d oldsTempR %5.4f sTempR %5.4f depth %5.4f \n", i, *(st->oldsTempR + i), sTempR[i], ( (i + 1) * deltaX) ); // *(oldsTempR + i) is equivalent to writing oldsTempR[i]
		
	if(toDebugSoilTemp == 1) printf("\nlayer temp values: \n");
	if(toDebugSoilTemp == 1) for( i=0; i < nlyrs; i++)
					printf("i %d oldTemp %5.4f sTemp %5.4f depth %5.4f  \n", i, oldsTemp[i], sTemp[i], st->depths[i]);
					
	// updating the values of yesterdays temperature regression for the next time the function is called...
	for( i=0; i < nRgr + 1; i++)
		st->oldsTempR[i] = sTempR[i];
		
	if(toDebugSoilTemp == 1) exit(0); // terminates the program, make sure to take this out later
}