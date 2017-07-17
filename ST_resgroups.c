/********************************************************/
/********************************************************/
/*  Source file: resgroups.c
 /*  Type: module
 /*  Application: STEPPE - plant community dynamics simulator
 /*  Purpose: This is the module that executes the growth
 *           functions and otherwise manages the bookkeeping
 *           at the group level.
 /*  History:
 /*     (6/15/2000) -- INITIAL CODING - cwb
 /*
 /********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "myMemory.h"
#include "rands.h"
#include "generic.h"

#include "ST_functions.h"
#include "sxw_funcs.h"
extern Bool UseSoilwat;

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
void species_Update_Estabs(SppIndex sp, IntS num);

/*------------------------------------------------------*/
/* Modular functions only used only one or two specific */
/* places; that is, they are not generally useful       */
/* anyplace in the model, but they have to be declared. */
void rgroup_Grow(void);
void rgroup_Establish(void);
void rgroup_IncrAges(void);
void rgroup_PartResources(void);
void rgroup_ResPartIndiv(void);
void rgroup_DropSpecies(SppIndex sp);
void rgroup_AddSpecies(GrpIndex rg, SppIndex sp);
void rgroup_Extirpate(GrpIndex rg);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static RealF _ppt2resource(RealF ppt, GroupType *g);
static void _res_part_extra(Bool isextra, RealF extra, RealF size[]);
static GroupType *_create(void);
static void _extra_growth(GrpIndex rg);
static void _add_annual_seedprod(SppIndex sp, RealF pr);
static RealF _get_annual_maxestab(SppIndex sp);
static RealF _add_annuals(const GrpIndex rg, const SppIndex sp, const RealF lastyear_relsize);

/****************** Begin Function Code ********************/
/***********************************************************/

void rgroup_PartResources(void)
{
	/*======================================================*/
	/* PURPOSE */
	/* Partition resources for this year among the resource
	 groups.  The allocation happens in three steps: basic
	 group allocation (which is done here), partitioning of
	 extra resources (in _res_part_extra()), and further
	 partitioning to individuals in the group (_ResPartIndiv()).

	 See COMMENT 1. at the end of this file for the algorithm.

	 For future documenting reference, note that the one
	 requirement regarding slope/intercept is that the
	 Globals.ppt.avg * slope + intercept == 1.0 for all
	 groups.
	 */
	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 12/15/2000            */
	/* This code was transformed from an earlier version's
	 Env_Partition() and moved here because I'd like to
	 rewrite the model as a C++ program and the RGroups
	 will become a class.

	 8-Apr-2003 - cwb -- Well, after ALL this time, it finally
	 dawned on me to put the resource calculations
	 into a spreadsheet and figure them out correctly.
	 Affected lines are g->res_avail, and remains.

	 15-May-03 (cwb) Bill suggested partitioning resource by
	 mm instead of proportionally, however, the results
	 aren't precisely the same as the proportional
	 method, so I'm making code for both methods.
	 /*------------------------------------------------------*/

	GrpIndex rg;
	SppIndex sp;
	RealF resource, /* amt of "resource" == 1 when ppt is avg */
	xtra_base = 0., /* pooled extra resource up to 1.0 */
	xtra_obase = 0., /* pooled resource > 1.0 */
	size_base[MAX_RGROUPS], /* total res. contrib to base, all groups */
	size_obase[MAX_RGROUPS]; /* total res. contrib. if xtra_obase */

	Bool noplants = TRUE;
	const Bool do_base = FALSE, /* monikers for _res_part_extra() */
	do_extra = TRUE;
	GroupType *g; /* shorthand for RGroup[rg] */
	int i;

	/*----------------------------------------------------*/

	/* ----- distribute basic (minimum) resources */
	ForEachGroup(rg)
	{
		g = RGroup[rg];

		//if (g->max_age == 1)
			//g->relsize = _add_annuals(rg, 1.0, no_seeds);

		/*this piece of the code is only used when SOILWAT is NOT running*/
#ifdef STEPWAT
		if (!UseSoilwat)
		{ /* use by-mm method */
#endif
		resource = _ppt2resource(Env.ppt, g);
		g->res_required = g->relsize / g->max_density;
		g->res_avail = fmin(1., fmin(g->res_required, resource));
		xtra_base += fmax(0., fmin(1., resource) - g->res_avail)
				* g->min_res_req;
		xtra_obase += fmax( 0., resource - 1.) * g->min_res_req;

		size_base[rg] = g->relsize * g->min_res_req;
		size_obase[rg] = (g->use_extra_res) ? size_base[rg] : 0.;

		/*this is how res_required and res_avail are set if SOILWAT is running*/
#ifdef STEPWAT
	}
	else
	{
		/* trick for later pr calc */

		g->res_required = RGroup_GetBiomass(rg);
		g->res_avail = SXW_GetTranspiration(rg); //1.0;
		//PR limit can be really high see mort no_resource I limit to the groups estab indiv because that is what we can kill.
		//This was at 10 but ten might be to low or high in some cases.
		if(!ZRO(g->res_avail) && g->res_required / g->res_avail > g->estabs)
		{
			g->res_required = g->estabs;
			g->res_avail = 1;
		}
		//A check
		if(ZRO(g->res_avail) && g->res_required > 0)
		{
			g->res_required = g->estabs;
			g->res_avail = 1;
			LogError(logfp, LOGWARN, "RGroup %s : res_avail is Zero and res_required > 0", g->name);
		}
        //KAP:Previously, a seperate set of code was used to determine resource availability for annual species
        //Now, resource paritioning for annuals occurs in the same way as for perennials, by getting the RGroup Biomass
        //and the resource_cur from SXW_GetTranspiration
		//Annuals seem to have a artificial limit of 20. We do Annuals here differently
		if(g->max_age == 1)
		{
			//{	ForEachGroupSpp(sp,rg,i)
			//	g->res_required += Species[sp]->mature_biomass * .75;}
			//g->res_avail = SXW_GetTranspiration(rg);
			if(!ZRO(g->res_avail) && g->res_required / g->res_avail > 20)
			{
				g->res_required = 20;
				g->res_avail = 1;
			}
			if(ZRO(g->res_avail) && g->res_required > 0)
			{
				g->res_required = 20;
				g->res_avail = 1;
			}
		}
	}
#endif
	  /* If relsize>0, reset noplants from TRUE to FALSE and if noplants=TRUE, exit from the loop */
		if (GT(g->relsize, 0.))
			noplants = FALSE;

	} /* End ForEachGroup(rg) */

	if (noplants)
		return;

	/*these functions are not used if using SOILWAT,
	 extra resource partitioning does not occur when running SOILWAT*/
		_res_part_extra(do_base, xtra_base, size_base);
		_res_part_extra(do_extra, xtra_obase, size_obase);
	
	/* reset annuals' "true" relative size here */
    //KAP: formely, this call established annual species. We have moved annual establishment to the Rgroup_Establish function,
    //where all other resource groups establish (e.g. perennials). This function is no longer required.
	ForEachGroup(rg)
	{
		g = RGroup[rg];
		g->pr = ZRO(g->res_avail) ? 0. : g->res_required / g->res_avail;
		//if (g->max_age == 1)
		//{
		//	g->relsize = _add_annuals(rg, g->pr, add_seeds);
		//}
	}

	rgroup_ResPartIndiv();

}

static RealF _add_annuals(const GrpIndex rg, const SppIndex sp, const RealF lastyear_relsize) {
    /*======================================================*/
    /* check regen_ok flag.  if true, apply establishment and
     * add to seedbank.  Otherwise, add 0 to seedbank and skip
     * adding plants this year.  
     * reset group size to account for additions, or 0 if none.
     */
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
    ( (sizeof(x) == sizeof(float)) \
  ? ((x)>-xF_DELTA && (x)<xF_DELTA) \
  : ((x)>-xD_DELTA && (x)<xD_DELTA) )
    IntU i;
    RealF x,
    newsize, sumsize;
    GroupType *g;
    SpeciesType *s;
    Bool forced;

    g = RGroup[rg];
    assert(g->max_age == 1);

    s = Species[sp];
    newsize = 0.0;
    forced = FALSE;

    /* force addition of new propagules */
    _add_annual_seedprod(sp, lastyear_relsize);
    forced = TRUE;

    x = 0.;
     if (RandUni() <= s->seedling_estab_prob) {
        x = (g->regen_ok) ? _get_annual_maxestab(sp) : 0.;
        //printf("g->regen_ok: %d  , x = %.5f %\n", g->regen_ok, x);
    }

    if (GT(x, 0.)) {
        if (ZERO(lastyear_relsize)) {
            newsize = RandUniRange(1, x);
            //printf("newsize   =%0.5f \n", newsize);
        }
        else {
            //printf("lastyear_relsize=%0.5f \n", lastyear_relsize);
            newsize = x * lastyear_relsize;
            //printf("newsize   =%0.5f \n", newsize);
        }
    }
    return newsize;
}

static RealF _get_annual_maxestab(SppIndex sp) {
    /*======================================================*/
    /* Get the maximum number of viable seeds from the seedbank that can 
     establish this year*/
    IntU i;
    RealF sum = 0.;
    SpeciesType *s = Species[sp];

    for (i = 1; i <= s->viable_yrs; i++)
        sum += s->seedprod[i - 1] / pow(i, s->exp_decay);

    return sum;
}

static void _add_annual_seedprod(SppIndex sp, RealF lastyear_relsize) {
    /*======================================================*/
    /* Add seeds to the seedbank this year and increment viable years for seeds.*/

    SpeciesType *s = Species[sp];
    IntU i;


    //incrementing the number of viable years
    for (i = s->viable_yrs - 1; i > 0; i--) {
        if (i == 1 || i == 2) {
            //printf("Species name=%s , old Array values for index i=%d, value=%.5f \n", s->name, i, s->seedprod[i]);
        }

        s->seedprod[i] = s->seedprod[i - 1];
    }

    //printf("Species name=%s ,old Array array 0 index i=%d, value =%.5f \n", s->name, i, s->seedprod[i]);

    //If the current year is year 1 of the simulation, then the number of seeds added is a random number draw between
    //1 and the maximum number of seedlings that can establish. Otherwise, this year's seed production is a function of the maximum 
    //number of seedlings that can establish and last year's species relative size.
    if (Globals.currYear == 1) {
        s->seedprod[i] = RandUniRange(1, s->max_seed_estab);
        //printf("Species name=%s ,currYear =1 so new calculated value s->seedprod[%d]= %.5f , s->max_seed_estab =%d\n", s->name, i, s->seedprod[i], s->max_seed_estab);
    } else {
        s->seedprod[i] = s->max_seed_estab * lastyear_relsize;
        //printf("Species name=%s ,currYear=%d  so new calculated value s->seedprod[%d]= %.5f , s->max_seed_estab =%d, lastyear_relsize=%.5f\n", s->name, Globals.currYear, i, s->seedprod[i], s->max_seed_estab, lastyear_relsize);
    }
}

static RealF _ppt2resource(RealF ppt, GroupType *g) {
    /*======================================================*/
    /* PURPOSE */
    /* Converts ppt in mm to a resource index for a group
     pointed to by g.
     Represents EQNs 2,3,4

     /* SCOPE
     Local routine called from rgroup_PartResources.

     /* HISTORY
     /* Chris Bennett @ LTER-CSU 2/12/2001

     /*------------------------------------------------------*/

    return (ppt * g->ppt_slope[Env.wet_dry] + g->ppt_intcpt[Env.wet_dry]);

}

/***********************************************************/
static void _res_part_extra(Bool isextra, RealF extra, RealF size[])
{
	/*======================================================*/
	/* PURPOSE */
	/* Determines if some groups have unused resources and
	 redistributes them to other groups that can use them
	 (if any).

	 See COMMENT 2 at the end of this file for the algorithm.

	 /* SCOPE
	 Local routine called from rgroup_PartResources.

	 /* HISTORY
	 /* Chris Bennett @ LTER-CSU 12/21/2000
	 Removed from rgroup_PartResources() to simplify that
	 routine.

	 15-May-03 (cwb) Adding code to accomodate resource-by-mm
	 when used with Soilwat, esp, see inclusion of 'space'.

	 16-July-16 (kap) These functions are not called when using SOILWAT,
	 despite the fact that their is separate code below for
	 when SOILWAT is running.
	 */

	/*------------------------------------------------------*/

	GrpIndex rg;
	GroupType *g; /* shorthand for RGroup[rg] */
	RealF req_prop, /* group's prop'l contrib to the total requirements */
	sum_size = 0., space; /* placeholder for 1 if by-mm or min_res_req otherwise */

	ForEachGroup(rg)
		sum_size += size[rg];

	ForEachGroup(rg)
	{
		g = RGroup[rg];
		if (ZRO(g->relsize))
			continue;
		if (isextra && !g->use_extra_res)
			continue;
		space = (UseSoilwat) ? 1.0 : g->min_res_req;

		req_prop = size[rg] / sum_size;

		if (isextra && g->use_extra_res && GT(g->xgrow, 0))
			g->res_extra = req_prop * extra / space;
		else
			g->res_avail += req_prop * extra / space;

	}

}

/***********************************************************/
void rgroup_ResPartIndiv(void)
{
	/*======================================================*/
	/* PURPOSE */
	/* The PR value used in the growth loop is now at the
	 individual level rather than the group level, so the
	 resource availability for the group is divided between
	 individuals of the group, largest to smallest.  Species
	 distinctions within the group are ignored.

	 The partitioning of resources to individuals is done
	 proportionally based on size so that individuals can
	 have their own PR value.  The reasoning is that larger
	 individuals should get more resources than their
	 less-developed competitors.

	 See COMMENT 3 at the end of this file for the algorithm.

	 /* HISTORY
	 /* Chris Bennett @ LTER-CSU 12/21/2000
	 Removed from rgroup_PartResources() to simplify that
	 routine.

	 7-Jul-02 - made this a global function to be called
	 by the STEPWAT code.  No other changes.

	 2-Mar-03 - changed indivs and res_prop to be local and
	 dynamically allocated.
	 - changed the premise of the routine.  see
	 algorithm comments in Comment 3 for more info.
	 */

	/*------------------------------------------------------*/

	GrpIndex rg;
	GroupType *g; /* shorthand for RGroup[rg] */
	IndivType **indivs, /* dynamic array of indivs in RGroup[rg] */
	*ndv; /* shorthand for the current indiv */
	IntS numindvs, n;
	RealF x, /* temporary multiplier */
	base_rem; /* remainder of resource after allocating to an indiv */

	/* -- apportion each group's resources to individuals */
	ForEachGroup(rg)
	{
		g = RGroup[rg];
		//if (g->max_age == 1 ) continue;  /* annuals don't have indivs */
		// removed to give annuals individuals (TEM 10-27-2015)
		if (!g->est_count)
			continue;

		/* --- allocate the temporary group-oriented arrays */
		indivs = RGroup_GetIndivs(rg, SORT_D, &numindvs);

		/* ---- assign indivs' availability, not including extra ---- */
		/*      resource <= 1.0 is for basic growth, >= extra         */
		/*      amount of extra, if any, is kept in g->res_extra      */
		/*    base_rem = fmin(1, g->res_avail);  */
		base_rem = g->res_avail;
		for (n = 0; n < numindvs; n++)
		{
			ndv = indivs[n];

			ndv->res_required = (ndv->relsize / g->max_spp) / g->estabs;
			if (GT(g->pr, 1.))
			{
				ndv->res_avail = fmin(ndv->res_required, base_rem);
				base_rem = fmax(base_rem - ndv->res_avail, 0.);
			}
			else
			{
				ndv->res_avail = ndv->grp_res_prop * g->res_avail;
			}
		}

		base_rem += fmin(0, g->res_avail - 1.0);

		/* --- compute PR, but on the way, assign extra resource */
		for (n = 0; n < numindvs; n++)
		{
			ndv = indivs[n];

			if (g->use_extra_res)
			{
				/* polish off any remaining resource not allocated to extra */
				ndv->res_avail +=
						ZRO(base_rem) ? 0. : ndv->grp_res_prop * base_rem;
				/* extra resource gets assigned here if applicable */
				if (GT(g->res_extra, 0.))
				{
					x = 1. - ndv->relsize;
					ndv->res_extra = (1. - x) * ndv->grp_res_prop
							* g->res_extra;
					ndv->res_avail += x * ndv->grp_res_prop * g->res_extra;
				}

			}

			/* ---- at last!  compute the PR value, or dflt to 100  */
			ndv->pr =
					GT(ndv->res_avail, 0.) ?
							ndv->res_required / ndv->res_avail : 100.;
		}

		Mem_Free(indivs);

	} /* end ForEachGroup() */

}

/***********************************************************/
void rgroup_Grow(void)
{
	/*======================================================*/
	/* PURPOSE */
	/* Main loop to grow all the plants.
	 */
	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */

	/*  7-Nov-03 (cwb) Annuals have a completely new method of
	 *     propogation, growth, and death.  Most notably, they
	 *     aren't subject to the growth mechanism here, so they
	 *     are now excluded from this code.  All of the annual
	 *     generation occurs in the PartResources() function.
	 *  7/13/2016 This comment by cwb is now outdated. We had added
	 functionality so that annuals now grow on an annaul basis.

	 /*------------------------------------------------------*/

	IntU j;
	GrpIndex rg;
	SppIndex sp;
	GroupType *g;
	SpeciesType *s;
	const RealF OPT_SLOPE = .95;
	RealF growth1, /* growth of one individual*/
	sppgrowth, /* sum of growth for a species' indivs */
	rate1, /* rate of growth for an individual*/
	tgmod, gmod; /* growth factor modifier*/
	IndivType *ndv; /* temp pointer for current indiv */
	/*------------------------------------------------------*/
	ForEachGroup(rg)
	{
		g = RGroup[rg];
		//if (g->max_age == 1) continue; /* annuals already taken care of */
		// removed to allow annuals to grow (TEM 10-27-2015))

		if (!g->est_count)
			continue;

		/* can't grow succulents if a wet year, so skip this group */
		if (g->succulent && Env.wet_dry == Ppt_Wet)
			continue;

		/*--------------------------------------------------*/
		/* for each non-annual species */
		/* grow individuals and increment size */
		/* all groups are either all annual or all perennial */
		ForEachEstSpp(sp, rg, j)
		{
			s = Species[sp];

			sppgrowth = 0.0;
			if (!Species[sp]->allow_growth)
				continue;

			/* Modify growth rate by temperature
			 calculated in Env_Generate() */
			if (s->tempclass != NoSeason)
				tgmod = Env.temp_reduction[s->tempclass];

			/*------------------------------------------------*/
			/* now grow the individual plants of current species*/
			ForEachIndiv(ndv, s)
			{

				/* modify growth rate based on resource availability*/
				/* deleted EQN 5 because it's wrong. OPT_SLOPE was formerly 0.05, which resulted in gmod values that were too low */
				/* Now the values fro gmod range between 0.05 and 0.99 similiar to Coffin and Lauenroth 1990 */

				gmod = 1.0 - OPT_SLOPE * min(1.0, ndv->pr);
				if (GT(ndv->pr, 1.0))
					gmod /= ndv->pr;

				gmod *= tgmod;

				if (ndv->killed && RandUni() < ndv->prob_veggrow)
				{
					/* if indiv appears killed it was reduced due to low resources */
					/* last year. it can veg. prop. this year but couldn't last year.*/
					growth1 = s->relseedlingsize
							* RandUniRange(1, s->max_vegunits);
					rate1 = growth1 / ndv->relsize;
					ndv->killed = FALSE;

				}
				else
				{
					/* normal growth: modifier times optimal growth rate (EQN 1)in Coffin and Lauenroth 1990 */
					rate1 = gmod * s->intrin_rate * (1.0 - ndv->relsize);
					growth1 = rate1 * (1.0 - ndv->relsize);
				}
				ndv->relsize += growth1;
				ndv->growthrate = rate1;
				sppgrowth += growth1;
			} /*END ForEachIndiv */

			Species_Update_Newsize(sp, sppgrowth);

		} /* ENDFOR j (for each species)*/

		_extra_growth(rg);

	} /* END ForEachGroup(rg)*/

}

/***********************************************************/
static void _extra_growth(GrpIndex rg)
{
	/*======================================================*/
	/* PURPOSE */
	/* When there are resources beyond the minimum necessary
	 for "optimal" growth, ie, use of eqn 5, the extra
	 resources are converted to superfluous growth that
	 only counts for the current year and is removed at
	 the beginning of the next year.

	 */
	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 11/8/2000            */
	/* 3/14/01 - made this a local function called from
	 Grow() and put killextragrowth() in Mort_*
	 3/27/01 - extra growth now occurs in individuals
	 rather than at the species level.  Also accounts
	 for the scaling by min_res_req.

	 * 7-Nov-03 (cwb) Annuals don't get extra growth.  Everything
	 *     is accounted for by PR in PartResources().

	 */

	/*------------------------------------------------------*/
	Int j;
	RealF extra, indivpergram;
	GroupType *g;
	SpeciesType *s;
	IndivType *ndv;
	SppIndex sp;

	//if ( RGroup[rg]->max_age == 1) return;
	// removed to allow extra growth in annuals (TEM 10-27-2015))
	if (ZRO(RGroup[rg]->xgrow))
		return;
	if (!RGroup[rg]->use_extra_res)
		return;

	g = RGroup[rg];

	ForEachEstSpp(sp, rg, j)
	{
		s = Species[sp];

		indivpergram = 1.0 / s->mature_biomass;
		ForEachIndiv(ndv, s)
		{
			extra = ndv->res_extra * g->min_res_req * Env.ppt * g->xgrow;
			s->extragrowth += extra * indivpergram;
		}
		Species_Update_Newsize(s->sp_num, s->extragrowth);
	}

}
/***********************************************************/
void rgroup_Establish(void) {
    /*======================================================*/
    /* PURPOSE */
    /* Determines which and how many species can establish
     in a given year.

     For each species in each group, check that a uniform
     random number between 0 and 1 is less than the species'
     establishment probability.
     a) If so, return a random number of individuals,
     up to the maximum allowed to establish for the
     species.  This is the number of individuals in
     this species that will establish this year.
     b) If not, continue with the next species.

     */

    /* HISTORY */
    /* Chris Bennett @ LTER-CSU 6/15/2000            */
    /*   cwb (3/7/01) - Complete rewrite that simplifies the
     process and attempts to make the total number of
     seedlings established match the probability of
     establishment.
     cwb (3/30) - Actually the probability of establishment
     emulates the occurrance of microsite conditions that
     allow establishment.  There may be more seedlings
     established than indicated by the probability.

     * 7-Nov-03 (cwb) Adding the new algorithm to handle annuals.
     *   It's more complicated than before (which didn't really
     *   work) so annuals are now added in the PartResources()
     *   function.  Only perennials are added here.
     *
     *   Also, there's now a parameter to define the start year
     *   of establishment for perennials.
     * 
     * KAP: Annual establishment now occurs here instead of in PartResources
     */
    /*------------------------------------------------------*/
    IntS i, num_est; /* number of individuals of sp. that establish*/
    GrpIndex rg;
    SppIndex sp;
    GroupType *g;
    RealF lastyear_relsize;

    /* Cannot establish if plot is still in disturbed state*/
    if (Plot.disturbed > 0) {
        ForEachGroup(rg)
        RGroup[rg]->regen_ok = FALSE;
        return; /* skip regen for all */
    }

    ForEachGroup(rg) {
        g = RGroup[rg];
        if (!g->use_me)
            continue;

        g->regen_ok = TRUE; /* default */

        if (Globals.currYear < RGroup[rg]->startyr) {
            g->regen_ok = FALSE;

        } else ///if ( g->max_age == 1 ) {
            /* see similar logic in mort_EndOfYear() for perennials */
            /// if ( GT( g->killfreq, 0.) ) {
            ///   if ( LT(g->killfreq, 1.0) ) {
            ///     if (RandUni() <= g->killfreq)
            ///       g->regen_ok = FALSE;
            ///   } else if ( (Globals.currYear - g->startyr) % (IntU)g->killfreq == 0) {
            ///     g->regen_ok = FALSE;
            ///   }
            /// }
            ///} else
            //above removed allow annuals to establish with other species (TEM 10-27-2015)
        {
            ForEachGroupSpp(sp, rg, i) {
                if (!Species[sp]->use_me)
                    continue;
                if (!Species[sp]->allow_growth)
                    continue;

                if (Species[sp]->max_age == 1) {
                    printf("Globals.currYear = %d, call to _add_annuals sp=%d Species[sp]->lastyear_relsize : %.5f \n", Globals.currYear, sp, Species[sp]->lastyear_relsize);
                    num_est = _add_annuals(rg, sp, Species[sp]->lastyear_relsize);
                    //printf("num_est for annuals=%d \n",num_est);
                }                    
                else {

                    num_est = Species_NumEstablish(sp);
                }

                if (num_est) {
                    /* printf("%d %d %d %d\n",
                     Globals.currIter, Globals.currYear, sp, num_est); */

                    Species_Add_Indiv(sp, num_est);
                    species_Update_Estabs(sp, num_est);
                }
            }
        }
    }
}

/***********************************************************/
void rgroup_IncrAges(void)
{
	/*======================================================*/
	/* PURPOSE */
	/*  Increment ages of individuals in a resource group.
	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */
	/*   8-Nov-03 (cwb) added check for annuals */

	/*------------------------------------------------------*/
	Int j;
	GrpIndex rg;
	SppIndex sp;
	IndivType *ndv;

	ForEachGroup(rg)
	{
		if (RGroup[rg]->max_age == 1)
			continue;

		ForEachEstSpp( sp, rg, j)
		{
			ForEachIndiv( ndv, Species[sp])
			{
				ndv->age++;
				if (ndv->age > Species[ndv->myspecies]->max_age)
				{
					LogError(logfp, LOGWARN,
							"%s grown older than max_age (%d > %d). Iter=%d, Year=%d\n",
							Species[ndv->myspecies]->name, ndv->age,
							Species[ndv->myspecies]->max_age, Globals.currIter,
							Globals.currYear);
				}
			}
		}
	}

}

/***********************************************************/
void RGroup_Update_Newsize(GrpIndex rg)
{
	/*======================================================*/
	/* PURPOSE */
	/*   Relative size for a group is 1.0 if all the group's
	 *   species are established and size 1.0;
	 */

	/* HISTORY */
	/*   Chris Bennett @ LTER-CSU 5-Apr-2003            */
	/*     Called from Species_Update_Newsize(), but safe to call
	 *     from anywhere.
	 *   7-Nov-03 (cwb) Added condition for annuals.  Annuals
	 *     don't use the linked list of indiv objects because they
	 *     come and go with each time step, so we just add their
	 *     estimated size to the group size and move on.  See also
	 *     Species_UpdateNewSize() and other calls with *annual*
	 *     somewhere in them.
	 */

	/*------------------------------------------------------*/
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
( (sizeof(x) == sizeof(float)) \
  ? ((x)>-xF_DELTA && (x)<xF_DELTA) \
  : ((x)>-xD_DELTA && (x)<xD_DELTA) )

	Int n;
	SppIndex sp;
	IndivType **indivs;
	IntS numindvs;
	RealF sumsize = 0.0;

	/* first get group's relative size adjusted for num indivs */
	/* ie, groupsize=1 when 1 indiv of each species is present */
	/* ie each indiv is an equivalent contributor, not based on biomass */
	ForEachEstSpp( sp, rg, n)
		sumsize += Species[sp]->relsize;

	if (RGroup[rg]->est_count < 0)
		RGroup[rg]->est_count = 0;

	if (RGroup[rg]->est_count == 0 || LT(sumsize, 0.0))
	{
		RGroup[rg]->relsize = 0.0;
	}
	else
	{
		//For calculating rgroup relSize, sumsize should be divide by no of current established species in rgroup rather than total no of species in rgroup.
		RGroup[rg]->relsize = sumsize / (RealF) RGroup[rg]->est_count;
	}

	///if (RGroup[rg]->max_age != 1) {
	/* compute the contribution of each indiv to the group's size */
	indivs = RGroup_GetIndivs(rg, SORT_0, &numindvs);
	for (n = 0; n < numindvs; n++)
		indivs[n]->grp_res_prop = indivs[n]->relsize / sumsize;
	Mem_Free(indivs);
	///}

	/* double check some assumptions */
	if (RGroup[rg]->est_count < 0)
		RGroup[rg]->est_count = 0;
	if (ZERO(RGroup[rg]->relsize))
		RGroup[rg]->relsize = 0.0;

#undef xF_DELTA
#undef xD_DELTA
#undef ZERO
}

/***********************************************************/
RealF RGroup_GetBiomass(GrpIndex rg)
{
	/*======================================================*/
	/* PURPOSE */
	/*   Convert relative size to biomass for a resource group.
	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */
	/*   10-Apr-03 - cwb - return sum of plant sizes times density
	 *       for the per-plot biomass.
	 *
	 *   8-Dec-03 (cwb) actually, what we want is the total biomass
	 *       of the group on the plot, so multiplying by density
	 *       doesn't make sense.  density limits the rsize
	 *       which already scales the per-plot values.
	 */

	/*------------------------------------------------------*/
	Int j;
	SppIndex sp;
	RealF biomass = 0.0;

	if (RGroup[rg]->est_count == 0)
		return 0.0;
	ForEachEstSpp( sp, rg, j)
	{
		biomass += Species[sp]->relsize * Species[sp]->mature_biomass;
	}

	return biomass;
}

/***********************************************************/
GrpIndex RGroup_Name2Index(const char *name)
{
	/*======================================================*/
	/* PURPOSE */
	/*   Simple utility to find the index of a group name string
	 *
	 *   returns group index number if name is a valid group name
	 *   otherwise returns 0.  This works because the
	 *   ForEachGroup macro starts at 1.
	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */

	/*------------------------------------------------------*/
	GrpIndex i, rg = -1;

	ForEachGroup(i)
	{
		if (strcmp(name, RGroup[i]->name) == 0)
		{
			rg = i;
			break;
		}
	}
	return (rg);
}

/***********************************************************/
static GroupType *_create(void)
{
	/*======================================================*/

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */

	/*------------------------------------------------------*/
	GroupType *p;

	p = (GroupType *) Mem_Calloc(1, sizeof(GroupType), "__Create");

	return (p);

}

/***********************************************************/
GrpIndex RGroup_New(void)
{
	/*======================================================*/
	/* PURPOSE */
	/* Create a new resource group (life form) object and give
	 it the next consecutive identifier.

	 Initialization is performed in parm_RGroup_Init().

	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */

	/*------------------------------------------------------*/
	GrpIndex i = (GrpIndex) Globals.grpCount;

	if (++Globals.grpCount > MAX_RGROUPS)
	{
		LogError(logfp, LOGFATAL, "Too many groups specified (>%d)!\n"
				"You must adjust MAX_RGROUPS and recompile!",
		MAX_RGROUPS);
	}

	RGroup[i] = _create();
	return (GrpIndex) i;
}

/***********************************************************/
void rgroup_DropSpecies(SppIndex sp)
{
	/*======================================================*/
	/* PURPOSE */
	/* When a species associated with a resource group dies
	 out, it is dropped from the group so it will not be
	 processed unnecessarily.
	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */

	/*------------------------------------------------------*/

	IntS i, j;
	GrpIndex rg;
	Bool f = FALSE;

	rg = Species[sp]->res_grp;
	ForEachEstSpp2(rg, i)
	{
		if (RGroup[rg]->est_spp[i] == sp)
		{
			f = TRUE;
			break;
		}
	}

	/* close up the array around the dead spp */
	if (f)
	{ /* note the < symbol not <= */
		for (j = i; j < RGroup[rg]->est_count; j++)
		{
			RGroup[rg]->est_spp[j] = RGroup[rg]->est_spp[j + 1];
		}
		RGroup[rg]->est_count--;
	}
}

/***********************************************************/
void rgroup_AddSpecies(GrpIndex rg, SppIndex sp)
{
	/*======================================================*/
	/* PURPOSE */
	/*   When a species associated with a resource group becomes
	 established, it is added to the list of species and
	 otherwise linked to the group.
	 */
	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000            */

	/*------------------------------------------------------*/
	Int i;
	Bool f = FALSE;

	ForEachEstSpp2( rg, i)
	{
		if (RGroup[rg]->est_spp[i] == sp)
		{
			f = TRUE;
			break;
		}
	}
	if (!f)
		RGroup[rg]->est_spp[RGroup[rg]->est_count++] = sp;
}

/**************************************************************/
void rgroup_Extirpate(GrpIndex rg)
{
	/*======================================================*/
	/* PURPOSE */
	/* Kill a group catastrophically, meaning, kill all
	 individuals, remove their biomass (relsize), and
	 don't let them regenerate ever again.
	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 7/10/2000
	 cwb - 11/27/00 - rewrote the algorithm to be more
	 efficient, ie, fixed the quick and dirty code
	 from before.
	 */

	/*------------------------------------------------------*/

	SppIndex sp, i;

	ForEachGroupSpp(sp, rg, i)
	{
		Species_Kill(sp, 5);
		Species[sp]->seedling_estab_prob = 0.0;
		/*TMartyn 5.26.2015 - added the following section because annual biomass
		 values were not being updated to 0 even if seedling establishment
		 probability was chaged to 0. The following code says if the max age of
		 * a species is 1 it will change the seed establishment to 0 and will
		 * vastly increase the exp_decay value so that the seeds in the
		 * seed bank will decay within 1 year.  This is a work around and we
		 * should look more into this to make biomass values absolutely 0 at the
		 * year of expiration implementation. */
		if (Species[sp]->max_age == 1)
		{
			Species[sp]->max_seed_estab = 0.0;
			Species[sp]->exp_decay = 1000000;

		}
	}

	RGroup[rg]->extirpated = TRUE;

}

/**************************************************************/
void RGroup_Kill(GrpIndex rg)
{
	/*======================================================*/
	/* PURPOSE */
	/* Kill all individuals of all species in the group,
	 but allow them to regenerate.
	 */

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 11/27/2000            */

	/*------------------------------------------------------*/

	//printf("inside RGroup_Kill() rg=%d, RGroup[rg]->proportion_killed=%f \n",rg,RGroup[rg]->proportion_killed);
	Int i;

	ForEachEstSpp2( rg, i)
		Species_Proportion_Kill(RGroup[rg]->est_spp[i], 6,
				RGroup[rg]->proportion_killed);
}

/**********************************************************/
IndivType **RGroup_GetIndivs(GrpIndex rg, const char sort, IntS *num)
{
	/*======================================================*/
	/* PURPOSE */
	/* put all the individuals of the group into a list
	 and optionally sort the list.

	 - Returns an allocated list of indivs; calling routine must
	 free when appropriate.
	 - Puts the number of individuals in *num.


	 /* HISTORY */
	/* Chris Bennett @ LTER-CSU 12/15/2000
	 *
	 *  2-Mar-03 - cwb - removed requirement for list to be
	 *      pre-allocated.  Added code to allocate and return list.
	 /*------------------------------------------------------*/

	IntS j, i = 0;
	size_t i_size = sizeof(IndivType **);
	SppIndex sp;
	IndivType *ndv, **nlist;

	ForEachEstSpp(sp, rg, j)
		ForEachIndiv(ndv, Species[sp])
			i++;
	nlist = Mem_Calloc(i, i_size, "Rgroup_GetIndivs(nlist)");

	i = 0;
	ForEachEstSpp(sp, rg, j)
		ForEachIndiv(ndv, Species[sp])
			nlist[i++] = ndv;

	*num = i;
	if (i > 0 && sort)
		Indiv_SortSize(sort, (size_t) i, nlist);

	return nlist;
}

#ifdef DEBUG_MEM
#include "myMemory.h"
/*======================================================*/
void RGroup_SetMemoryRefs( void)
{
	/* when debugging memory problems, use the bookkeeping
	 code in myMemory.c
	 This routine sets the known memory refs in this module
	 so they can be  checked for leaks, etc.  All refs will
	 have been cleared by a call to ClearMemoryRefs() before
	 this, and will be checked via CheckMemoryRefs() after
	 this, most likely in the main() function.

	 EVERY dynamic allocation must be noted here or the
	 check will fail (which is the point, to catch unknown
	 or missing pointers to memory).
	 */
	GrpIndex rg;

	ForEachGroup(rg)
	{
		NoteMemoryRef(RGroup[rg]);
		NoteMemoryRef(RGroup[rg]->kills); /* this is set it params() */
	}

}

#endif

/* COMMENT 1 - Algorithm for rgroup_PartResources() */
/*
 * Assign minimum resources to each group based on relative
 * size (which can be greater than 1.0).  Resource is computed
 * by scaling the PPT according to the slope and intercept
 * parameters in the groups input file.  Average PPT should
 * yield a resource value of 1.0.  Any resource above 1.0 is
 * committed to extra resource and is allocated in
 * _PartResExtra().
 *
 * Resource required is equal to the relative size of the group;
 * resource available is equivalent to the relative size truncated
 * at 1.0, ie, this is how a group becomes resource limited.
 * For example, a group of 1 three-quarter-sized individual
 * requires min_res_req times .75 times the resource (scaled PPT);
 * for an average PPT year this produces a resource availability
 * of  0.75.  The minimum availability needed for maintenance is
 * equivalent to the relative size <= 1.0.  If the size is < 1,
 * the plant is too small to use all of the potentially available
 * resource, so the difference is added to the pool of extra
 * resource.  If the size is > 1.0, then the plant is too large
 * to use the minimum available. On the other hand, if this is a
 * below average PPT year and the size is 1.0 (or more) the group
 * only gets min_res_req times resource < 1.0 , so groups with
 * size > space * resource will have stretched resources.
 *
 * However if this is an above average PPT year, the groups start
 * out with no more than min_res_req amount of the resource even
 * if the requirements (and size) are higher.  In any case,
 * _PartResExtra() will determine if there is enough extra
 * from other smaller groups and or above average PPT to make up
 * the difference.
 */

/*  COMMENT 3 - Algorithm for rgroup_ResPartIndiv()
 *  2-Mar-03
 *
 Originally (before 2-Mar-03), the premise for the entire scheme
 of growth modification dependent upon the ratio of required to
 available resources was probably faulty in that availability was
 limited by the size of the indiv (or group).  In retrospect,
 this doesn't make sense.  Rather the availability is defined by
 the amount of resource, irrespective of whether plants are
 available to use it, ie, min(resource, 1.0). Generally this can
 be interpreted as "free range" and any able-bodied plant is
 welcome to take what it can.  For individuals within a group
 though, it's necessary to define availability as the total
 available resource assigned to indivs based on their
 proportional size.  All the available resource is "credited" to
 the individuals based on their size (a competitive feature).

 Whereas in the past, size-limited available resource was partitioned
 according to a "cup"-based method (meaning each individual was
 like a cup of relsize filled in decreasing order until resource
 was depleted), now that isn't necessary because, and this is
 very important, small groups (< resource) have ample resource
 to meet each indiv's minimum need so all indivs get enough.

 However, for large (>= resource) groups, indivs still get available
 resource by the cup method.  For example, if this is a dry year
 with full requirements, the larger indivs get as much as they
 need (relsize amt) until the resource runs out, after which the
 smallest indivs are considered to be in severe resource deprivation
 for that year.

 For groups that are allowed to use extra resource, the pooled
 extra resource is allocated proportionally.  If a plant is
 allowed to exhibit extra (superfluous yearly) growth, the amount
 assigned is directly related to the size of the plant, ie,
 smaller plants get more of the extra resource applied to
 persistent growth (relsize) whereas larger plants get more applied
 to extra (leafy) growth.

 This general approach should work equally well whether the resource
 comes from SOILWAT or the MAP equation because it's based on the
 group PR value which is determined before this routine.  Note,
 though, that the implication is that both methods of resource
 creation imply that availability is not tied to relsize.

 */

