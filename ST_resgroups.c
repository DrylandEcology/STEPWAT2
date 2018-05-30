/********************************************************/
/*  Source file: resgroups.c
    Type: module
    Application: STEPPE - plant community dynamics simulator
    Purpose: This is the module that executes the growth
            functions and otherwise manages the bookkeeping
            at the group level.
   History:
      (6/15/2000) -- INITIAL CODING - cwb
 
* =================================================== *
*                INCLUDES / DEFINES                   *
* --------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "myMemory.h"
#include "rands.h"
#include "generic.h"
#include "sw_src/filefuncs.h"
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
static void _res_part_extra(RealF extra, RealF size[]);
static GroupType *_create(void);
static void _extra_growth(GrpIndex rg);
static void _add_annual_seedprod(SppIndex sp, RealF pr);
static RealF _get_annual_maxestab(SppIndex sp);
static RealF _add_annuals(const GrpIndex rg, const RealF g_pr,
		const Bool add_seeds);

/****************** Begin Function Code ********************/
/***********************************************************/

void rgroup_PartResources(void)
{
/* Partition resources for this year among the resource
groups. The allocation happens in three steps: basic
group allocation, partitioning of extra resources in _res_part_extra(),
and further partitioning to individuals in the group _ResPartIndiv().
See COMMENT 1. at the end of this file for the algorithm.*/

	GrpIndex rg;
	Bool noplants = TRUE;
	const Bool no_seeds = FALSE;
	GroupType *g; /* shorthand for RGroup[rg] */

	/* ----- distribute basic (minimum) resources */
	ForEachGroup(rg)
	{
		g = RGroup[rg];

		//next two lines have been deleted in annual_final branch
		if (g->max_age == 1)
			g->relsize = _add_annuals(rg, 1.0, no_seeds);

		g->res_required = RGroup_GetBiomass(rg);
		g->res_avail = SXW_GetTranspiration(rg); 
		//printf("g->res_avail = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_avail); 
        //printf("g->res_required = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_required); 
        
		//A check and reset of res_required and res_avail, this should never happen
		if(ZRO(g->res_avail) && g->res_required > 0)
		{
			g->res_required = g->estabs;
			g->res_avail = 1;
			LogError(logfp, LOGWARN, "RGroup %s : res_avail is Zero and res_required > 0", g->name);
		}
        
		/* Annuals seem to have a artificial limit of 20. We do Annuals here differently. 
		This will be re-evaluated when annual_final is merged in. I would opt to remove,
        but additional testing is required with the new annual code. */
		if(g->max_age == 1)
		{
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

	  /* If relsize>0, reset noplants from TRUE to FALSE and if noplants=TRUE, exit from the loop */
		if (GT(g->relsize, 0.))
			noplants = FALSE;

	} /* End ForEachGroup(rg) */

	if (noplants)
		return;

	//calculate PR at the group level: resources required/resources available  
	ForEachGroup(rg)
	{
		g = RGroup[rg];
		g->pr = ZRO(g->res_avail) ? 0. : g->res_required / g->res_avail;
		//printf("g->pr = %f\n,Group = %s \n",RGroup[rg]->name,  g->pr);
	}

 	//partition resources to individuals and determine 'extra' resources
	rgroup_ResPartIndiv();

}

static RealF _add_annuals(const GrpIndex rg, const RealF g_pr,
		const Bool add_seeds)
{
	/*======================================================*/
	/* if add_seeds==FALSE, don't add seeds to the seedbank or
	 * add species, but do calculations required to return a
	 * temporary group relative size for resource allocation.
	 * if add_seeds==TRUE, we should have done the above and now
	 * we really want to add to the biomass and seedbank.
	 *
	 * check regen_ok flag.  if TRUE, apply establishment and
	 * add to seedbank.  Otherwise, add 0 to seedbank and skip
	 * adding plants this year.  We also check probability of establishment to
	 * account for introduction of propagules, in which case
	 * we can't add to the seedbank twice.
	 *
	 * reset group size to account for additions, or 0 if none.
	 */
	const RealF PR_0_est = 20.0; /* PR at which no seeds can establish */
	IntU i;
	RealF x, estabs, /* number of estabs predicted by seedbank */
	pr_inv, /* inverse of PR as the limit of estabs */
	newsize, sumsize;
	SppIndex sp;
	GroupType *g;
	SpeciesType *s;
	Bool forced;

	g = RGroup[rg];

	assert(g->max_age == 1);

	if (!g->use_me)
		return 0.0;

	sumsize = 0.0;
	ForEachGroupSpp(sp,rg,i)
	{
		s = Species[sp];
		newsize = 0.0;
		forced = FALSE;
		// printf("add_seeds=%d \n",add_seeds);
		if (!s->use_me)
			continue;
		//if add seeds is F and a random number drawn from the uniform distribution is less than prob of seed establishment
		if (!add_seeds && RandUni() <= s->seedling_estab_prob)
		{
			/* force addition of new propagules */
			//if regen_ok is TRUE, pass the g_pr parameter, if regen_ok is False, pass -1
			_add_annual_seedprod(sp, (g->regen_ok) ? g_pr : -1.);
			forced = TRUE;
		}

		//if regen_ok is T then x = values coming from _get_annual_maxestab, if regen ok is F then
		//x=0.
		x = (g->regen_ok) ? _get_annual_maxestab(sp) : 0.;
		if (GT(x, 0.))
		{
			estabs = (x - (x / PR_0_est * g_pr)) * exp(-g_pr);
			pr_inv = 1.0 / g_pr;
			newsize = fmin(pr_inv, estabs);
			if (add_seeds)
			{
				rgroup_AddSpecies(rg, sp);
				Species_Update_Newsize(sp, newsize);
			}
		}

		if (add_seeds && !forced)
			_add_annual_seedprod(sp, (ZRO(x)) ? -1. : g_pr);

		sumsize += newsize;
	}

	/*if (g->regen_ok == TRUE) {
	 printf("Regen_ok %u\n",rg);}*/

	/*will print out a flag of Regen_ok with the resourse group number (numbered
	 * in the order of turned on groups in the "rgroup.in" file */

	return (sumsize / (RealF) g->max_spp);

}

static RealF _get_annual_maxestab(SppIndex sp)
{
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

static void _add_annual_seedprod(SppIndex sp, RealF pr)
{
	/*======================================================*/
	/* Add seeds to the seedbank this year and increment viable years for seeds.
         * negative pr means add 0 seeds to seedbank */

	SpeciesType *s = Species[sp];
	IntU i;
// incrementing the number of viable years
	for (i = s->viable_yrs - 1; i > 0; i--)
	{
		s->seedprod[i] = s->seedprod[i - 1];
	}

	//compare pr (actual value or -1) to zero, if less than zero, then make it zero. If greater than 0
	//than multiple max_seed_estab * exp(-pr)
	s->seedprod[i] = LT(pr, 0.) ? 0. : s->max_seed_estab * exp(-pr);
	//s->seedprod[i] = LT(pr, 0.) ? 0. : s->max_seed_estab * 1/pr ;
	//s->seedprod[i] = LT(pr, 0.) ? 0. : s->max_seed_estab * s->relsize;
}

/***********************************************************/
static void _res_part_extra(RealF extra, RealF size[]) {
    /*======================================================*
     * PURPOSE *
     * Partitions "extra" resources to other groups that can use them
     * (if any).*/

    GrpIndex rg;
    GroupType *g; /* shorthand for RGroup[rg] */
    RealF req_prop, /* group's prop'l contrib to the total requirements */
            sum_size = 0.; /* summed sizes of all functional groups */

    ForEachGroup(rg) {
        g = RGroup[rg];
        if (ZRO(g->relsize))
            continue;
        if (g->est_count == 0)
            continue;

        //where size = size_obase, which is biomass(g/m2) or 0 if the group can't use resources       
        sum_size += size[rg];
        //printf("size[rg]  = %f\n,Rgroup = %s \n",RGroup[rg]->name, size[rg]);
        //printf("sum_size  = %f\n,Rgroup = %s \n",RGroup[rg]->name, sum_size);
    } /* End ForEachGroup(rg) */
    //printf("sum_size  = %f\n", sum_size);

    ForEachGroup(rg) {
        g = RGroup[rg];
        if (ZRO(g->relsize))
            continue;
        if (g->est_count == 0)
            continue;

        // checking to make sure not dividing by 0
        if (sum_size == 0.)
            req_prop = 0.;
        else
            req_prop = size[rg] / sum_size;

        if (g->use_extra_res)
            g->res_extra = req_prop * extra;

        else
            g->res_extra = 0.;
        //printf("res_extra = %f\n,Rgroup = %s \n",RGroup[rg]->name,g->res_extra);

    } /* End ForEachGroup(rg) */

}
/***********************************************************/
void rgroup_ResPartIndiv(void) {
    /* The PR value used in the growth loop is now at the
     individual level rather than the group level, so the
     resource availability for the group is divided between
     individuals of the group, largest to smallest. Species
     distinctions within the group are ignored.The partitioning of resources
     to individuals is done proportionally based on size so that individuals
     should get more resources than their less-developed competitors.
     Chris Bennett @ LTER-CSU 12/21/2000*/

    GrpIndex rg;
    GroupType *g; /* shorthand for RGroup[rg] */
    SppIndex sp;
    Int j;
    IndivType **indivs, /* dynamic array of indivs in RGroup[rg] */
            *ndv; /* shorthand for the current indiv */
    IntS numindvs, n;
    RealF base_rem = 0., /* remainder of resource after allocating to an indiv */
            xtra_obase = 0., /* summed extra resources across all groups */
<<<<<<< HEAD
            size_obase[MAX_RGROUPS] = {0}; /* total res. contrib. if xtra_obase */
    
=======
            size_base[MAX_RGROUPS] = {0}, /* total res. contrib to base, all groups */
    size_obase[MAX_RGROUPS] = {0}; /* total res. contrib. if xtra_obase */
    const Bool do_extra = TRUE; /* monikers for extra resource partitioning */

>>>>>>> parent of 26e6f3d... Resolving issue #117
    /* -- apportion each group's normal resources to individuals */
    ForEachGroup(rg) {
        g = RGroup[rg];
        if (g->est_count == 0)
            continue;

        /* --- allocate the temporary group-oriented arrays */
        indivs = RGroup_GetIndivs(rg, SORT_D, &numindvs);

        /* Set resources available at the group level to base_rem */
        base_rem = g->res_avail;
        //printf("g->res_avail = %f\n, Group = %s \n", RGroup[rg]->name, g->res_avail);

        ForEachEstSpp(sp, rg, j) {
            for (n = 0; n < numindvs; n++) {
                ndv = indivs[n];

                /* Calculate resources required for each individual in terms of biomass */
                ndv->res_required = ndv->relsize * Species [sp]->mature_biomass;
                //printf("ndv->res_required = %f\n, Species = %s \n", Species[sp]->name, ndv->res_required);
                //printf("ndv->relsize = %f\n, Species = %s \n", Species[sp]->name, ndv->relsize);

                /* Calculate resources available for each individual based on res_required */
                ndv->res_avail = fmin(ndv->res_required, base_rem);
                //printf("ndv->res_avail = %f\n", ndv->res_avail);

                /* Remaining extra resource for each individual (if any) */
                base_rem = fmax(base_rem - ndv->res_avail, 0.);
                //printf("base_rem = %f\n", base_rem);
            }
        } /* end ForEachEstSpp() */

        //sum "extra" resources for all functional groups
        xtra_obase += base_rem;
        //printf("xtra_obase = %f\n", xtra_obase);

        /*check to see if each functional group can use extra resources,
        if so, then pass the biomass of the FG, otherwise pass 0 in _res_part_extra */
        size_base[rg] = RGroup_GetBiomass(rg);
        size_obase[rg] = (g->use_extra_res) ? size_base[rg] : 0.;
        //printf("size_obase = %f\n", size_obase[rg]);

        Mem_Free(indivs);

    } /* end ForEachGroup() */

    //assign extra resources to functional groups that can use them
    _res_part_extra(xtra_obase, size_obase);

    /* now loop back through all individuals in each group, assign extra resource
     * and calculate PR at the individual level. This extra resource will be used to
     * increment plant sizes. Extra resources that apply to superfluous biomass increment
     * that will ultimately be killed at the end of the year is assigned in _extra_growth */
    ForEachGroup(rg) {
        g = RGroup[rg];
        if (g->est_count == 0)
            continue;
        
        // Check to see if there are extra resources and if the group can use them
        if (!g->use_extra_res)
            continue;
        if (ZRO(g->res_extra))
            continue;
        
        //printf("g->res_extra = %f\n, RGroup= %s \n", RGroup[rg]->name, g->res_extra);

        /* --- allocate the temporary group-oriented arrays */
        indivs = RGroup_GetIndivs(rg, SORT_D, &numindvs);

        ForEachEstSpp(sp, rg, j) {

            /* --- compute PR, but on the way, assign extra resource */
            for (n = 0; n < numindvs; n++) {
                ndv = indivs[n];

                    //printf("ndv->res_avail before  = %f\n", ndv->res_avail);

                    // if individuals already have the resources they require do not assign extra
                    if (ndv->res_avail == ndv->res_required) {
                        ndv->res_extra = 0.0;
                        //printf("ndv->res_extra = %f\n", ndv->res_extra);    
                    }    
                        /* assign extra resource as the difference of what is required
                         * by the individual - what was assigned to the individual through
                         * partitioning of normal resources */
                    else {
                        ndv->res_extra = ndv->res_required - ndv->res_avail;
                        //printf("ndv->res_required = %f\n", ndv->res_required);
                        //printf("ndv->res_avail = %f\n", ndv->res_avail);
                        //printf("ndv->res_extra = %f\n", ndv->res_extra);

                        /* Updated resources available for each individual to include extra*/
                        ndv->res_avail += ndv->res_extra;
                        //printf("ndv->res_avail after = %f\n", ndv->res_avail);

                        /* Remaining extra resources */
                        g->res_extra = fmax(g->res_extra - ndv->res_extra, 0.);
                        //printf("g->res_extra in loop = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_extra); 
                    }
                
                /* Calculate the PR value, or dflt to 100, used in growth module */
                ndv->pr = GT(ndv->res_avail, 0.) ? ndv->res_required / ndv->res_avail : 100.;
                //printf("ndv->pr  = %f\n", ndv->pr);
            }
        } /* end ForEachEstSpp() */

        //printf("g->res_extra after = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_extra); 

        Mem_Free(indivs);

    } /* end ForEachGroup() */
}
/***********************************************************/
void rgroup_Grow(void) {
    /*======================================================*
     * PURPOSE *
     * Main loop to grow all the plants.
     *
     * HISTORY *
     * Chris Bennett @ LTER-CSU 6/15/2000            
     *  7-Nov-03 (cwb) Annuals have a completely new method of
     *     propogation, growth, and death.  Most notably, they
     *     aren't subject to the growth mechanism here, so they
     *     are now excluded from this code.  All of the annual
     *     generation occurs in the PartResources() function.
     *  7/13/2016 This comment by cwb is now outdated. We had added
     * functionality so that annuals now grow on an annual basis.*/

    IntU j;
    GrpIndex rg;
    SppIndex sp;
    GroupType *g;
    SpeciesType *s;
    const RealF OPT_SLOPE = .05; /*from Coffin and Lauenroth 1990, EQN 5 */
    RealF growth1, /* growth of one individual */
            sppgrowth, /* sum of growth for a species' indivs */
            rate1, /* rate of growth for an individual */
            tgmod, /* temperature growth factor modifier */
            gmod; /* growth factor modifier */
    IndivType *ndv; /* temp pointer for current indiv */

    ForEachGroup(rg) {
        g = RGroup[rg];
        //if (g->max_age == 1) continue; /* annuals already taken care of */
        // removed to allow annuals to grow (TEM 10-27-2015))

        if (g->est_count == 0)
            continue;

        /* can't grow succulents if a wet year, so skip this group */
        if (g->succulent && Env.wet_dry == Ppt_Wet)
            continue;

        /* grow individuals and increment size */
        ForEachEstSpp(sp, rg, j) {
            s = Species[sp];

            sppgrowth = 0.0;
            if (!Species[sp]->allow_growth)
                continue;

            /* Modify growth rate by temperature calculated in Env_Generate() */
            if (s->tempclass != NoSeason)
                tgmod = Env.temp_reduction[s->tempclass];

            /* now grow the individual plants of current species*/
            ForEachIndiv(ndv, s) {
                /* Modify growth rate based on resource availability. Values for gmod range 
                between 0.05 and 0.95 similiar to Coffin and Lauenroth 1990 */
                gmod = 1.0 - OPT_SLOPE;

                if (GT(ndv->pr, 1.0))
                    gmod /= ndv->pr;

                gmod *= tgmod;

                if (ndv->killed && RandUni() < ndv->prob_veggrow) {
                    /* if indiv appears killed it was reduced due to low resources */
                    /* last year. it can veg. prop. this year but couldn't last year.*/
                    growth1 = s->relseedlingsize * RandUniRange(1, s->max_vegunits);
                    rate1 = growth1 / ndv->relsize;
                    ndv->killed = FALSE;
                    //printf("growth1 killed  = %f\n", growth1);
                } else {
                    /* normal growth: modifier times optimal growth rate (EQN 1)in Coffin and Lauenroth 1990 */
                    rate1 = gmod * s->intrin_rate * (1.0 - ndv->relsize);
                    growth1 = rate1 * ndv->relsize;
                }
                //printf("growth1  = %f\n", growth1);
                //printf("old ndv->relsize  = %f\n,Species = %s \n", Species[sp]->name, ndv->relsize);

                ndv->relsize += growth1;
                //printf("new ndv->relsize  = %f\n, Species = %s \n", Species[sp]->name,ndv->relsize);

                ndv->growthrate = rate1;

                sppgrowth += growth1;
                //printf("sppgrowth = %f\n, Species = %s \n", Species[sp]->name,sppgrowth);

            } /*END ForEachIndiv */

            Species_Update_Newsize(sp, sppgrowth);

        } /* ENDFOR j (for each species)*/

        _extra_growth(rg);

    } /* END ForEachGroup(rg)*/
}

/***********************************************************/
static void _extra_growth(GrpIndex rg) {
    /*======================================================*/
    /* PURPOSE */
    /* When there are resources beyond the minimum necessary
     for "optimal" growth, ie, use of eqn 5, the extra
     resources are converted to superfluous growth that
     only counts for the current year and is removed at
     the beginning of the next year. */
    /* HISTORY */
    /* Chris Bennett @ LTER-CSU 11/8/2000            */
    /* 3/14/01 - made this a local function called from
     Grow() and put killextragrowth() in Mort_*
     3/27/01 - extra growth now occurs in individuals
     rather than at the species level. */

    Int j;
    RealF extra, indivpergram;
    GroupType *g;
    SpeciesType *s;
    IndivType *ndv;
    SppIndex sp;

    g = RGroup[rg];

    if (!RGroup[rg]->use_extra_res)
        return;
    if (ZRO(g->res_extra))
        return;

    //printf("g->res_extra growth = %f\n,Group = %s \n",RGroup[rg]->name, g->res_extra); 

    ForEachEstSpp(sp, rg, j) {
        s = Species[sp];
        
        /* Clear extra for each species*/
        Species[sp]->extragrowth = 0.0;

        /* Calculate the proportion of maximum species biomass that is represented by 1 unit */
        indivpergram = 1.0 / s->mature_biomass;

        ForEachIndiv(ndv, s) {
            /* Clear extra for each individual*/
            extra = 0.0;

            /* Calculate the extra resource available to each individual based on size */
            ndv->res_extra = ndv->grp_res_prop * g->res_extra;
            //printf("ndv->res_extra = %f\n,Species = %s \n", Species[sp]->name, ndv->res_extra);
            //printf("ndv->grp_res_prop = %f\n,Species = %s \n", Species[sp]->name, ndv->grp_res_prop);

            /* Calculate the increment in size due to res_extra for each individual 
            Note, we have not actually updated the ndv->relsize here */
            extra = ndv->res_extra * g->xgrow * indivpergram;
            //printf("extra = %f\n,Species = %s \n", Species[sp]->name, extra);
            //printf("indivpergram = %f\n,Species = %s \n", Species[sp]->name, indivpergram);

            /* Check to make sure extra does not result in ndv->relsize > 1 (full-sized individual) */
            extra = GT(extra, 1 - ndv->relsize) ? 1 - ndv->relsize : extra;

            /* Sum extra growth of all individuals for each species */
            Species[sp]->extragrowth += extra;
            //printf("s->extragrowth  = %f\n, Species = %s \n", Species[sp]->name,s->extragrowth);

        } /*END ForEachIndiv */

        //printf("s->relsize before  = %f\n, Species = %s \n", Species[sp]->name, Species[sp]->relsize);
        //printf("s->extragrowth out of loop  = %f\n, Species = %s \n", Species[sp]->name,s->extragrowth);

        Species_Update_Newsize(sp, Species[sp]->extragrowth);
        //printf("s->relsize after = %f\n, Species = %s \n", Species[sp]->name, Species[sp]->relsize);
    } /* ENDFOR j (for each species)*/
}

/***********************************************************/
void rgroup_Establish(void)
{
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

	/* Cannot establish if plot is still in disturbed state*/
	if (Plot.disturbed > 0)
	{
		ForEachGroup(rg)
			RGroup[rg]->regen_ok = FALSE;
		return; /* skip regen for all */
	}

	ForEachGroup(rg)
	{
		g = RGroup[rg];
		if (!g->use_me)
			continue;

		g->regen_ok = TRUE; /* default */

		if (Globals.currYear < RGroup[rg]->startyr)
		{
			g->regen_ok = FALSE;

		}
		else  ///if ( g->max_age == 1 ) {
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

			ForEachGroupSpp(sp,rg,i)
			{
				if (!Species[sp]->use_me)
					continue;
				if (!Species[sp]->allow_growth)
					continue;
            //KAP:The below RandUni function has been added to now only establish annuals here, but allow establishment to be stochastic.
            //A random number is drawn from the uniform distribution and if that number is equal to or less than the probability of establishment,
            //then annual species will establish via the get_annual_maxestab function.
				if (Species[sp]->max_age == 1)
				{
				if (RandUni() <= Species[sp]->seedling_estab_prob)
				{
					num_est = _get_annual_maxestab(sp);
				}
				else
				{
				num_est = 0;
				}
					  //printf("num_est for annuals=%d \n",num_est);

					// above inserted to establish individuals for annuals
					// num_est for individuals is the number called from the seedbank in
					//     _get_annual_maxestab() (TEM 10-27-2015)
				}
				else
				{

					num_est = Species_NumEstablish(sp);
				}

				if (num_est)
				{
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
	numindvs = 0; // set to 0 so no problems passing to function
	/* compute the contribution of each indiv to the group's size */
	indivs = RGroup_GetIndivs(rg, SORT_0, &numindvs);
	for (n = 0; n < numindvs; n++)
		indivs[n]->grp_res_prop = indivs[n]->relsize / sumsize;
	//Mem_Free(indivs); // dont call free on variable unless it was initialized with malloc or calloc
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

	 HISTORY
	 Chris Bennett @ LTER-CSU 12/15/2000
	 *
	 *  2-Mar-03 - cwb - removed requirement for list to be
	 *  pre-allocated.  Added code to allocate and return list.
	 *------------------------------------------------------*/

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
