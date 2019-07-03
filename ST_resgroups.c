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

extern
  pcg32_random_t resgroups_rng;

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
void copy_rgroup(const GroupType* src, GroupType* dest);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static void _res_part_extra(RealF extra, RealF size[]);
static GroupType *_create(void);
static void _extra_growth(GrpIndex rg);
static void _add_annual_seedprod(SppIndex sp, RealF lastyear_relsize);
static RealF _get_annual_maxestab(SppIndex sp);
static RealF _add_annuals(const GrpIndex rg, const SppIndex sp, const RealF lastyear_relsize);

/****************** Begin Function Code ********************/

/***********************************************************/

void rgroup_PartResources(void) {
    /* Partition resources for this year among the resource
    groups. The allocation happens in three steps: basic
    group allocation, partitioning of extra resources in _res_part_extra(),
    and further partitioning to individuals in the group _ResPartIndiv(). */

    GrpIndex rg;
    Bool noplants = TRUE;
    GroupType *g; /* shorthand for RGroup[rg] */

    /* ----- distribute basic (normal) resources */
    ForEachGroup(rg) {
        g = RGroup[rg];
        
        //Set resources available and resources required for each functional group
        g->res_required = RGroup_GetBiomass(rg);
        g->res_avail = SXW_GetTranspiration(rg);
        //printf("g->res_avail = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_avail);
        //printf("g->res_required = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_required);

        /* A check, which should not happen. Res_required is set to the number of
        established individuals for the species */
        if (ZRO(g->res_avail) && g->res_required > 0) {
            //printf("'rgroup_PartResources': Group = %s(%d) error with res (relsize = %f, pr = %f): \n",
            //    g->name, rg, g->relsize, g->pr);
            //printf("\tbefore correction: res_avail = %f, res_required = %f\n",
            //    g->res_avail, g->res_required);

            g->res_required = 0.0;
            Int i;
            ForEachEstSpp2(rg, i) {
                g->res_required += Species[g->est_spp[i]]->est_count;
            }

            g->res_avail = 1;

           // printf("\tafter correction: res_avail = %f, res_required = %f\n",
           //     g->res_avail, g->res_required);

            LogError(logfp, LOGWARN, "RGroup %s : res_avail is Zero and res_required > 0", g->name);
        }

        /* Calculate PR at the functional group level: resources required/resources available */
        g->pr = ZRO(g->res_avail) ? 0. : g->res_required / g->res_avail;
        //printf("g->pr = %f\n,Group = %s \n",RGroup[rg]->name,  g->pr);

        /* If relsize>0 and individuals are established, reset noplants from TRUE to FALSE */
        if (GT(getRGroupRelsize(rg), 0.))
            noplants = FALSE;

    } /* End ForEachGroup(rg) */

    /* If noplants=TRUE, exit from the loop */
    if (noplants)
        return;

    //Partition resources to individuals and determine 'extra' resources
    rgroup_ResPartIndiv();
}

static RealF _add_annuals(const GrpIndex rg, const SppIndex sp, const RealF lastyear_relsize) {
    /*======================================================*/
    /* Establishment for annual species, which includes a function that modifies
     * the seedbank through seed addition and deletion, a function that determines
     * the amount of viable seed in the seedbank if regen_ok is true, a random draw
     * from the beta distribution that is used to determine the number of seeds
     * that will emerge as seedlings this year, and removal of those germinated
     * seeds from the seedbank.*/

    IntU i, num_est; //number of individuals that will establish this year
    RealF viable_seeds; //number of viable seeds in the seedbank
    float var; //random number draw from beta distribution for calculation of num_est
    GroupType *g;
    SpeciesType *s;

    g = RGroup[rg];
    s = Species[sp];

    /*Increment viable years for seeds and implement the decay process (seed mortality).
    * Then add seeds to the seedbank*/
    _add_annual_seedprod(sp, lastyear_relsize);

    /*Get viable seeds from seed bank*/
    viable_seeds = (g->regen_ok) ? _get_annual_maxestab(sp) : 0;

    /* Create a beta random number draw based on alpha and beta for each species
     * (calculated based on mean (s->seedling_estab_prob and variance (s->var)) */
     var = RandBeta(Species[sp]->alpha, Species[sp]->beta, &resgroups_rng);

    /*Determine number of seedlings to add. If the number of seeds calculated
     * from the random draw is larger than max_seed_estab, max_seed_estab is used instead*/
    if (g->extirpated)
    {
        num_est = 0;
    }
    else
    {
        num_est = min(viable_seeds * var, s->max_seed_estab);
        //printf("Species name=%s , num_est   =%u \n",s->name,  num_est);
    }


    /*Multiple the proportion of seeds in each viable year array by the total
    number of seeds that germinated as seedlings and subtract those seeds from
    the relevant seedprod array.*/
       for (i = 0; i < s->viable_yrs; i++) {
        //printf("Species name=%s , old calculated value s->seedprod[%hu]= %d \n", s->name, i, s->seedprod[i]);
        s->seedprod[i] =  s->seedprod[i] -  round(num_est * s->seedprod[i] / viable_seeds);
        //printf("Species name=%s , so new calculated value s->seedprod[%hu]= %d \n", s->name, i, s->seedprod[i]);
       }
        return num_est;
}

static RealF _get_annual_maxestab(SppIndex sp) {
    /*======================================================*/
    /* Get the maximum number of viable seeds from the seedbank that can
 establish this year*/
    IntU i;
    RealF sum = 0.; //sum of the viable seedbank
    SpeciesType *s = Species[sp];

    for (i = 0; i < s->viable_yrs; i++) {
        sum += s->seedprod[i];
    }
    //printf("sum =%f \n",sum);

    return sum;
}

static void _add_annual_seedprod(SppIndex sp, RealF lastyear_relsize) {
    /*======================================================*/
    /* Increase the age of seeds and then add seeds to the seedbank.*/

    SpeciesType *s = Species[sp];
    IntU i;

    /*Age of seeds is increased by one year, seed mortality occurs, and seeds produced in the
    previous year are added to the seedbank at relative age 0 */
    for (i = s->viable_yrs - 1; i > 0; i--) {
        //printf("Species name=%s , seedprod before decay s->seedprod[%hu]= %d \n", s->name, i, s->seedprod[i]);
        s->seedprod[i] = s->seedprod[i - 1] / pow(i, s->exp_decay);
        //printf("Species name=%s , seedprod after decay s->seedprod[%hu]= %d \n", s->name, i, s->seedprod[i]);
    }

    // printf("Species name=%s ,old Array array 0 index i=%u, value =%hu \n", s->name, i, s->seedprod[i]);

    /* If the current year is year 1 of the simulation, then the number of seeds
     * added is a random number draw between 1 and the maximum number of seedlings
     * that can establish. Otherwise, this year's seed production is a function
     * of the number of seeds produced per unit biomass multiplied by species biomass
     * (maximum species biomass * last year's species relative size). */
    if (Globals->currYear == 1) {
        s->seedprod[0] = RandUniIntRange(1, s->max_seed_estab, &resgroups_rng);
        //printf("Species name=%s ,currYear =1 so new calculated value s->seedprod[%u]= %hu , s->max_seed_estab =%hu\n", s->name, i, s->seedprod[i], s->max_seed_estab);

    } else {
        s->seedprod[0] = (IntU) (s->pseed * s->mature_biomass * lastyear_relsize);
        //printf("Species name=%s ,currYear=%hu  so new calculated value s->seedprod[%u]= %hu , s->max_seed_estab =%hu, lastyear_relsize=%.5f\n", s->name, Globals->currYear, i, s->seedprod[i], s->max_seed_estab, lastyear_relsize);
    }
}

/***********************************************************/
static void _res_part_extra(RealF extra, RealF size[]) {
    /*======================================================*
     * PURPOSE *
     * Partitions "extra" resources to other groups that can use them (if any). 
     * HISTORY *
     * Updated by KAP 5/2018 */

    GrpIndex rg;
    GroupType *g; /* shorthand for RGroup[rg] */
    RealF req_prop, /* group's prop'l contrib to the total requirements */
            sum_size = 0.; /* summed size of all functional groups that can use extra resources */
    
    /* Determine the summed size of functional groups that can use extra resources */
    ForEachGroup(rg) {
        g = RGroup[rg];
        if (ZRO(getRGroupRelsize(rg)))
            continue;
        if (g->est_count == 0)
            continue;

        // Where size = size_obase, which is biomass(g/m2) or 0 if the group can't use resources
        sum_size += size[rg];
        //printf("size[rg]  = %f\n,Rgroup = %s \n",RGroup[rg]->name, size[rg]);
        //printf("sum_size  = %f\n,Rgroup = %s \n",RGroup[rg]->name, sum_size);

    } /* End ForEachGroup(rg) */

    //printf("sum_size  = %f\n", sum_size);

    ForEachGroup(rg) {
        g = RGroup[rg];
        if (ZRO(getRGroupRelsize(rg)))
            continue;
        if (g->est_count == 0)
            continue;

        /* Check to avoid dividing by 0 */
        if (sum_size == 0.)
            req_prop = 0.;

        /* Calculate proportional biomass of each group out of the total biomass
         * of all functional groups that can use extra resources */
        else
            req_prop = size[rg] / sum_size;

        /* If the group can use extra resources, divide out extra based on
         * proportional biomass */
        if (g->use_extra_res)
            g->res_extra = req_prop * extra;

        /* If the group can't use extra resources, set res_extra to 0 */
        else
            g->res_extra = 0.;
        //printf("res_extra = %f\n,Rgroup = %s \n",RGroup[rg]->name,g->res_extra);

    } /* End ForEachGroup(rg) */

}
/***********************************************************/
void rgroup_ResPartIndiv(void) {
    /* Resources at the functional group level are divided among individuals in 
     * the group, largest to smallest. Species distinctions within the group are 
     * ignored. The partitioning of resources to individuals is done proportionally 
     * based on size so that large individuals get more resources than their 
     * smaller competitors.
     * HISTORY:
     * Chris Bennett @ LTER-CSU 12/21/2000
     * Updated by KAP 5/2018 */

    GrpIndex rg;
    GroupType *g; /* shorthand for RGroup[rg] */
    SppIndex sp;
    Int j;
    IndivType **indivs, /* dynamic array of indivs in RGroup[rg] */
            *ndv; /* shorthand for the current indiv */
    IntS numindvs, n;
    RealF base_rem = 0., /* remainder of resource after allocating to an indiv */
            xtra_obase = 0., /* summed extra resources across all groups */
            *size_base, /* biomass of the functional group */
    		*size_obase; /* biomass of functional groups that can use extra resources */

    size_base = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF), "rgroup_ResPartIndiv");
    size_obase = (RealF *)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(RealF), "rgroup_ResPartIndiv");
    
    /* Divide each group's normal resources to individuals */
    ForEachGroup(rg) {
        g = RGroup[rg];
        if (g->est_count == 0)
            continue;

        /* Allocate the temporary group-oriented arrays */
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

        /* Sum "extra" resources for all functional groups */
        xtra_obase += base_rem;
        //printf("xtra_obase = %f\n", xtra_obase);

        /* Check to see if each functional group can use extra resources. If so,
         * then pass the biomass of the FG, otherwise pass 0 into _res_part_extra */
        size_base[rg] = RGroup_GetBiomass(rg);
        size_obase[rg] = (g->use_extra_res) ? size_base[rg] : 0.;
        //printf("size_obase = %f\n", size_obase[rg]);

        Mem_Free(indivs);

    } /* end ForEachGroup() */

    /* Assign extra resources to functional groups that can use them */
    _res_part_extra(xtra_obase, size_obase);

    /* Now loop back through all individuals in each group, assign extra resource
     * and calculate PR at the individual level. This extra resource will be used to
     * increment plant sizes. Extra resources that remain after this step which are 
     * applied to superfluous biomass increment that will ultimately be killed at the 
     * end of the year is assigned in _extra_growth */
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

        /* Allocate the temporary group-oriented arrays */
        indivs = RGroup_GetIndivs(rg, SORT_D, &numindvs);

        ForEachEstSpp(sp, rg, j) {

            /* Compute PR at the individual level and assign extra resources */
            for (n = 0; n < numindvs; n++) {
                ndv = indivs[n];

                    //printf("ndv->res_avail before  = %f\n", ndv->res_avail);

                    /* If individuals already have the resources they require do 
                     * not assign extra */
                    if (ndv->res_avail == ndv->res_required) {
                        ndv->res_extra = 0.0;
                        //printf("ndv->res_extra = %f\n", ndv->res_extra);
                    }
                    
                    /* Assign extra resource as the difference of what is required
                    * by the individual using what was assigned to the individual 
                    * through partitioning of normal resources */
                    else {
                        ndv->res_extra = fmin(ndv->res_required - ndv->res_avail, g->res_extra);
                        //printf("ndv->res_required = %f\n", ndv->res_required);
                        //printf("ndv->res_avail = %f\n", ndv->res_avail);
                        //printf("ndv->res_extra = %f\n", ndv->res_extra);

                        /* Updated resources available for each individual to include extra */
                        ndv->res_avail += ndv->res_extra;
                        //printf("ndv->res_avail after = %f\n", ndv->res_avail);

                        /* Remaining extra resources to be used for superfluous growth */
                        g->res_extra = fmax(g->res_extra - ndv->res_extra, 0.);
                        //printf("g->res_extra in loop = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_extra);
                    }

                /* Calculate the PR value, or dflt to 100. Used in growth module */
                ndv->pr = GT(ndv->res_avail, 0.) ? ndv->res_required / ndv->res_avail : 100.;
                //printf("ndv->pr  = %f\n", ndv->pr);
            }
        } /* end ForEachEstSpp() */

        //printf("g->res_extra after = %f\n,Group = %s \n",RGroup[rg]->name,  g->res_extra);

        Mem_Free(indivs);

    } /* end ForEachGroup() */
    
    Mem_Free(size_obase);
    Mem_Free(size_base);
}

/***********************************************************/
void rgroup_Grow(void) {
    /*======================================================*
     * PURPOSE *
     * Main loop to allow for plant growth.
     * HISTORY *
     * Chris Bennett @ LTER-CSU 6/15/2000
     *  7-Nov-03 (cwb) Annuals have a completely new method of
     *     propogation, growth, and death.  Most notably, they
     *     aren't subject to the growth mechanism here, so they
     *     are now excluded from this code.  All of the annual
     *     generation occurs in the PartResources() function.
     *  7/13/2016 KAP: This comment by cwb is now outdated. We had added
     * functionality so that annuals now grow on an annual basis.
     * Updated by KAP 5/2018 */

    IntU j;
    GrpIndex rg;
    SppIndex sp;
    GroupType *g;
    SpeciesType *s;
    const RealF OPT_SLOPE = .05; /* from Coffin and Lauenroth 1990, EQN 5 */
    RealF growth1, /* growth of one individual */
            sppgrowth, /* sum of growth for a species' indivs */
            rate1, /* rate of growth for an individual */
            tgmod, /* temperature growth factor modifier */
            gmod; /* growth factor modifier */
    IndivType *ndv; /* temp pointer for current indiv */

    ForEachGroup(rg) {
        g = RGroup[rg];

        if (g->est_count == 0)
            continue;

        /* Succulents don't grow if conditions are wet */
        if (g->succulent && Env->wet_dry == Ppt_Wet)
            continue;

        /* Increment size of each individual in response to normal resources */
        ForEachEstSpp(sp, rg, j) {
            s = Species[sp];

            sppgrowth = 0.0;
            if (!Species[sp]->allow_growth)
                continue;

            /* Modify growth rate by temperature calculated in Env_Generate() */
            if (s->tempclass != NoSeason)
                tgmod = Env->temp_reduction[s->tempclass];

            /* Now increase size of the individual plants of current species */
            ForEachIndiv(ndv, s) {
                /* Growth rate (gmod) initially set to 0.95. Values for gmod range
                between 0.05 and 0.95 similar to Coffin and Lauenroth 1990 */
                gmod = 1.0 - OPT_SLOPE;
                
                /* Reduction in gmod if PR > 1 (resource limitation) */
                if (GT(ndv->pr, 1.0))
                    gmod /= ndv->pr;
                
                /* Further modification of gmod based on this year's temperature */
                gmod *= tgmod;
                
                /* For clonal species, if the individual appears killed it was 
                 * reduced due to low resources last year. It can reproduce vegetatively 
                 * this year but couldn't last year. */
                if (ndv->killed && RandUni(&resgroups_rng) < ndv->prob_veggrow) {
                    growth1 = s->relseedlingsize * RandUniIntRange(1, s->max_vegunits, &resgroups_rng);
                    rate1 = growth1 / ndv->relsize;
                    ndv->killed = FALSE;
                    //printf("growth1 killed  = %f\n", growth1);
                } else {
                    /* Normal growth: modifier times optimal growth rate (EQN 1)
                     * in Coffin and Lauenroth 1990 */
                    rate1 = gmod * s->intrin_rate * (1.0 - ndv->relsize);
                    growth1 = rate1 * ndv->relsize;
                }
                //printf("growth1  = %f\n", growth1);
                //printf("old ndv->relsize  = %f\n,Species = %s \n", Species[sp]->name, ndv->relsize);

                ndv->relsize += growth1;
                //printf("new ndv->relsize  = %f\n, Species = %s \n", Species[sp]->name,ndv->relsize);

		// Save the increment in size due to normal resources for use in the grazing module.
		ndv->normal_growth = growth1;

                ndv->growthrate = rate1;

                sppgrowth += growth1;
                //printf("sppgrowth = %f\n, Species = %s \n", Species[sp]->name,sppgrowth);

            } /*END ForEachIndiv */

        } /* ENDFOR j (for each species)*/
        
        /* Implement growth due to extra resources which will support superfluous 
         * increases in biomass that are ultimately removed at the end of the year */
        _extra_growth(rg);

    } /* END ForEachGroup(rg)*/
}

/***********************************************************/
static void _extra_growth(GrpIndex rg) {
    /*======================================================*/
    /* PURPOSE */
    /* When there are resources beyond the minimum necessary for "normal" growth, 
     * the extra resources are converted to superfluous growth in the current year 
     * and is removed at the end of the year in _kill_extra_growth. */
    /* HISTORY */
    /* Chris Bennett @ LTER-CSU 11/8/2000 */
    /* Updated by KAP 5/2018 */

    Int j;
    RealF extra_ndv, indivpergram;
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
            extra_ndv = 0.0;

            RGroup_Update_GrpResProp(rg);
            /* Calculate the extra resource available to each individual based on size */
            ndv->res_extra = ndv->grp_res_prop * g->res_extra;
            //printf("ndv->res_extra = %f\n,Species = %s \n", Species[sp]->name, ndv->res_extra);
            //printf("ndv->grp_res_prop = %f\n,Species = %s \n", Species[sp]->name, ndv->grp_res_prop);

            /* Calculate the increment in size due to res_extra for each individual
            Note, we have not actually updated the ndv->relsize here */
            extra_ndv = ndv->res_extra * g->xgrow * indivpergram;
            //printf("extra = %f\n,Species = %s \n", Species[sp]->name, extra);
            //printf("indivpergram = %f\n,Species = %s \n", Species[sp]->name, indivpergram);

            /* Check to make sure extra does not result in ndv->relsize > 1 (full-sized individual) */
            extra_ndv = GT(extra_ndv, 1 - ndv->relsize) ? 1 - ndv->relsize : extra_ndv;

            /* Sum extra growth of all individuals for each species */
            Species[sp]->extragrowth += extra_ndv;
            //printf("s->extragrowth  = %f\n, Species = %s \n", Species[sp]->name,s->extragrowth);

        } /*END ForEachIndiv */

    } /* ENDFOR j (for each species)*/
}


/***********************************************************/
void rgroup_Establish(void) {
    /*======================================================*/
    /* PURPOSE */
    /* Determines which and how many species can establish in a given year. For
     * each species in each perennial functional group, check that a uniform
     * random number between 0 and 1 is less than the species' establishment
     * probability. a) If so, return a random number of individuals up to the
     * maximum allowed to establish for the species. This is the number of individuals
     * in this species that will establish this year. b) If not, continue with
     * the next species. Establishment for species of annual functional groups
     * occurs differently. See notes at the top of _add_annuals */

    /* HISTORY */
    /* Chris Bennett @ LTER-CSU 6/15/2000 */
    /* The probability of establishment emulates the occurrence of microsite
     * conditions that allow for establishment.
     * 7-Nov-03 (cwb) Adding the new algorithm to handle annuals. It's more
     * complicated than before (which didn't really work) so annuals are now added
     * in the PartResources()function.  Only perennials are added here. Also, there's
     * now a parameter to define the start year of establishment for perennials.
     * KAP: Annual establishment now occurs here instead of in PartResources*/
    /*------------------------------------------------------*/

    IntS i, num_est; /* number of individuals of sp. that establish*/
    GrpIndex rg;
    SppIndex sp;
    GroupType *g;

    /* Cannot establish if plot is still in disturbed state*/
    if (Plot->disturbed > 0) {
        ForEachGroup(rg)
        RGroup[rg]->regen_ok = FALSE;
        return; /* skip regen for all */
    }

    /* If the functional group is turned off, continue */
    ForEachGroup(rg) {
        g = RGroup[rg];
        if (!g->use_me)
            continue;

        g->regen_ok = TRUE; /* default */

        if (Globals->currYear < RGroup[rg]->startyr) {
            g->regen_ok = FALSE;

        } else {

            ForEachGroupSpp(sp, rg, i) {

                /* If the species is turned off, continue */
                if (!Species[sp]->use_me)
                    continue;
                if (!Species[sp]->allow_growth)
                    continue;

                /* Establishment for species that belong to annual functional groups*/
                if (Species[sp]->max_age == 1) {
                    //printf("Globals->currYear = %hu, call to _add_annuals sp=%d Species[sp]->lastyear_relsize : %.5f \n", Globals->currYear, sp, Species[sp]->lastyear_relsize);
                    num_est = _add_annuals(rg, sp, Species[sp]->lastyear_relsize);
                    //  printf("g->seedbank annuals=%d \n",g->seedbank);
                }

                    /* Establishment for species that belong to perennial functional groups*/
                else {
                    num_est = Species_NumEstablish(sp);
                }

                if (num_est) {
                    /* printf("%d %d %d %d\n",
                     Globals->currIter, Globals->currYear, sp, num_est); */

                    Species_Add_Indiv(sp, num_est);
                    species_Update_Estabs(sp, num_est);
                }

            } /* end ForEachGroupSpp() */
        }
    } /* end ForEachGroup() */
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
							Species[ndv->myspecies]->max_age, Globals->currIter,
							Globals->currYear);
				}
			}
		}
	}
}

/************************************************************/
/* Sums relsize for all individuals in all species in RGroup rg.
   param rg = RGroup index.
   Return: RGroup relsize. */
RealF getRGroupRelsize(GrpIndex rg){
    Int n;
	SppIndex sp;
	double sum = 0.0;

    ForEachEstSpp( sp, rg, n){
		sum += getSpeciesRelsize(sp);
    }

    if(RGroup[rg]->est_count > 0){
        return (RealF) sum / (RealF) RGroup[rg]->est_count;
    } else {
        return 0;
    }
}

/***********************************************************/
/* Update the grp_res_prop field of every individual in the RGroup. 
   param rg = RGroup index */
void RGroup_Update_GrpResProp(GrpIndex rg)
{
	Int n;
	IndivType **indivs;
	IntS numindvs;
	RealF sumsize = 0.0;

	/* first get group's relative size adjusted for num indivs */
	/* ie, groupsize=1 when 1 indiv of each species is present */
	/* ie each indiv is an equivalent contributor, not based on biomass */
	sumsize = getRGroupRelsize(rg);

    /*printf("Group = %s, sumsize = %f, est_count = %d, relsize = %f\n",
        RGroup[rg]->name, sumsize, RGroup[rg]->est_count, RGroup[rg]->relsize);
    */

	numindvs = 0; // set to 0 so no problems passing to function
	/* compute the contribution of each indiv to the group's size */
	indivs = RGroup_GetIndivs(rg, SORT_0, &numindvs);
	for (n = 0; n < numindvs; n++)
		indivs[n]->grp_res_prop = indivs[n]->relsize / sumsize;

	Mem_Free(indivs);

	/* double check some assumptions */
	if (RGroup[rg]->est_count < 0)
		RGroup[rg]->est_count = 0;
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
		biomass += getSpeciesRelsize(sp) * Species[sp]->mature_biomass;
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

	p = (GroupType *) Mem_Calloc(1, sizeof(GroupType), "_create");
        p->name = (char *) Mem_Calloc(SuperGlobals.max_groupnamelen + 1, sizeof(char), "_create");
        p->est_spp = (SppIndex *) Mem_Calloc(SuperGlobals.max_spp_per_grp, sizeof(SppIndex), "_create");
        p->species = (SppIndex *) Mem_Calloc(SuperGlobals.max_spp_per_grp, sizeof(SppIndex), "_create");
        
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
	GrpIndex i = (GrpIndex) Globals->grpCount;

	if (++Globals->grpCount > SuperGlobals.max_rgroups)
	{
		LogError(logfp, LOGFATAL, "Too many groups specified (>%d)!\n"
				"You must adjust MAX_RGROUPS in maxrgroupspecies.in!",
		SuperGlobals.max_rgroups);
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
	}

	RGroup[rg]->extirpated = TRUE;

}

void copy_rgroup(const GroupType* src, GroupType* dest){
    int i; 
    SppIndex sp;

    // This would be very bad.
    if(src == dest){
        return;
    }

    /* -------------- Copy any arrays -------------- */
    Mem_Free(dest->kills);
    dest->kills = (IntUS*) Mem_Calloc(GrpMaxAge(src->grp_num), sizeof(IntUS), "copy_rgroup: kills");
    for(i = 0; i < GrpMaxAge(src->grp_num); ++i){
        dest->kills[i] = src->kills[i];
    }

    Mem_Free(dest->est_spp);
    dest->est_spp = (SppIndex*) Mem_Calloc(src->est_count, sizeof(SppIndex), "copy_rgroup: est_spp");
    for(i = 0; i < src->est_count; ++i){
        dest->est_spp[i] = src->est_spp[i];
    }

    /* ------------- Copy all fields --------------- */
    dest->cheatgrass_coefficient = src->cheatgrass_coefficient;
    dest->depth = src->depth;
    dest->est_annually = src->est_annually;
    dest->est_count = src->est_count;
    dest->estabs = src->estabs;
    dest->extirp = src->extirp;
    dest->extirpated = src->extirpated;
    dest->grazingfreq_startyr = src->grazingfreq_startyr;
    dest->grazingfrq = src->grazingfrq;
    dest->grp_num = src->grp_num;
    dest->ignition = src->ignition;
    dest->killfreq = src->killfreq;
    dest->killfreq_startyr = src->killfreq_startyr;
    dest->killyr= src->killyr;
    dest->max_age = src->max_age;
    dest->max_bmass = src->max_bmass;
    dest->max_density = src->max_density;
    dest->max_per_sqm = src->max_per_sqm;
    dest->max_spp = src->max_spp;
    dest->max_spp_estab = src->max_spp_estab;
    dest->max_stretch = src->max_stretch;
    dest->min_res_req = src->min_res_req;
    dest->mm_extra_res = src->mm_extra_res;
    strcpy(dest->name, src->name);
    dest->ppt_intcpt[0] = src->ppt_intcpt[0];
    dest->ppt_intcpt[1] = src->ppt_intcpt[1];
    dest->ppt_intcpt[2] = src->ppt_intcpt[2];
    dest->ppt_slope[0] = src->ppt_slope[0];
    dest->ppt_slope[1] = src->ppt_slope[1];
    dest->ppt_slope[2] = src->ppt_slope[2];
    dest->pr = src->pr;
    dest->prescribedfire = src->prescribedfire;
    dest->proportion_grazing = src->proportion_grazing;
    dest->proportion_killed = src->proportion_killed;
    dest->proportion_recovered = src->proportion_recovered;
    dest->regen_ok = src->regen_ok;
    dest->res_avail = src->res_avail;
    dest->res_extra = src->res_extra;
    dest->res_required = src->res_required;
    dest->rgroupFractionOfVegTypeBiomass = src->rgroupFractionOfVegTypeBiomass;
    dest->slowrate = src->slowrate;
    dest->startyr = src->startyr;
    dest->succulent = src->succulent;
    dest->use_extra_res = src->use_extra_res;
    dest->use_me = src->use_me;
    dest->use_mort = src->use_mort;
    dest->veg_prod_type = src->veg_prod_type;
    dest->wild_fire_slope = src->wild_fire_slope;
    dest->wildfire = src->wildfire;
    dest->xgrow = src->xgrow;
    dest->yrs_neg_pr = src->yrs_neg_pr;

    /* ---------------- Copy Species Array ----------------- */
    for(i = 0; i < SuperGlobals.max_spp_per_grp; ++i){
        dest->species[i] = src->species[i];
    }
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
