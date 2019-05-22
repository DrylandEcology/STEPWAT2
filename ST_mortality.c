/********************************************************/
/********************************************************/
/*  Source file: mortality.c
 *  Type: module
 *  Application: STEPPE - plant community dynamics simulator
 *  Purpose: This collection of routines implements all of
 *           the mortality functions.  mort_Main() is the
 *           entry point. All the other functions implement
 *           a specific mortality or are support routines. */
/*  History */
/*     (6/15/2000) -- INITIAL CODING - cwb */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "sw_src/filefuncs.h"
#include "generic.h"
#include "rands.h"
#include "myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/pcg/pcg_basic.h"

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
void rgroup_Extirpate( GrpIndex rg) ;
Bool indiv_Kill_Partial( MortalityType code,
                          IndivType *ndv,
                          RealF killamt);
void indiv_Kill_Complete( IndivType *ndv, int killType);
void check_sizes(const char *); /* found in main */
void _delete(IndivType *ndv);

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
void mort_Main( Bool *killed);
void mort_EndOfYear( void);
void proportion_Recovery(void);
void grazing_EndOfYear( void);
void rgroup_DropSpecies(SppIndex sp);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static void _pat( const SppIndex sp);
static void _mound( const SppIndex sp);
static void _burrow( const SppIndex sp);
static void _succulents( const SppIndex sp);
static void _slow_growth( const SppIndex sp);
static void _no_resources( GrpIndex rg);
static void _age_independent( const SppIndex sp);
static void _stretched_clonal( GrpIndex rg, Int start, Int last,
                           IndivType *nlist[]);
//Nov 4th 15 -AT - Made below two function non-static as they also getting call from main.c
void _kill_annuals(void);
void _kill_extra_growth(void);
void _kill_maxage(void);


/************ File-Level Variable Declarations *************/
/***********************************************************/
Bool *_SomeKillage;
/* flag: some plant was reduced and PR is affected. */
/* 7/5/01  - currently flag is set but unused. */
extern
  pcg32_random_t mortality_rng; //declared in ST_main.c


/***********************************************************/
/***********************************************************/
void mort_Main( Bool *killed) {
/*======================================================*/
/* PURPOSE */
/* This routine contains all the references to mortality
   that occurs during the growing season.  (See mort_EndOfYear()
   for other routines.)  It reduces "plant-space" according to
   resource availability, age, growth rate, and disturbance
   conditions.  Also, succulents are "mortified" if it's a wet
   year.

   An important consideration is the fact that the age of new
   plants starts at 1.  Consider that for it to establish, it
   has already survived the massive mortality of year 0; this
   rate of year 0 mortality, among other factors of survival,
   is captured by the probability of establishment (see
   rgroup_Establish()).  However, for consistency with other
   arrays, the age-related arrays are indexed from 0.  Thus,
   age is base1 but the arrays are base0.

   ALGORITHM
   The outline of the steps is:
     _no_resources(group)  - eqns 7, 8, 9
     _age_independent(spp) - eqn 14
     _slow_growth(spp) - slow growth rate, constant probability
     _succulents(sp) - if wet year
     Mort based on disturbances.

   More specifically:
   - Process each group separately.
   - If resources required for a group are greater than available
     (PR > 1) for the maximum years stretching is allowed, kill
     plants according the rules for insufficient resources.  If
     PR < 1, reset number of stretched years to zero.  This means
     that even one year of adequate resources nulls any number of
     years of stretched resources.
   - For each species, if the user requested the mortality functions
     to be executed, run the age-independent and slow-growth
     routines.
   - If the current species is a succulent and this is a wet year,
     execute the mortality function for succulents.  That is,
     reduce each plant by some constant proportion, defined by
     eqn 16.  The parameters are defined in the group-parms and
     the reduction amount is computed in the Env module based
     on precipitation.
   - execute disturbance effects, if any.
*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*   11/5/00 - moved NoResources() from growth routine
 *              to here.
 *   7-Nov-03 (cwb) No need to apply this mortality system
 *             to annuals.  That is now done in mort_EndOfYear()
 *             so the statistics can be accumulated. */
/*------------------------------------------------------*/

  Int j;
  GrpIndex rg;
  SppIndex sp;
  GroupType *g;

  *_SomeKillage = FALSE;

  ForEachGroup(rg) {
    g = RGroup[rg];
    if (g->est_count == 0) continue;
    /* annuals are not subject to these sources of mortality and instead die in _kill_annuals */
    if (g->max_age == 1) continue;  

    /* kill plants if low resources for consecutive years */
    /* increment yrs_neg_pr if pr > 1, else zero it. */
    /* one good year cancels all previous bad years. */
    if ( GT(g->pr, 1.0) ) {
       if (++g->yrs_neg_pr >= g->max_stretch)
          _no_resources( rg);

    } else {
      g->yrs_neg_pr = 0;
    }

    ForEachEstSpp(sp,rg,j) {

      /* Implement mortality types 1 and 2: age independent mortality and slow-growth 
      * mortality if rgroups are susceptible to those xres (use_mort) = 1 */
      if ( g->use_mort ) {
        _age_independent( sp );

        _slow_growth( sp );

      }
      /* Implement mortality of succulents if this year's PPT is above the PPT threshold that triggers succulent mortality */
      if (g->succulent
          && Env->wet_dry == Ppt_Wet
          && RandUni(&mortality_rng) <= Succulent->prob_death )
        _succulents( sp );

      /* Finally, implement mortality due to fecal pats, ant mounds, or animal burrows 
      * (disturbances originally conceptualized for the shortgrass steppe that result in 
      * plant mortality */
      switch (Plot->disturbance) {
        case FecalPat:
             _pat( sp);
             break;
        case AntMound:
             _mound( sp);
             break;
        case Burrow:
             _burrow( sp);
             break;
        case NoDisturb: // Added this case to prevent a compiler warning.
             break;
        case LastDisturb:// Added this case to prevent a compiler warning.
             break;
        default:
             break;
      }
      
    } /* end for each species*/

  } /* end ForEachGroup(rg) */

  killed = _SomeKillage;
}

/***********************************************************/
void mort_EndOfYear(void) {
    /*======================================================*/
    /* PURPOSE */
    /* Implements killing of plants through wildfire, prescribed fire, or another event
     * that would result in mortality in single year (killyr) or
     * extirpation (extirp) where plants are killed and do not return. */

    GrpIndex rg;
    GroupType *g = NULL;
    SppIndex sp;
    IntU j;
    RealF fire_possibility, random_number, biomass_cheatgrass;
    char *cheatgrass_name = "brte";
    int i = 0;
    Bool prescribed_fire_on = FALSE;
    
    /* Check species index number from the beginning to all the species in
     *  species.in , if the species name == checkname then get the biomass and stop the loop*/
    for (i = 0; i < Globals->sppCount; i++) { /* if species name = checkname = brte then get the biomass of brte(cheatgrass)*/
        if (strcmp(cheatgrass_name, Species[i]->name) == 0) {
            biomass_cheatgrass = Species_GetBiomass(i); /* calculate biomass of cheatgrass*/
            g = RGroup[Species[i]->res_grp];
            break;
        }
    }
    
    /* Set a random number outside of the loop to make sure the kill probability for each functional group is the same */
    random_number = RandUni(&mortality_rng);

    //determine if prescribed fire is on for any group. If TRUE, we do NOT want to simulate cheatgrass wildfire.
    ForEachGroup(rg){
      if(RGroup[rg]->killfreq > 0){
        prescribed_fire_on = TRUE;
        break;
      }
    }

    // in this for loop "g" refers to the RGroup of cheatgrass. RGroup[rg] refers
    // to the current iteration's RGroup.
    ForEachGroup(rg) {
      if (Globals->currYear < RGroup[rg]->startyr) {
        continue;
      }

      RGroup[rg]->prescribedfire = 0;
      RGroup[rg]->wildfire = 0;

      // If these conditions are true we want to simulate wildfire based on cheatgrass abundance
      if(!prescribed_fire_on && g != NULL){
        /* ------------------------- WILDFIRE BASED ON CHEATGRASS BIOMASS------------------------- */
        // Calculate fire_possibility
        if (g->ignition == 0) { 
          /* If ignition == 0, no wildfire occurs */
          fire_possibility = 0;
        } else if (biomass_cheatgrass < g->ignition) {
          /* If cheatgrass biomass is less than the biomass required for wildfire ignition, wildfire probability is very low*/
          fire_possibility = .01;
        } else { 
          /* Otherwise a wildfire probability is calculated, which increases with cheatgrass biomass*/
          fire_possibility = g->cheatgrass_coefficient + g->wild_fire_slope * biomass_cheatgrass;

          // Cap fire_possibility at 1. This isn't needed from an algorithmic perspective,
          // but a value greater than 1 does not make sense as a probability.
          if(fire_possibility > 1){
            fire_possibility = 1;
          }
        }

        // If a wildfire occurs this year
        if (random_number <= fire_possibility) {
          RGroup[rg]->killyr = Globals->currYear;
          /* Increase the number of wildfires that have occurred across all iterations in this year by 1 */
          RGroup[rg]->wildfire = 1;
        }
        /* ------------------------- END WILDFIRE BASED ON CHEATGRASS BIOMASS ------------------------- */

      } else if(Globals->currYear >= RGroup[rg]->killfreq_startyr) { // Otherwise simulate prescribed fire

        if(RGroup[rg]->killfreq < 1){
          /* --------------------- STOCHASTIC PRESCRIBED FIRE -------------------- */
          if(random_number <= RGroup[rg]->killfreq) {
            RGroup[rg]->killyr = Globals->currYear;

           /* Increase the number of prescribed fires that have occurred across all iterations in this year by 1 */
            RGroup[rg]->prescribedfire = 1;
          } 
          /* ------------------- END STOCHASTIC PRESCRIBED FIRE ----------------- */

        } else if (((Globals->currYear - RGroup[rg]->killfreq_startyr) % (IntU) RGroup[rg]->killfreq) == 0) {
          /* ------------------------ PRESCRIBED FIRE AT A FIXED RETURN INTERVAL ----------------------- */
          RGroup[rg]->killyr = Globals->currYear;
          /* Calculate the prescribed fire counts */
          RGroup[rg]->prescribedfire = 1;
          /* ------------------------ END PRESCRIBED FIRE AT A FIXED RETURN INTERVAL ------------------- */
        }
      }
      
      /* Kill all individuals of the functional group and don't let them re-establish */
      if (Globals->currYear == RGroup[rg]->extirp) {
          rgroup_Extirpate(rg);
          
      /* If the current year is a kill year, implement mortality */    
      } else if (Globals->currYear == RGroup[rg]->killyr) {
          RGroup_Kill(rg);
      }

      /* If the current year is a fire year, then remove extra_growth here
      instead of in _kill_extra_growth called in ST_main.c. Otherwise, 
      biomass will be non-zero in a fire year with complete killing */
      if (Globals->currYear == RGroup[rg]->killyr) {
        if (!RGroup[rg]->use_extra_res){
          continue;
        }

        ForEachGroupSpp(sp, rg, j) {
          /* If the species is turned off, continue */
          if (!Species[sp]->use_me){
            continue;
          }

          if (ZRO(Species[sp]->extragrowth)) continue;

          Species[sp]->extragrowth = LE(Species[sp]->extragrowth, Species[sp]->relsize) ? Species[sp]->extragrowth : Species[sp]->relsize;
          Species_Update_Newsize(sp, -Species[sp]->extragrowth);
          Species[sp]->extragrowth = 0.0;
        }
      }
    }  
}

void grazing_EndOfYear( void){

	/*======================================================*/
    /* PURPOSE */
    /* Implements grazing at the frequency and intensity that is specified in rgroup.in */
    /* HISTORY */
	/* 1st Nov 2015 -AT  -Added Species grazing EndOfYear  */
	/*======================================================*/
	
	GrpIndex rg;
	GroupType *g;

	ForEachGroup(rg)
	{
		IntU grazingyr =0;
		g = RGroup[rg];

		//printf("inside grazing_EndOfYear() year=%d, rgroupName=%s, grazingfreq_startyr=%d, grazingfreq=%d, proportionGrazing=%f startYear=%d \n",Globals->currYear,g->name,g->grazingfreq_startyr,g->grazingfrq,g->proportion_grazing, RGroup[rg]->startyr);

		if (Globals->currYear < RGroup[rg]->startyr)
		{
			/* Grazing cannot occur for an RGroup[rg] until the year that RGroup[rg] is turned on */
			continue;
		}

		if ((Globals->currYear >=g->grazingfreq_startyr) && (g->grazingfrq > 0))
		{
			if (g->grazingfrq < 1.0)
			{
				if (RandUni(&mortality_rng) <= g->grazingfrq)
				{
					grazingyr = Globals->currYear;
				}

			}
			else if (((Globals->currYear - g->grazingfreq_startyr) % (IntU) g->grazingfrq) == 0)
			{
				grazingyr = Globals->currYear;
			}

		}

		//Implement grazing if this year is a year where grazing should occur
		if (Globals->currYear == grazingyr)
		{
			//printf( "currYear is equal to grazingYear so will iterate all the Species for doing grazing, RGroup[g]->est_count =%d \n",RGroup[rg]->est_count);
			Int i;
			ForEachEstSpp2( rg, i)
			{
				if (!Species[RGroup[rg]->est_spp[i]]->use_me)
				{
					continue;
				}
				//printf( "year=%d calling Species_Proportion_Grazing()  rgroupName=%s, est_count =%d,grazingfreq_startyr=%d, grazingfreq=%d, proportionGrazing=%f \n",Globals->currYear,g->name,RGroup[rg]->est_count,g->grazingfreq_startyr,g->grazingfrq,g->proportion_grazing);
				
				/* Remove plant biomass to implement grazing using the proportion_grazing specified in inputs */
				Species_Proportion_Grazing(RGroup[rg]->est_spp[i],RGroup[rg]->proportion_grazing );
			}
		}
	}
}

void proportion_Recovery(void) {
    /*======================================================*/
    /* PURPOSE */
    /* Implements recovery of biomass that represents re-sprouting after a fire. This
     * is controlled by proportion_recovered, specified in inputs and can be turned
     * on or off for each functional group, depending on their capacity to resprout. */
    /* HISTORY */
    /* 1st Nov 2015 -AT -Added Species Proportion Recovery  */
    /*======================================================*/

    GrpIndex rg;
    SppIndex sp;

    ForEachGroup(rg) {

        if (Globals->currYear < RGroup[rg]->startyr) {
            /* Recovery of biomass after fire cannot occur for an RGroup[rg] until the year 
            * that RGroup[rg] is turned on */
            continue;
        }
		
        // Implement recovery of biomass after fire that represents re-sprouting
        if (Globals->currYear == RGroup[rg]->killyr) {
            Int i;

            //printf("'before proportion_recovery': Group = %s, relsize = %f, est_count = %d\n",
            //RGroup[rg]->name, RGroup[rg]->relsize, RGroup[rg]->est_count);

            ForEachEstSpp(sp, rg, i) {

                /* Annuals have already been killed in _kill_annuals and are not
                 * subject to proportion recovery after fire */
                if (Species[sp]->max_age == 1)
                    continue;

                //printf("'before proportion_recovery': Species = %s, relsize = %f, est_count = %d\n",
                // Species[sp]->name, Species[sp]->relsize, Species[sp]->est_count);

                Species_Proportion_Recovery(RGroup[rg]->est_spp[i], 6,
                        RGroup[rg]->proportion_recovered,
                        RGroup[rg]->proportion_killed);

                //printf("'after proportion_recovery': Species = %s, relsize = %f, est_count = %d\n",
                //Species[sp]->name, Species[sp]->relsize, Species[sp]->est_count);
            }
            //printf("'after proportion_recovery': Group = %s, relsize = %f, est_count = %d\n",
            // RGroup[rg]->name, RGroup[rg]->relsize, RGroup[rg]->est_count);
        }
    }
}

/***********************************************************/
static void _pat( const SppIndex sp) {
/*======================================================*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/
    Int i, k=-1;
    IndivType *p, **kills;
    
    kills = (IndivType **)Mem_Calloc(Globals.max_indivs_per_spp, sizeof(IndivType *), "_pat");

    /* ---------------------------------------------*/
    /* Generate kill list, depending on sensitivity */
    /* ---------------------------------------------*/
    if ( Plot->pat_removed) {
      /* get list of seedlings and annuals*/
      ForEachIndiv(p, Species[sp]) {
        if ( p->age == 1 || Species[sp]->disturbclass == VerySensitive)
          kills[++k] = p;
      }
      for( i=0; i<= k; i++)
        indiv_Kill_Complete( kills[i], 1);

    } else { /* kill according to disturbance class*/
      switch ( Species[sp]->disturbclass) {
        case VerySensitive:
        case Sensitive:
             Species_Kill(sp,1);
             k=1;
             break;
        case Insensitive:
        case VeryInsensitive:
          /* unaffected*/
      default:
            break;
      }
    }

    if (k >= 0) *_SomeKillage = TRUE;
    
    Mem_Free(kills);
}


/***********************************************************/
static void _mound( const SppIndex sp) {
/*======================================================*/
/* Ant mounds kill all but the hardiest plants. In C&L-90
 * that would be the succulents.
 */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

    Bool k = FALSE;

    switch ( Species[sp]->disturbclass) {
      case VerySensitive:
      case Sensitive:
      case Insensitive:
           Species_Kill(sp,2);
           k = TRUE;
           break;
      case VeryInsensitive:
        /* unaffected*/
    default:
         break;
    }


    if (k) *_SomeKillage = TRUE;
}


/***********************************************************/
static void _burrow( const SppIndex sp) {
/*======================================================*/
/* Kills all individuals on the plot if a burrow occurs.
 */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

    Bool k=FALSE;


    switch ( Species[sp]->disturbclass) {
      case VerySensitive:
      case Sensitive:
      case Insensitive:
      case VeryInsensitive:
           Species_Kill(sp,3);
           k = TRUE;
    }

    if (k ) *_SomeKillage = TRUE;
}


/***********************************************************/
static void _succulents( const SppIndex sp) {
/*======================================================*/

/* HISTORY */
/*   Chris Bennett @ LTER-CSU 6/15/2000            */
/*   cwb - 2-Dec-02 -- While adding SOILWAT code I found
 *      an old bug, ie while looping through the list of
 *      individuals, indiv_Kill_Partial could actually
 *      kill the whole plant, thus causing it to be deleted
 *      from the list, which in turn caused an error in
 *      the ForEachIndiv() code.  Now, indiv_Kill_Partial()
 *      returns FALSE if the amount to kill is greater
 *      than the size of the plant.  This routine was
 *      modified to make a new list of the dead plants
 *      and remove them properly. */
/*------------------------------------------------------*/

  IndivType *p,
            **kills;
  RealF killamt = Succulent->reduction;
  int i, k=0;
  
  kills = (IndivType **)Mem_Calloc(Globals.max_indivs_per_spp, sizeof(IndivType *), "_succulents");

  ForEachIndiv (p, Species[sp]) {
    if ( GT(p->relsize, killamt) )
      indiv_Kill_Partial( Slow, p, killamt);
    else
      kills[k++] = p;
  }

  for (i=0; i < k; i++)
    indiv_Kill_Complete(kills[i], 7);


  if (Species[sp]->est_count) *_SomeKillage = TRUE;
  
  Mem_Free(kills);
}


/***********************************************************/
static void _slow_growth( const SppIndex sp) {
/*======================================================*/
/* Kill plants based on a probability if the growth rate
   is less than the "slow rate" which is defined by the
   user in the group-level parameters (max_slow) and in
   the species-level parameters (max_rate). The slow rate
   is growthrate <= max_slow * max_rate.

   Increment the counter for number of years of slow growth.
   If the number of years of slow growth is greater than
   max_slow (defined in species.in), draw a random number
   and test it against the probability of mortality.  C&L'90
   defines this value as a constant, but it might be better
   to define it in the groups or species parameters.

   Of course, annuals aren't subject to this mortality,
   nor are new plants. */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

  Int n, k=-1;
  RealF pm = 0.368, /* probability of mortality*/
        slowrate;
  IndivType *ndv,
            **kills;
  
  kills = (IndivType **)Mem_Calloc(Globals.max_indivs_per_spp, sizeof(IndivType *), "_slow_growth");

  slowrate = RGroup[Species[sp]->res_grp]->slowrate
           * Species[sp]->max_rate;

  ForEachIndiv (ndv, Species[sp]) {
    if ( ndv->age == 1) continue;
    if (ndv->growthrate <= slowrate) {
      ndv->slow_yrs++;
      /* add to kill list if pm met*/
      if ( ndv->slow_yrs >= Species[sp]->max_slow
           && RandUni(&mortality_rng) <= pm)
         kills[++k] = ndv;
    } else
      ndv->slow_yrs = max( ndv->slow_yrs -1, 0);

  }

  for( n=0; n <= k; n++ )
    indiv_Kill_Complete(kills[n], 8);

  if (k >= 0) *_SomeKillage = TRUE;
  
  Mem_Free(kills);
}

/***********************************************************/
static void _age_independent( const SppIndex sp) {
/*======================================================*/

/* Kills individuals in a species by the age-independent function (eqn 14) in C&L'90
   assuming that AGEMAX was defined. */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/* 5/22/01 (cwb) - Annuals are not killed here. Slso, skip species with max_age==0 (longest lived). */
/*------------------------------------------------------*/
  
  Int n, k=-1;

  RealF pn, /* probability of mortality by year n (eqn 14)*/
       a;
  /* need a kill list because the bookkeeping in */
  /* Indiv_Kill() would confound kill-as-you-go*/
  IndivType **kills,
            *ndv;
  /*----------------------------------------------------*/

  assert(SppMaxAge(sp) != 0);

  if (SppMaxAge(sp) == 1) return;

  kills = (IndivType **) Mem_Calloc(Species[sp]->est_count,
                                   sizeof(IndivType *),
                                   "_age_independent(kills)");

  ForEachIndiv (ndv, Species[sp]) {
    a = (RealF)ndv->age / SppMaxAge(sp);
    pn = pow(SppMaxAge(sp), a -1)        /* EQN 14 */
         - (a * Species[sp]->cohort_surv);
    /* add to kill list if pn met*/
    if (RandUni(&mortality_rng) <= pn)
      kills[++k] = ndv;
  }

  for( n=0; n <= k; n++ ) {
    indiv_Kill_Complete(kills[n], 9);
  }

  if (k >= 0) *_SomeKillage = TRUE;

  Mem_Free(kills);

}

/***********************************************************/
static void _no_resources( GrpIndex rg) {
/*======================================================*/
/* use EQN 7, 8, 9 prior to growing
 * Resource limitation results in plant mortality
 * (ie, individuals or portions of individuals). In reality, this
 * would happen gradually over the season, but in this
 * yearly-time-step model it has to be done either before
 * or after the growth routine.  C&L 1990 (p241) note that
 * growth rates are also reduced--this is done in the main
 * growth loop by setting gmod = 1/PR.

 * Note to make it here, the rgroup's PR MUST BE > 1.0., which is checked in mort_Main().

 * This routine also calls _stretched_clonal() to kill
 * additional amounts of clonal plants (if any), which also
 * happens due to insufficient resources. */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000
 *     8/20/01  - replaced some loops with RGroup_GetIndivs()
 *         because that's why the function was created.
*/
/*------------------------------------------------------*/

  IntS i,
       n,   /* number of individuals in group */
       nk;  /* number of plants to kill */

  /* to-be-sorted list of indivs*/
  IndivType **indv_list;


  /*----------------------------------------------------*/
  /* get a sorted list of individuals in this rgroup*/
  indv_list = RGroup_GetIndivs(rg, SORT_A, &n);

  /*----------------------------------------------------*/
  /* kill until nk reached   (EQN 7)    */
  nk = (IntS) ( (n * (1.0 - 1.0/RGroup[rg]->pr)));
  for( i=0; i < nk; i++)
    indiv_Kill_Complete(indv_list[i], 10);

  if (nk) *_SomeKillage = TRUE;

  /* Check to see if this group's resources have been stretched,
   * and commit mortality of clonal plants which get additional
   * killamt if they are not killed in the preceeding loop.
   * Pick up next largest individual from where we left off in
   * previous loop by reusing i without resetting it to 0.
   * i comes out of the loop pointing to the next living plant.
   * If for some reason all plants were killed, _stretched_clonal
   * exits before doing anything. */
  _stretched_clonal( rg, i, n-1, indv_list);

  Mem_Free( indv_list);

}

/**************************************************************/
static void _stretched_clonal( GrpIndex rg, Int start, Int last,
                           IndivType *nlist[]) {
/*======================================================*/
/* Kill portions of clonal individuals when resources are limited */
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

  Int i,
      y,  /* number of years of stretched resources*/
      np, /* number of clonal plants in this resource group*/
      nk; /* number of clonal plants to kill if pm met*/
  RealF pm; /* Probability of mortality (eqn 8)*/

  /* These are used if reducing proportionally (pm not met)*/
  RealF total_size,
       indiv_size,
       total_reduction,
       indiv_reduction;

  IndivType **clist; /* list of clonal individuals */

  clist = (IndivType **)Mem_Calloc(Globals.max_indivs_per_spp, sizeof(IndivType *), "_stretched_clonal");
  
  /* get a list of remaining clonal plants, still ranked by size */
  for( np=-1, i=start; i <= last; i++) {
    if (Species[nlist[i]->myspecies]->isclonal)
      clist[++np] = nlist[i];
  }
  if (np < 0)
  {
    Mem_Free(clist);
    return;  /* Exit if no clonals remain alive in this rgroup */
  }
    
  y = RGroup[rg]->yrs_neg_pr;

  if (y >= RGroup[rg]->max_stretch) {
    pm = .04 * y * y;  /* EQN 8 from Coffin and Lauenroth (1990) */

    if (RandUni(&mortality_rng) <= pm ) {  /* kill on quota basis */
      /* must be more than 10 plants for any to survive ? */
      /* if so, then use ceil(), otherwise, use floor() */
      nk = (Int) floor(((RealF) (np+1) * 0.9)); /* EQN 9 from Coffin and Lauenroth (1990) */

      /* Kill until we reach quota or number of plants*/
      nk = min( nk, (np+1));
      for( i = 0; i < nk; i++) {
        indiv_Kill_Complete(clist[i], 11);
      }

      if (nk >= 0) *_SomeKillage = TRUE;

    } else {  /* reduce inverse-proportionally */

      total_reduction = 1.0 / RGroup[rg]->pr;

      /* Making sure PR will always be > 1 here */
      if (total_reduction > 1.0)
        LogError(logfp, LOGFATAL,
            "PR too large in Mort_StretchClonal()\n");

      /* sum up relsizes for total size*/
      for(i = 0, total_size = 0.0;
          i <= np;
          total_size += clist[i++]->relsize);

      /* the 0.8 is a "magic" correction to prevent the
         reduction from being too large; the assumption
         is that the plants are hardy enough to survive
         without dying beyond what is necessary to make
         resources required exactly <= availability. */
      total_reduction *= 0.8;
      for( i=0; i<= np; i++ ) {
        indiv_size = clist[i]->relsize / total_size;
        indiv_reduction = indiv_size * total_reduction;
        
        /* always succeeds if magic number < 1.0 */
        indiv_Kill_Partial( NoResources,
                            clist[i],
                            indiv_reduction);

      }
      if (np >= 0) *_SomeKillage = TRUE;

    } /* end if pm*/
  } /* end if y >= 1*/
  
  Mem_Free(clist);
}

/***********************************************************/
void _kill_annuals( void) {
/*======================================================*/
/* PURPOSE */
/* Loop through all species and kill the annual species.  This
   routine should be called at the end of the year after
   all growth happens and statistics are calculated and
   we don't need to know about the annuals any more.

   The assumption, of course, is that all of the annual
   species that are established are indeed one year old.
   See the discussion at the top of this file and in
   indiv_create() for more details. */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 3/14/2001 */
/* New function that kills all annual individuals (TEM 10-27-2015)) */
/*------------------------------------------------------*/

  GrpIndex rg;
  SppIndex sp;
  Int i;

  ForEachGroup(rg) {
    if (RGroup[rg]->max_age == 1) {
      for(i=RGroup[rg]->est_count, sp=RGroup[rg]->est_spp[i-1]; i>0; sp=RGroup[rg]->est_spp[(--i) - 1]){
               Species_Annual_Kill(sp, 4);             
          }
      }
    }

}

/***********************************************************/
void _kill_extra_growth(void) {
    /*======================================================*
     * PURPOSE *
     * Remove superfluous growth due to extra resources accumulated during the 
     * growing season. This should be done after all the statistics are accumulated 
     * for the year. 
     * HISTORY *
     * Updated by KAP 5/2018 */
    /*------------------------------------------------------*/
    
    IntU j;
    GrpIndex rg;
    SppIndex sp;
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
                ( (sizeof(x) == sizeof(float)) \
                                ? ((x)>-xF_DELTA && (x)<xF_DELTA) \
                                                : ((x)>-xD_DELTA && (x)<xD_DELTA) )

    ForEachGroup(rg) {

        ForEachGroupSpp(sp, rg, j) {
            /* If the species is turned off, continue */
            if (!Species[sp]->use_me)
                continue;

            /* Annuals have already been killed and relsize for annuals set to 
             * 0.0 in _kill_annuals */
            if (Species[sp]->max_age == 1)
                continue;
            //printf("s->extragrowth kill before  = %f\n", Species[sp]->extragrowth);

            /* Check that extragrowth <= s->relsize, otherwise relsize will become 
             * negative. If not, then reset to s->relsize. If the current year 
             * is a fire year, return, as killing of extragrowth has already occurred in Mort_EndofYear */
            if (!ZERO(Species[sp]->extragrowth) && Globals->currYear != RGroup[rg]->killyr) {
                Species[sp]->extragrowth = LE(Species[sp]->extragrowth, Species[sp]->relsize) ? Species[sp]->extragrowth : Species[sp]->relsize;

            // Sets relsize to 0 then sums up all the the individuals' relsizes and adds them back.
            // This has two effects:
            // 1: removes extragrowth from relsize by recalculating without it.
            // 2: Removes any small differences between the sum of individual relsizes and the species
            //    relsize. This ensures that Plot_Initialize in ST_main.c will always set relsize to 0. 
            IndivType *p = Species[sp]->IndvHead;
            Species[sp]->relsize = 0;
            while(p)
            {
                Species[sp]->relsize += p->relsize;
                p = p->Next;
            }
            RGroup_Update_Newsize(Species[sp]->res_grp); //need to call this to sync species and rgroup.

                //printf("s->relsize kill after  = %f\n", Species[sp]->relsize);
                Species[sp]->extragrowth = 0.0;
            }
            
            /* Now FINALLY remove individuals that were killed because of fire or grazing and set 
             * relsizes to 0, and remove the Species if the following cases are true */
            if (ZERO(Species[sp]->relsize) || LT(Species[sp]->relsize, 0.0)) {
                // printf("s->relsize in _kill_extra_growth check1 before = %f\n", Species[sp]->relsize);
                Species[sp]->relsize = 0.0;
                RGroup[rg]->relsize = 0.0;
                // printf("s->relsize in _kill_extra_growth check1 after = %f\n", Species[sp]->relsize);

                IndivType *p1 = Species[sp]->IndvHead, *t1;
                while (p1) {
                    t1 = p1->Next;
                    _delete(p1);
                    p1 = t1;
                }
                rgroup_DropSpecies(sp);
            }
        }
    }
#undef xF_DELTA
#undef xD_DELTA
#undef ZERO    
}

/******************************************************************************/
void _kill_maxage(void) {
/******************************************************************************/
/* PURPOSE:
 * Kill plants once they reach their maximum age.
 *
 * HISTORY:
 * Created by Frederick Pierson on 4/6/2019. */
/******************************************************************************/
    
    SppIndex s;
    IndivType *i;
    
    ForEachSpecies(s) {
        ForEachIndiv(i, Species[s]) {
            if (i->age == Species[s]->max_age) {
                indiv_Kill_Complete(i, 12);
            }
        }
    }
}
