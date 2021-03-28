/**
 * \file ST_mortality.c
 * \brief Implements all of the mortality functions.  
 * 
 * mort_Main() is the entry point. All the other functions implement
 * a specific type of mortality or are support routines.
 *  
 * \author
 *     Chandler Haukap\n
 *     Kyle Palmquist\n
 *     Chris Bennett\n
 *     Ashish Tiwari
 * 
 * \date 23 August 2019
 * 
 * \ingroup MORTALITY
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */
#include <stdlib.h>
#include <string.h>

#include "ST_mortality.h"
#include "sw_src/filefuncs.h"
#include "sw_src/rands.h"
#include "sw_src/myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/pcg/pcg_basic.h"
#include "sxw_vars.h"

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

void rgroup_DropSpecies(SppIndex sp);

/* ---- Local function declarations. Treat these functions as private. ----- */
static void _pat( const SppIndex sp);
static void _mound( const SppIndex sp);
static void _burrow( const SppIndex sp);
static void _succulents( const SppIndex sp);
static void _slow_growth( const SppIndex sp);
static void _no_resources( GrpIndex rg);
static void _age_independent( const SppIndex sp);
static void _stretched_clonal( GrpIndex rg, Int start, Int last, 
                              IndivType *nlist[]);
void _updateCheatgrassPrecip(int year);
double _getCheatgrassCover(double biomass);
double _getWildfireProbability(double percentCover);
Bool _simulateWildfire(double cheatgrassCover);
Bool _simulatePrescribedFire(void);
double _getCheatgrassBiomass(void);

/********************** Private mortality objects ****************************/
/**
 * \brief [Mortality](\ref MORTALITY)'s private precipitation information.
 * 
 * \author Chandler Haukap
 * \date 13 January 2020
 * \ingroup MORTALITY_PRIVATE
 */
CheatgrassPrecip* cheatgrassPrecip = 0;

/**
 * \brief TRUE if a function in ST_mortality.c killed an individual.
 * \ingroup MORTALITY_PRIVATE
 */
Bool *_SomeKillage = 0;

/**
 * \brief Performs mortality that occurs during the growing season.
 * 
 *  This function reduces "plant-space" according to
 *  resource availability, age, growth rate, and fecal pats, 
 *  ant mounds, and animal burrows. Also, succulents are 
 *  killed if it's a wet year.
 * 
 *  This function is not to be confused with mort_EndOfYear().
 *
 *  An important consideration is the fact that the age of new
 *  plants starts at 1.  Consider that for it to establish, it
 *  has already survived the massive mortality of year 0; this
 *  rate of year 0 mortality, among other factors of survival,
 *  is captured by the probability of establishment (see
 *  rgroup_Establish()).  However, for consistency with other
 *  arrays, the age-related arrays are indexed from 0.  Thus,
 *  age is base1 but the arrays are base0.
 * 
 *  The outline for this function is:\n
 *    1) Death by limited resources using eqns 7, 8, 9.\n
 *    2) Age independent mortality using eqn 14.\n
 *    3) Slow growth mortality.\n
 *    4) Succulent wet season mortality.\n
 *    5) Mortality due to fecal pats, ant mounds or animal burrows.\n
 *  With all equations from Coffin and Lauenroth 1990.
 *
 *  More specifically:
 *  - Process each group separately.
 *  - If resources required for a group are greater than available
 *    (PR > 1) for the maximum years stretching is allowed, kill
 *    plants according the rules for insufficient resources.  If
 *    PR < 1, reset number of stretched years to zero.  This means
 *    that even one year of adequate resources nulls any number of
 *    years of stretched resources.
 *  - For each species, if the user requested the mortality functions
 *    to be executed, run the age-independent and slow-growth
 *    routines.
 *  - If the current species is a succulent and this is a wet year,
 *    execute the mortality function for succulents.  That is,
 *    reduce each plant by some constant proportion, defined by
 *    eqn 16.  The parameters are defined in the group-parms and
 *    the reduction amount is computed in the Env module based
 *    on precipitation.
 *  - Execute disturbances from fecal pats, ant mounds, and animal 
 *    burrows, if any.
 * 
 *  Note that other disturbances (i.e. fire and grazing) are implemented 
 *  elsewhere.
 * 
 * \param killed Reference to a boolean. Set to TRUE if one or more
 *               individuals were killed in any species. FALSE otherwise.
 * 
 * \sideeffect Individuals are killed. Individuals are stored in a
 *             linked list in Species[sp]->IndivHead, so expect the 
 *             list to be modified for all sp.\n
 *             _SomeKillage will be set to TRUE if any individuals die.
 * 
 * \sa mort_EndOfYear()
 * \sa rgroup_Establish()
 * 
 * \ingroup MORTALITY
 */
void mort_Main( Bool *killed) {

  Int j;
  GrpIndex rg;
  SppIndex sp;
  GroupType *g;

  *_SomeKillage = FALSE;

  ForEachGroup(rg) {
    g = RGroup[rg];
    if (g->est_count == 0) continue;
    /* annuals are not subject to these sources of mortality and instead die in killAnnuals */
    if (g->max_age == 1) continue;

    /* Calculate PR at the functional group level: resources required/resources available */
    g->pr = ZRO(g->res_avail) ? 0. : g->res_required / g->res_avail;

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

/**
 * \brief Simulates prescribed fire for all [groups](\ref RGROUP).
 * 
 * Fire is simulated based on the flags read in from inputs.
 * 
 * \return TRUE if prescribed fire is on for at least one [group](\ref RGROUP).
 *         FALSE if prescribed fire is off for all [groups](\ref RGROUP).
 *         Note that the return values does _not_ indicate whether the 
 *         _current_ year is a fire year.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup MORTALITY_PRIVATE
 */
Bool _simulatePrescribedFire(void) {
  GrpIndex rg;
  double randomNumber = RandUni(&mortality_rng);
  Bool prescribedFireOn = FALSE;

  ForEachGroup(rg) {
    if((Globals->currYear < RGroup[rg]->killfreq_startyr) || 
       (RGroup[rg]->killfreq_startyr == 0)) {
      continue;
    }
    
    if (RGroup[rg]->killfreq < 1) {
      prescribedFireOn = TRUE;
      /* -------------------- STOCHASTIC PRESCRIBED FIRE ------------------- */
      if (randomNumber <= RGroup[rg]->killfreq) {
        RGroup[rg]->killyr = Globals->currYear;
        /* Mark this year as a prescribed fire year for the statistics. */
        RGroup[rg]->prescribedfire = 1;
      }
      /* ------------------ END STOCHASTIC PRESCRIBED FIRE ----------------- */
    } else if (((Globals->currYear - RGroup[rg]->killfreq_startyr) % 
               (IntU) RGroup[rg]->killfreq) == 0) {
      prescribedFireOn = TRUE;
      /* ------------ PRESCRIBED FIRE AT A FIXED RETURN INTERVAL ----------- */
      RGroup[rg]->killyr = Globals->currYear;
      /* Mark this year as a prescribed fire year for the statistics. */
      RGroup[rg]->prescribedfire = 1;
      /* --------- END PRESCRIBED FIRE AT A FIXED RETURN INTERVAL ---------- */
    }
  }

  return prescribedFireOn;
}

#ifdef STDEBUG
Bool (*simulatePrescribedFire)(void) = _simulatePrescribedFire;
#endif

/***********************************************************/
/**
 * \brief Simulates fires and extirpation.
 * 
 * Implements killing of plants through wildfire, prescribed fire, or another event
 * that would result in mortality in single year (killyr) or
 * extirpation (extirp) where plants are killed and do not return.
 * 
 * If wildfire is simulated it will affect all individuals in all species in all rgroups.
 * If prescribed fire is simulated it is possible to specify inputs such that fire will 
 * only affect one rgroup.
 * 
 * \sideeffect If fire is simulated RGroup[rg]->killyr will be set to the current year.\n
 *             If fire is simulated all individuals in RGroup[rg] will be killed.\n
 *             If prescribed fire is simulated RGroup[rg]->prescribedfire will be set to 1.\n
 *             If wildfire is simulated RGroup[rg]->wildfire will be set to 1.\n
 *             If extirpation is performed all species in RGroup[rg] will be dropped.\n
 *             
 * \sa mort_Main()
 * \sa rgroup_Extirpate()
 * 
 * \ingroup MORTALITY
 */
void mort_EndOfYear(void) {
  GrpIndex rg;
  SppIndex sp;
  IntU j;
  RealF cheatgrassBiomass;

  /* Update the precipitation parameters that determine cheatgrass-driven
   * wildfire. This function call is commented out for the time being until the
   * struct that it updates (cheatgrassPrecip) can be incorporated into the 
   * rest of the code. */
  // _updateCheatgrassPrecip(Globals->currYear);
  /* printf("%d:\tMeanSpring = %f\t ThisSpring = %f\t MeanWinter = %f\t LastWinter = %f\n",
         Globals->currYear, cheatgrassPrecip->springMean, cheatgrassPrecip->currentSpring, 
         cheatgrassPrecip->winterMean, cheatgrassPrecip->lastWinter); */

  /* Reset prescribedfire and wildfire for every RGroup. */
  ForEachGroup(rg){
    RGroup[rg]->prescribedfire = 0;
    RGroup[rg]->wildfire = 0;
  }

  if(UseCheatgrassWildfire) {
    cheatgrassBiomass = _getCheatgrassBiomass();
    // If cheatgrass was found in the Species array.
    if(cheatgrassBiomass >= 0) {
      _simulateWildfire(cheatgrassBiomass);
    }
  } else {
    _simulatePrescribedFire();
  }

  // For all RGroups determine if this year is a kill year. If it was, implement
  // killing.
  ForEachGroup(rg) {
    /* Kill all individuals of the functional group and don't let them re-establish */
    if (Globals->currYear == RGroup[rg]->extirp) {
      rgroup_Extirpate(rg);
    } else if (Globals->currYear == RGroup[rg]->killyr) {
      RGroup_Kill(rg);
    }

    /* If the current year is a fire year, then remove extra_growth here
      instead of in killExtraGrowth called in ST_main.c. Otherwise,
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
        Species[sp]->extragrowth = 0.0;
      }
    }
  }
}

/**
 * \brief Implements grazing at the frequency and intensity that is specified in rgroup.in
 * 
 * This is a straightforward function. It checks if the current year is a grazing year for
 * every rgroup. If it is, it performs grazing on all species in the given group. 
 * 
 * Most of the time we want grazing to occur for all groups in the same year. However, this
 * function checks each RGroup individually, so any number of groups could be grazed in any
 * given year.
 * 
 * \sideeffect If RGroup[rg] is grazed every individual's biomass in every species in RGroup[rg]
 *            will be reduced based on RGroup[rg]->proportion_grazing.
 * 
 * \sa Species_Proportion_Grazing() this function is called to modify any species that need
 *                                  grazing.
 * 
 * \ingroup MORTALITY
 */
void grazing_EndOfYear( void){
	GrpIndex rg;
	GroupType *g;

	ForEachGroup(rg)
	{
		IntU grazingyr =0;
		g = RGroup[rg];

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
			Int i;
			ForEachEstSpp2( rg, i)
			{
				if (!Species[RGroup[rg]->est_spp[i]]->use_me)
				{
					continue;
				}
        
				/* Remove plant biomass to implement grazing using the proportion_grazing specified in inputs */
				Species_Proportion_Grazing(RGroup[rg]->est_spp[i],RGroup[rg]->proportion_grazing );
			}
		}
	}
}

/**
 * \brief Initialize the variables necessary for cheatgrass-driven wildfire.
 * 
 * This function can be called multiple times, and should be called before
 * every iteration of the program. It will allocate memory the first time it is
 * run, then every subsequent time it will reset the fields to 0.
 * 
 * \sideeffect
 *     This function will allocate memory the first time it is run.
 * 
 * \sa CheatgrassPrecip
 * \author Chandler Haukap
 * \date 13 January 2020
 * \ingroup MORTALITY
 */
void initCheatgrassPrecip(void) {
  /* If cheatgrassPrecip hasn't been allocated */
  if(!cheatgrassPrecip){
    cheatgrassPrecip = Mem_Calloc(1, sizeof(CheatgrassPrecip),
                                  "initCheatgrassPrecip: cheatgrassPrecip");
  }

  /* Reset all fields to 0 */
  cheatgrassPrecip->currentSpring = 0;
  cheatgrassPrecip->lastWinter = 0;
  cheatgrassPrecip->prevSprings[0] = 0;
  cheatgrassPrecip->prevSprings[1] = 0;
  cheatgrassPrecip->prevSprings[2] = 0;
  cheatgrassPrecip->currentSpring = 0;
  cheatgrassPrecip->springMean = 0;
  cheatgrassPrecip->winterMean = 0;
  cheatgrassPrecip->lastOctThruDec = 0;
  cheatgrassPrecip->thisOctThruDec = 0;
  cheatgrassPrecip->thisJanThruMar = 0;
}

/**
 * \brief Free all memory allocated inside the 
 *        [mortality module](\ref MORTALITY).
 * 
 * \sideeffect
 *     All local variables that have been dynamically allocated will be 
 *    deallocated.
 * 
 * \author Chandler Haukap
 * \date 14 January 2020
 * \ingroup MORTALITY
 */
void freeMortalityMemory(void) {
  if(cheatgrassPrecip){
    Mem_Free(cheatgrassPrecip);
  }
}

/**
 * \brief Load a different CheatgrassPrecip struct into the 
 *        [mortality](\ref MORTALITY) module.
 * 
 * This is particularly useful in [gridded mode](\ref GRID) where each 
 * [cell](\ref CellType) has it's own precipitation information.
 * 
 * As far as allocation goes, you can either allocate the new 
 * \ref CheatgrassPrecip struct yourself, or, call this function then
 * call \ref initCheatgrassPrecip.
 * 
 * \author Chandler Haukap
 * \date 13 January 2020
 * \ingroup MORTALITY
 */
void setCheatgrassPrecip(CheatgrassPrecip* newCheatgrassPrecip) {
  cheatgrassPrecip = newCheatgrassPrecip;
}

/**
 * \brief Returns a pointer to \ref cheatgrassPrecip.
 * 
 * This is necessary in [gridded mode](\ref GRID) in order to allocate more
 * than 1 \ref cheatgrassPrecip.
 * 
 * \sa setCheatgrassPrecip
 * \author Chandler Haukap
 * \date 14 January 2020
 * \ingroup MORTALITY
 */
CheatgrassPrecip* getCheatgrassPrecip(void) {
  return cheatgrassPrecip;
}

/**
 * \brief Updates \ref cheatgrassPrecip with the current year's precipitation.
 * 
 * Call this function before determining if a cheatgrass-driven wildfire should
 * occur for the given year.
 * 
 * Note that \ref cheatgrassPrecip must be allocated before calling this
 * function, either by calling \ref initCheatgrassPrecip or by allocating a new
 * \ref CheatgrassPrecip variable then calling \ref setCheatgrassPrecip.
 * 
 * \param year is the year in which this function is being called. It is
 *             necessary when calculating the running precipitation averages.
 * 
 * \sideeffect
 *     Every field in \ref cheatgrassPrecip will be updated to reflect this
 *     year's values.
 * 
 * \author Chandler Haukap
 * \date 13 January 2020
 * \ingroup MORTALITY_PRIVATE
 */
void _updateCheatgrassPrecip(int year) {
  int i;
  /* Note that arrays are 0 indexed in C. Therefore, when I say
   * SXW->ppt_monthly[0] I am refering to January's precipitation. */

  /* --------------------- Update the running averages --------------------- */
  if(year > 1) {
    // We only have a data point for the spring after year 1.
    cheatgrassPrecip->springMean = 
        get_running_mean(year - 1, cheatgrassPrecip->springMean, 
                         cheatgrassPrecip->currentSpring);

    if(year > 2) {
      // We only have a data point for the winter after year 2 because we need
      // The Oct - Dec values from 2 years previously
      cheatgrassPrecip->winterMean = 
          get_running_mean(year - 2, cheatgrassPrecip->winterMean,
                           cheatgrassPrecip->lastWinter);
    }
  }

  /* --------------- Shift the Spring values back one place ---------------- */
  cheatgrassPrecip->prevSprings[2] = cheatgrassPrecip->prevSprings[1];
  cheatgrassPrecip->prevSprings[1] = cheatgrassPrecip->prevSprings[0];
  cheatgrassPrecip->prevSprings[0] = cheatgrassPrecip->currentSpring;

  /* ----------- Calculate this year's mean Spring precipitation ----------- */
  cheatgrassPrecip->currentSpring = 0;
  for(i = 3; i < 6; ++i){
    cheatgrassPrecip->currentSpring += SXW->ppt_monthly[i];
  }
  cheatgrassPrecip->currentSpring /= 3;

  /* ------------- Calculate last year's Winter precipitation -------------- */
  cheatgrassPrecip->lastWinter = 0;
  cheatgrassPrecip->lastWinter += cheatgrassPrecip->lastOctThruDec;
  cheatgrassPrecip->lastWinter += cheatgrassPrecip->thisJanThruMar;
  cheatgrassPrecip->lastWinter /= 6;

  /* --------------- Shift the Winter values back one place ---------------- */
  // For the same reason we shift the Spring values back one position we need
  // to shift the winter values to keep them current.
  // October - December
  cheatgrassPrecip->lastOctThruDec = cheatgrassPrecip->thisOctThruDec;
  cheatgrassPrecip->thisOctThruDec = 0;
  for(i = 9; i < 12; ++i){
    cheatgrassPrecip->thisOctThruDec += SXW->ppt_monthly[i];
  }

  // January - March (this one isn't really a shift; its just a rewrite)
  cheatgrassPrecip->thisJanThruMar = 0;
  for(i = 0; i < 3; ++i){
    cheatgrassPrecip->thisJanThruMar += SXW->ppt_monthly[i];
  }
}

/**
 * \brief Recovers biomass that represents re-sprouting after a fire.
 * 
 * This is controlled by proportion_recovered, specified in inputs and can be turned
 * on or off for each functional group, depending on their capacity to resprout.
 * 
 * \sideeffect If this year was a fire year for RGroup[rg] every individual in every
 *             species in RGroup[rg] will be given a chance to recover some biomass. 
 * 
 * \sa Species_Proportion_Recovery() which is called by this function to determine
 *                                   how much biomass should be recovered.
 * 
 * \ingroup MORTALITY
 */
void proportion_Recovery(void) {
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

                /* Annuals have already been killed in killAnnuals and are not
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

/**
 * \brief Implements mortality by fecal pats.
 * 
 * This function is called my mort_Main(). If there is a pat on the plot
 * it kills any annual individuals and any sensitive individuals.
 * 
 * \param sp The index of the species to potentially kill.
 * 
 * \sideeffect If Species[sp]->disturbclass is VerySensitive or Sensitive 
 *             all individuals of the species will be killed.\n
 *             _SomeKillage is set to TRUE if at least one individual is
 *             killed.
 * 
 * \sa mort_Main() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _pat( const SppIndex sp) {
    Int i, k=-1;
    IndivType *p, **kills;
    
    kills = (IndivType **)Mem_Calloc(SuperGlobals.max_indivs_per_spp, sizeof(IndivType *), "_pat");

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


/**
 * \brief Implements mortality by ant mounds.
 * 
 * This function is called my mort_Main(). If there is a mound on the plot
 * it kills most species with the exception of succulents (Coffin and 
 * Lauenroth 1990).
 * 
 * \param sp The index of the species to potentially kill.
 * 
 * \sideeffect If Species[sp]->disturbclass is more sensitive than VeryInsensitive
 *             all individuals of the species will be killed.\n
 *             _SomeKillage is set to TRUE if at least one individual is
 *             killed.
 * 
 * \sa mort_Main() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _mound( const SppIndex sp) {
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


/**
 * \brief Implements mortality by animal burrows.
 * 
 * This function is called my mort_Main(). If there is a burrow on the plot
 * it kills all individuals in Species[sp]
 * 
 * \param sp The index of the species to potentially kill.
 * 
 * \sideeffect Species[sp]->disturbclass doesn't matter for burrows. All 
 *             individuals are killed.
 *             _SomeKillage is set to TRUE if at least one individual is
 *             killed.
 * 
 * \sa mort_Main() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _burrow( const SppIndex sp) {
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


/**
 * \brief Reduces succulent biomass based on reduction specified in inputs 
 *        due to precipitation.
 * 
 * This function will remove the same amount of biomass from all individuals 
 * in Species[sp]. If the removal removes all biomass the individual is killed.
 * 
 * \param sp the index in Species of the succulent.
 * 
 * \sideeffect Some individuals of Species[sp] will most likely be killed, and
 *             all individuals will loose some biomass.\n
 *             _SomeKillage set to TRUE if some individuals were killed.
 * 
 * \sa indiv_Kill_Partial(), which is called to reduce biomass (but not kill).
 * \sa indiv_Kill_Complete(), which is called to kill individuals.
 * \sa mort_Main() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _succulents( const SppIndex sp) {
  IndivType *p,
            **kills;
  RealF killamt = Succulent->reduction;
  int i, k=0;
  
  kills = (IndivType **)Mem_Calloc(SuperGlobals.max_indivs_per_spp, sizeof(IndivType *), "_succulents");

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


/**
 * \brief Implements slow growth mortality.
 * 
 * Kill plants if the growth rate is less than the "slow rate" which is 
 * defined by the user in the group-level parameters (max_slow) and in
 * the species-level parameters (max_rate). The slow rate is 
 * growthrate <= max_slow * max_rate.
 *
 * Increment the counter for number of years of slow growth.
 * If the number of years of slow growth is greater than
 * max_slow (defined in species.in), draw a random number
 * and compare it to the probability of mortality.  C&L'90
 * defines this value as a constant, but it might be better
 * to define it in the groups or species parameters.
 *
 * Of course, annuals aren't subject to this mortality,
 * nor are 1 year old plants.
 * 
 * \param sp The index of the species in the Species array.
 * 
 * \sideeffect Calculates ndv->slow_yrs for all individuals in
 *                  the species.\n 
 *              If ndv->slow_yrs passes the mortality threshold
 *                  the given individual is killed.\n
 *              _SomeKillage is set to TRUE is at least on individual
 *                  is killed.
 * 
 * \sa mort_Main() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _slow_growth( const SppIndex sp) {
  Int n, k=-1;
  RealF pm = 0.368, /* probability of mortality*/
        slowrate;
  IndivType *ndv,
            **kills;
  
  kills = (IndivType **)Mem_Calloc(SuperGlobals.max_indivs_per_spp, sizeof(IndivType *), "_slow_growth");

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

/**
 * \brief Implements age-independent mortality.
 * 
 * Kills individuals in a species by the age-independent function (eqn 14) in 
 * C&L'90 assuming that AGEMAX was defined.
 * 
 * Annuals are NOT killed here. Also, if Species[sp]->max_age == 0 (a species 
 * can live indefinitely) this function does nothing.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000  
 * 
 * \param sp The index of the species in the Species array.
 * 
 * \sideeffect Some individuals from the Species[sp]->IndivHead linked list 
 *             might be removed if this function kills them.\n
 *             _SomeKillage will be set to TRUE is an individual is killed.
 * 
 * \sa mort_Main() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _age_independent( const SppIndex sp) {
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

/**
 * \brief Kills individuals to simulate death by insufficient resources.
 * 
 * Note that this function does NOT check if resource availibility is low.
 * It kills an amount of individuals proportional to resource limitations,
 * but deciding if any individuals should be killed at all (whether or
 * not to call this function) takes place in mort_Main().
 * 
 * In reality, resource death would occur gradually over the growing
 * season. However, in this yearly-time-step model, it must be done either
 * before or after the growing season. 
 * 
 * This function takes an additional step with clonal species. It calls
 * _stretched_clonal() to kill additional amounts of clonal plants due 
 * to insufficient resources.
 * 
 * Note that to make it here, the rgroup's PR MUST BE > 1.0., which is 
 * checked in \ref mort_Main().
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \param rg The Index in RGroup of the group to kill.
 * 
 * \sideeffect Individuals killed by this function will be removed from
 *             the corresponding individual linked list.\n
 *             _SomeKillage set to TRUE if an individual was killed.
 * 
 * \sa _stretched_clonal() which is called to perfom additional mortality.
 * \sa mort_Main() where _no_resources() is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _no_resources( GrpIndex rg) {
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

/**
 * \brief Kill portions of clonal individuals when resources are limited.
 * 
 * Mortality is based on equations 8 and 9 of Coffin and Lauenroth (1990).
 * 
 * \param rg The group to perform mortality on.
 * \param start Index in nlist to start considering killing individuals.
 * \param last Index in nlist to stop considering killing.
 * \param nlist An array of IndivType pointers that are under consideration.
 * 
 * \sideeffect This function will kill individuals from nlist if there
 *             are any clonal individuals.\n
 *             _SomeKillage set to TRUE if an individual is killed.
 * 
 * \sa _no_resources() where this function is called.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
static void _stretched_clonal( GrpIndex rg, Int start, Int last,
                           IndivType *nlist[]) {

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

  clist = (IndivType **)Mem_Calloc(SuperGlobals.max_indivs_per_spp, sizeof(IndivType *), "_stretched_clonal");
  
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

/**
 * \brief Kill all individuals belonging to annual species
 * 
 * Loop through all species and kill the annual species.  This
 * routine should be called at the end of the year after
 * all growth happens and statistics are calculated and
 * we don't need any more information about the annuals.
 *
 * The assumption, of course, is that all of the annual
 * species that are established are indeed one year old.
 * See the discussion at the top of this file and in
 * indiv_create() for more details. 
 * 
 * \sideeffect All individuals in all annual species will
 *             be deleted.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
void killAnnuals( void) {
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

/**
 * \brief Remove superfluous growth due to extra resources accumulated during the 
 * growing season. 
 * 
 * This should be done after all the statistics are accumulated for the year. 
 * 
 * Updated by KAP 5/2018.
 * 
 * \sideeffect Species[sp]->extragrowth will be set to 0 for all sp.\n
 *             If this results in 0 biomass for a species it will also be
 *             dropped from the RGroup.
 * 
 * \ingroup MORTALITY_PRIVATE
 */
void killExtraGrowth(void) {
    IntU j;
    GrpIndex rg;
    SppIndex sp;

    ForEachGroup(rg) {

        ForEachGroupSpp(sp, rg, j) {
            /* If the species is turned off, continue */
            if (!Species[sp]->use_me)
                continue;

            /* Extra growth might have been set to zero in Mort_EndofYear if it is a fire year.
             * However, setting it to zero again is just as fast as checking. */
            Species[sp]->extragrowth = 0.0;

            /* Now FINALLY remove individuals that were killed because of fire or grazing and set 
             * relsizes to 0, and remove the Species if the following cases are true */
            if (getSpeciesRelsize(sp) <= 0.0) {
                // printf("s->relsize in killExtraGrowth check1 before = %f\n", Species[sp]->relsize);
                // printf("s->relsize in killExtraGrowth check1 after = %f\n", Species[sp]->relsize);

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
}

/**
 * \brief Kill plants once they reach their maximum age.
 * 
 * Created by Frederick Pierson on 4/6/2019.
 * 
 * \sideeffect Every individual in every species that has reached it's
 *             max age will be deleted (killed).
 * 
 * \ingroup MORTALITY_PRIVATE
 */
void killMaxage(void) {
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

/**
 * \brief Converts the biomass of cheatgrass to the % cover of cheatgrass.
 * 
 * This relationship between biomass and percent cover was derived by Maggie 
 * England.
 * 
 * \param biomass is the biomass of cheatgrass.
 * 
 * \return A double 0 and 100 representing the percent cover of cheatgrass.
 * 
 * \author Maggie England (derived the algorithm)
 * \author Chandler Haukap (implemented the code)
 * \date February 5 2020
 * \ingroup MORTALITY_PRIVATE
 */
double _getCheatgrassCover(double biomass) {
  double cover = biomass / (10.296 * Globals->plotsize);
  return (cover > 100) ? 100 : cover;
}

/**
 * \brief Calculates the probability of a cheatgrass-driven wildfire occuring.
 * 
 * This equation was derived by Maggie England.
 * 
 * \param percentCover is the percent of the total plot covered in cheatgrass.
 *                     A value between 0 and 100 is expected. I suggest using 
 *                     the \ref _getCheatgrassCover function.
 * 
 * \return A double between 0 and 1 representing the probability of a wildfire.
 * 
 * \author Maggie England (derived the equation)
 * \author Chandler Haukap (implemented the code)
 * \date February 5 2020
 * \ingroup MORTALITY_PRIVATE
 */
double _getWildfireProbability(double percentCover) {
  return 0.015 * pow(percentCover, 0.0649);
}

/**
 * \brief Returns the biomass of cheatgrass.
 * 
 * This function assumes that cheatgrass is present in the simulation and named
 * "brte".
 * 
 * \return The biomass of cheatgrass. 
 * \return -1 if cheatgrass cannot be found.
 * 
 * \author Chandler Haukap
 * \ingroup MORTALITY_PRIVATE
 */
double _getCheatgrassBiomass() {
  char *cheatgrassName = "brte";
  SppIndex sp;

  ForEachSpecies(sp) {
    /* if this species is cheatgrass then get the biomass */
    if (strcmp(cheatgrassName, Species[sp]->name) == 0) {
      return Species_GetBiomass(sp);
    }
  }

  return -1;
}

/**
 * \brief Simulates cheatgrass-driven wildfire.
 * 
 * \param cheatgrassBiomass the biomass of cheatgrass.
 * 
 * \return TRUE if a wildfire happens.
 * \return FALSE if no wildfire happens.
 * 
 * \author Chandler Haukap
 * \date February 6 2020
 * \ingroup MORTALITY_PRIVATE
 */
Bool _simulateWildfire(double cheatgrassBiomass) {
  GrpIndex rg;
  double percentCover = _getCheatgrassCover(cheatgrassBiomass);
  Bool wildfire = FALSE;
  
  if(RandUni(&mortality_rng) < _getWildfireProbability(percentCover)) {
    wildfire = TRUE;
    ForEachGroup(rg) {
      RGroup[rg]->wildfire = 1;
      RGroup[rg]->killyr = Globals->currYear;
    }
  }

  return wildfire;
}
