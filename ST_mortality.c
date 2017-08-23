/********************************************************/
/********************************************************/
/*  Source file: mortality.c
/*  Type: module
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: This collection of routines implements all of
 *           the mortality functions.  mort_Main() is the
 *           entry point. All the other functions implement
 *           a specific mortality or are support routines.
/*  History:
/*     (6/15/2000) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "generic.h"
#include "rands.h"
#include "myMemory.h"
#include "ST_steppe.h"
#include "ST_globals.h"

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
void rgroup_Extirpate( GrpIndex rg) ;
Bool indiv_Kill_Partial( MortalityType code,
                          IndivType *ndv,
                          RealF killamt);
void indiv_Kill_Complete( IndivType *ndv, int killType);
void check_sizes(const char *); /* found in main */

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
void mort_Main( Bool *killed);
void mort_EndOfYear( void);
void proportion_Recovery(void);
void grazing_EndOfYear( void);

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


/************ File-Level Variable Declarations *************/
/***********************************************************/
Bool _SomeKillage;
/* flag: some plant was reduced and PR is affected. */
/* 7/5/01  - currently flag is set but unused. */


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
               to here.

 *   7-Nov-03 (cwb) No need to apply this mortality system
 *             to annuals.  That is now done in mort_EndOfYear()
 *             so the statistics can be accumulated.

/*------------------------------------------------------*/

  Int j;
  GrpIndex rg;
  SppIndex sp;
  GroupType *g;

  _SomeKillage = FALSE;

  ForEachGroup(rg) {
    g = RGroup[rg];
    if (!g->est_count) continue;
    if (g->max_age == 1) continue;  /* annuals get theirs in EndOfYear() */

    /* mortify plants if low resources for consecutive years */
    /* increment yrs_neg_pr if pr > 1, else zero it. */
    /* one good year cancels all previous bad years. */
    if ( GT(g->pr, 1.0) ) {
       if (++g->yrs_neg_pr >= g->max_stretch)
          _no_resources( rg);
    } else {
      g->yrs_neg_pr = 0;
    }

    ForEachEstSpp(sp,rg,j) {

      /* Take care of mortality types 1 and 2*/
      if ( g->use_mort ) {
        _age_independent( sp );
        _slow_growth( sp );
      }
      /* now deal with succulents problems*/
      if (g->succulent
          && Env.wet_dry == Ppt_Wet
          && RandUni() <= Succulent.prob_death )
        _succulents( sp );

      /* finally, implement disturbance mortality*/
      switch (Plot.disturbance) {
        case FecalPat:
             _pat( sp);
             break;
        case AntMound:
             _mound( sp);
             break;
        case Burrow:
             _burrow( sp);
             break;
      }

    } /* end for each species*/

  } /* end ForEachGroup(rg) */

  *killed = _SomeKillage;
}

/***********************************************************/
void mort_EndOfYear( void)
{
/*======================================================*/
/* PURPOSE */
/* Implements killing of plants through fire or another event
 * that would result in mortality in single year (killyr) or 
 * extirpation (extirp) where plants are killed and do not return. */

	//printf("inside mort_EndOfYear() \n");
	GrpIndex rg;
	GroupType *g;
        RealF fire_possibility;
        RealF x_cheatgrass = Species[g->cheatgrass_index]->relsize * Species[g->cheatgrass_index]->mature_biomass; // calculate biomass of cheatgrass
        
        

       //printf("[Rui] x_cheatgrass: %f\n",x_cheatgrass);
        if (g->killfreq == 0)
                    { 
                         fire_possibility = 0;
                    }
        if (x_cheatgrass < g->ignition)
                    {
                       fire_possibility = 1.0 / Globals.runModelYears;
                     
                      
                    } 
        /* don't ignite wild fire based on cheatgrass until the cheatgrass could reach the ignition value */
                    if (x_cheatgrass > g->ignition){
                        fire_possibility = g->cheatgrass_coefficient + g->wild_fire_slope * x_cheatgrass;
                        }
                    if GT (g->killfreq, fire_possibility)
                    { 
                         fire_possibility = g->killfreq;
                    }
                     
          /*  if the input prescribed fire value < 1; wild fire and prescribed fire happen at the same time, 
           * while if the input prescribed fire value >1, only prescribed fire happens */
     
          printf("[Rui] fire_possibility: %f\n",fire_possibility); 
         
	ForEachGroup(rg)
	{
		if (Globals.currYear < RGroup[rg]->startyr)
		{
			/* don't start trying to kill or grow or do grazing until RGroup[rg]->startyr year */
			continue;
		}

		g = RGroup[rg];
                     
                
                    if ((Globals.currYear >= g->killfreq_startyr) && GT(fire_possibility, 0.))
                    {
			if (LT(fire_possibility, 1.0))
			{
				if (RandUni() <= fire_possibility)
				{
					g->killyr = Globals.currYear;
				}
			}
			else if (((Globals.currYear - g->killfreq_startyr) % (IntU) fire_possibility) == 0)
			{
				g->killyr = Globals.currYear;
			}

                    }

		if (Globals.currYear == RGroup[rg]->extirp)
		{
			rgroup_Extirpate(rg);
		}
		else if (Globals.currYear == RGroup[rg]->killyr)
		{
			RGroup_Kill(rg);
		}

	}

}

void grazing_EndOfYear( void){

	/*======================================================*/
    /* PURPOSE */
	/* Perform the sorts of grazing one might expect at end of year, it is based on grazing frequency
    /* HISTORY */
	/* 1st Nov 2015 -AT  -Added Species grazing EndOfYear  */
	/*======================================================*/

	GrpIndex rg;
	GroupType *g;

	ForEachGroup(rg)
	{
		IntU grazingyr =0;
		g = RGroup[rg];

		//printf("inside grazing_EndOfYear() year=%d, rgroupName=%s, grazingfreq_startyr=%d, grazingfreq=%d, proportionGrazing=%f startYear=%d \n",Globals.currYear,g->name,g->grazingfreq_startyr,g->grazingfrq,g->proportion_grazing, RGroup[rg]->startyr);

		if (Globals.currYear < RGroup[rg]->startyr)
		{
			/* don't start trying to grow or do grazing until RGroup[rg]->startyr year */
			continue;
		}

		if ((Globals.currYear >=g->grazingfreq_startyr) && GT(g->grazingfrq, 0.))
		{
			if (LT(g->grazingfrq, 1.0))
			{
				if (RandUni() <= g->grazingfrq)
				{
					grazingyr = Globals.currYear;
				}

			}
			else if (((Globals.currYear - g->grazingfreq_startyr) % (IntU) g->grazingfrq) == 0)
			{
				grazingyr = Globals.currYear;
			}

		}

		//rgroup grazing
		if (Globals.currYear == grazingyr)
		{	
			//printf( "currYear is equal to grazingYear so will iterate all the Species for doing grazing, RGroup[g]->est_count =%d \n",RGroup[rg]->est_count);
			Int i;
			ForEachEstSpp2( rg, i)
			{
				if (!Species[RGroup[rg]->est_spp[i]]->use_me)
				{
					continue;
				}
				//printf( "year=%d calling Species_Proportion_Grazing()  rgroupName=%s, est_count =%d,grazingfreq_startyr=%d, grazingfreq=%d, proportionGrazing=%f \n",Globals.currYear,g->name,RGroup[rg]->est_count,g->grazingfreq_startyr,g->grazingfrq,g->proportion_grazing);
				Species_Proportion_Grazing(RGroup[rg]->est_spp[i],RGroup[rg]->proportion_grazing );
			}
		}

	}

}

void proportion_Recovery(void)
{
	/*======================================================*/
    /* PURPOSE */
	/* Perform the sorts of proportion_Recovery one might expect at next year after the killing year
    /* HISTORY */
	/* 1st Nov 2015 -AT -Added Species Proportion Recovery  */
	/*======================================================*/

	GrpIndex rg;
	GroupType *g;

	ForEachGroup(rg)
	{

		if (Globals.currYear < RGroup[rg]->startyr)
		{
			/* don't start trying to grow  until RGroup[rg]->startyr year */
			continue;
		}

		g = RGroup[rg];

		
		//rgroup proportion recovery
		if (Globals.currYear == RGroup[rg]->killyr)
		{
			Int i;
			ForEachEstSpp2( rg, i)
			{
				Species_Proportion_Recovery(RGroup[rg]->est_spp[i], 6,
						RGroup[rg]->proportion_recovered,
						RGroup[rg]->proportion_killed);
			}
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
    IndivType *p, *kills[MAX_INDIVS_PER_SPP];

    /* ---------------------------------------------*/
    /* Generate kill list, depending on sensitivity */
    /* ---------------------------------------------*/
    if ( Plot.pat_removed) {
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

    if (k >= 0) _SomeKillage = TRUE;
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


    if (k) _SomeKillage = TRUE;
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

    if (k ) _SomeKillage = TRUE;
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
 *      and remove them properly.
/*------------------------------------------------------*/

  IndivType *p,
            *kills[MAX_INDIVS_PER_SPP];
  RealF killamt = Succulent.reduction;
  int i, k=0;

  ForEachIndiv (p, Species[sp]) {
    if ( GT(p->relsize, killamt) )
      indiv_Kill_Partial( Slow, p, killamt);
    else
      kills[k++] = p;
  }

  for (i=0; i < k; i++)
    indiv_Kill_Complete(kills[i], 7);


  if (Species[sp]->est_count) _SomeKillage = TRUE;
}


/***********************************************************/
static void _slow_growth( const SppIndex sp) {
/*======================================================*/
/* Kill plants based on a probability if the growth rate
   is less than the "slow rate" which is defined by the
   user in the group-level parameters (max_slow) and in
   the species-level parameters (max_rate).  the slow rate
   is growthrate <= max_slow * max_rate.

   Increment the counter for number of years of slow growth.
   If the number of years of slow growth is greater than
   max_slow (defined in species parms), get a random number
   and test it against the probability of death.  C&L'90
   defines this value as a constant, but it might be better
   to define it in the groups or species parameters.

   Of course, annuals aren't subject to this mortality,
   nor are new plants.

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

  Int n, k=-1;
  RealF pm = 0.368, /* probability of mortality*/
        slowrate;
  IndivType *ndv,
            *kills[MAX_INDIVS_PER_SPP];

  slowrate = RGroup[Species[sp]->res_grp]->slowrate
           * Species[sp]->max_rate;

  ForEachIndiv (ndv, Species[sp]) {
    if ( ndv->age == 1) continue;
    if (ndv->growthrate <= slowrate) {
      ndv->slow_yrs++;
      /* add to kill list if pm met*/
      if ( ndv->slow_yrs >= Species[sp]->max_slow
           && RandUni() <= pm)
         kills[++k] = ndv;
    } else
      ndv->slow_yrs = max( ndv->slow_yrs -1, 0);

  }

  for( n=0; n <= k; n++ )
    indiv_Kill_Complete(kills[n], 8);

  if (k >= 0) _SomeKillage = TRUE;
}

/***********************************************************/
static void _age_independent( const SppIndex sp) {
/*======================================================*/

/* kills possibly all individuals in a species
   by the age-independent function (eqn 14) in C&L'90
   assuming that AGEMAX was defined.
*/
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/* 5/22/01 (cwb) - don't kill annuals here or else they
      die too early in the process and stats don't work, use
      mort_EndOfYear() instead. also, skip species with
      max_age==0 (longest lived).
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
    if (RandUni() <= pn)
      kills[++k] = ndv;
  }

  for( n=0; n <= k; n++ ) {
    indiv_Kill_Complete(kills[n], 9);
  }


  if (k >= 0) _SomeKillage = TRUE;

  Mem_Free(kills);
}

/***********************************************************/
static void _no_resources( GrpIndex rg) {
/*======================================================*/
/* use EQN 7, 8, 9 prior to growing
 * Lack of resources require a reduction in "plant-space"
 * (ie, individuals or portions thereof) such that
 * plant-space balances resource-space.  In reality, this
 * would happen gradually over the season, but in this
 * yearly-time-step model it has to be done either before
 * or after the growth routine.  C&L 1990 (p241) note that
 * growth rates are also reduced--this is done in the main
 * growth loop by setting gmod = 1/PR.

 * Note that the group's PR MUST BE > 1.0 BY THIS TIME.
 * Normally this is tested for in mort_Main().

 * This routine also calls _stretched_clonal() to kill
 * additional amounts of clonal plants (if any), which also
 * happens due to insufficient resources.
 */

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

  if (nk) _SomeKillage = TRUE;

  /* Check to see if this group's resources have been stretched,
   * and commit mortality of clonal plants which get additional
   * killamt if they are not killed in the preceeding loop.
   * Pick up next largest individual from where we left off in
   * previous loop by reusing i without resetting it to 0.
   * i comes out of the loop pointing to the next living plant.
   * If for some reason all plants were killed, _stretched_clonal
   * bails before doing anything.
   */
  _stretched_clonal( rg, i, n-1, indv_list);

  Mem_Free( indv_list);

}

/**************************************************************/
static void _stretched_clonal( GrpIndex rg, Int start, Int last,
                           IndivType *nlist[]) {
/*======================================================*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

  Int i,
      y,  /* number of years of stretched resources*/
      np, /* number of clonal plants in this resource group*/
      nk; /* number of clonal plants to kill if pm met*/
  RealF pm; /* Probability of mortality (eqn 8)*/

  /* these are used if reducing proportionally (pm not met)*/
  RealF total_size,
       indiv_size,
       total_reduction,
       indiv_reduction;

  IndivType *clist[MAX_INDIVS_PER_SPP]; /* list of clonals*/

  /*-----------------------------------------*/
  /* get a list of remaining clonal plants ,*/
  /* still in order of size*/
  for( np=-1, i=start; i <= last; i++) {
    if (Species[nlist[i]->myspecies]->isclonal)
      clist[++np] = nlist[i];
  }
  if (np < 0)
    return;  /* Bail if no clonals remain alive in this group */

  /*-----------------------------------------*/
  y = RGroup[rg]->yrs_neg_pr;

  if (y >= RGroup[rg]->max_stretch) {

    pm = .04 * y * y;  /* EQN 8*/

  /*-----------------------------------------*/
    if (RandUni() <= pm ) {  /* kill on quota basis*/
      /* must be more than 10 plants for any to survive ?*/
      /* if so, then use ceil(), otherwise, use floor()*/
      nk = (Int) floor(((RealF) (np+1) * 0.9)); /* EQN 9*/

      /*  kill until we reach quota or number of plants*/
      nk = min( nk, (np+1));
      for( i = 0; i < nk; i++) {
        indiv_Kill_Complete(clist[i], 11);
      }

      if (nk >= 0) _SomeKillage = TRUE;

  /*-----------------------------------------*/
    } else {  /* reduce inverse-proportionally*/

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
      if (np >= 0) _SomeKillage = TRUE;

    } /* end if pm*/
  } /* end if y >= 1*/

}

/***********************************************************/
void _kill_annuals( void) {
/*======================================================*/
/* PURPOSE */
/* Loop through all species and kill the annuals.  This
   routine should be called at the end of the year after
   all growth happens and statistics are calculated and
   we don't need to know about the annuals any more.

   The assumption, of course, is that all of the annual
   species that are established are indeed one year old.
   See the discussion at the top of this file and in
   indiv_create() for more details.
*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 3/14/2001            */
/*
 *   7-Nov-03 (cwb) Note that annuals can't be mixed with
 *       perennials as I'd originally imagined they might.
 *       Still, the original code works, if not most
 *       efficiently.  Chances are I won't optimize it
 *       because I want to convert the whole thing to C++.

/*------------------------------------------------------*/

  GrpIndex rg;
  SppIndex sp;
  Int i;

  ForEachGroup(rg) {
    if (RGroup[rg]->max_age == 1) {
      for(i=RGroup[rg]->est_count, sp=RGroup[rg]->est_spp[i-1]; i>0; sp=RGroup[rg]->est_spp[(--i) - 1]){
               Species_Annual_Kill(sp, 4);
              //above calls new function and kill all annual individuals (TEM 10-27-2015))
          }
      }
    }

}


/***********************************************************/
void _kill_extra_growth( void) {
/*======================================================*/
/* PURPOSE */
/* Remove any extra growth accumulated during the growing
   season.  This should be done after all the statistics
   are accumulated for the year.

   Note that extra growth is not stored by the individuals
   but only by Species and RGroup.

  This should probably be split into separate functions
  so the actual killing of extra growth is one function
  that can be called from anywhere, and the looping is
  in another function with another name.
  That way, individual functions can kill specific
  extra growth without affecting others.
*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 10/26/2000            */

/*------------------------------------------------------*/

  IntU j;
  GrpIndex rg;
  SppIndex sp;


  ForEachGroup(rg) {
    if (!RGroup[rg]->use_extra_res) continue;
    ForEachEstSpp( sp, rg, j) {
      if ( ZRO(Species[sp]->extragrowth) ) continue;
      Species_Update_Newsize( sp, -Species[sp]->extragrowth);
      Species[sp]->extragrowth = 0.0;
    }
  }


}

