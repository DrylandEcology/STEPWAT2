/**
 * \file ST_indivs.c
 * \brief Manages plant individuals.
 * 
 * \author
 *     Kyle Palmquist\n
 *     Chandler Haukap\n
 *     Freddy Pierson\n 
 *     Chris Bennett\n
 *     Ashish Tiwari
 * 
 * \date 23 August 2019
 * 
 * \ingroup INDIVIDUAL
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdlib.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"


/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
void rgroup_DropSpecies( SppIndex sp) ;
void species_Update_Kills( SppIndex sp, IntS age );

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
Bool indiv_New( SppIndex sp);
void copy_individual(const IndivType* src, IndivType* dest);
Bool indiv_Kill_Partial( MortalityType code,
                          IndivType *ndv,
                          RealF killamt);
void indiv_Kill_Complete( IndivType *ndv, int killType);
void indiv_proportion_Kill( IndivType *ndv, int killType,RealF proportionKilled);
void indiv_proportion_Recovery( IndivType *ndv, int killType,RealF proportionRecovery,RealF proportionKilled);
void indiv_proportion_Grazing( IndivType *ndv, RealF proportionGrazing);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static IndivType *_create ( void);
void _delete (IndivType *ndv);

/***********************************************************/
/****************** Begin Function Code ********************/

/**
 * \brief Add an individual to the Species->IndivHead linked list.
 * 
 * \param sp The SppIndex of the species that this individual belongs to.
 * 
 * Use this routine to when a new plant is established.
 * Calls _create() to allocate the object, then the
 * individuals parameters are initialized. The Species
 * "object" contains a list of its individuals which is
 * updated here.  See also Indiv_Kill_Complete() for removing
 * individuals.
 *
 * Individuals are created at age == 1; the model assumes
 * one year has passed for a plant to become established.
 * The age is incremented at the end of both growth and
 * mortality for the current year so plants killed in their
 * first established growing season are 1 year old.  The
 * massive year 0 mortality one expects in nature is
 * accounted for by the probability of establishment.
 * But be aware of the base0 problem with age-related arrays.
 *
 * See rgroup_Establish() for establishment process.
 *
 * The Species obj keeps an array of pointers to allocated
 * indiv objects for normal access to the individuals. The
 * Head pointer is set to null in Species_New().
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \return TRUE
 * 
 * \sideeffect The new individual becomes the first individual in
 *             the Species linked list.
 * 
 * \sa Species
 * \sa SppIndex
 * \sa rgroup_Establish()
 * \sa Indiv_Kill_Complete()
 * 
 * \ingroup INDIVIDUAL
 */
Bool indiv_New( SppIndex sp) {
  IndivType *p;
  static int id=0;

  if (Species[sp]->est_count == SuperGlobals.max_indivs_per_spp) {
    LogError(logfp, LOGWARN, "Limit reached: %s is about to get %d "
                   "indivs (max=%d)\n", Species[sp]->name,
                   Species[sp]->est_count +1,
                   SuperGlobals.max_indivs_per_spp);
  }

  p = _create();
  p->id = id;
  p->myspecies = sp;
  p->killed = FALSE;
  p->age = 1;
  p->slow_yrs = 0;
  p->relsize = Species[sp]->relseedlingsize;

  /* link the new indiv to the head of the species list*/
  /* newer objects are closer to head */
  p->Next = Species[sp]->IndvHead;
  if (p->Next != NULL)
    p->Next->Prev = p;
  p->Prev = NULL;
  Species[sp]->IndvHead = p;

  // This functionality is unused, but it might be useful in the future.
  //sql for inserting new indiv
  //if(!UseGrid)
	//  insertIndiv(p);

  id++;
  return TRUE;
}

/** \brief Copy one individual's information to another individual. 
 * 
 *  \param src is the source [IndivType](\ref IndivType) to copy from.
 *  \param dest is the destination [IndivType](\ref IndivType) to copy to.
 * 
 *  Note: this does not modify either individual's linked list functionality.
 * 
 *  Both individuals MUST be allocated prior to calling this function. 
 */
void copy_individual(const IndivType* src, IndivType* dest){
  dest->id = src->id;
  dest->normal_growth = src->normal_growth;
  dest->pr = src->pr;
  dest->prob_veggrow = src->prob_veggrow;
  dest->prv_yr_relsize = src->prv_yr_relsize;
  dest->relsize = src->relsize;
  dest->res_avail = src->res_avail;
  dest->res_extra = src->res_extra;
  dest->res_required = src->res_required;
  dest->slow_yrs = src->slow_yrs;
  dest->yrs_neg_pr = src->yrs_neg_pr;
  dest->age = src->age;
  dest->growthrate = src->growthrate;
  dest->grp_res_prop = src->grp_res_prop;
  dest->killed = src->killed;
  dest->killedby = src->killedby;
  dest->mm_extra_res = src->mm_extra_res;
  dest->myspecies = src->myspecies;
}

/**
 * \brief Creates a new IndivType object.
 * 
 * Local routine creates object, initializes to zero and
 * returns a pointer to it. Returned indiv has no identity
 * until called by Indiv_New().
 * 
 * \return pointer to the individual.
 * 
 * \sa indiv_New() which calls this function.
 * 
 * \author Chris Bennett in 2000
 * 
 * \ingroup INDIVIDUAL_PRIVATE
 */
static IndivType *_create ( void) {
  IndivType *p;

  p = (IndivType *) Mem_Calloc( 1, sizeof(IndivType),
                                "indiv_create");
  return (p);

}

/**
 * \brief Partially kills a single individual.
 * 
 * Clonal plants can be partially killed. If so, the type
 * of mortality must be recorded as it determines when or
 * if the plant can vegetatively propogate.  Amount of
 * damage/killage/shrinkage is subtracted from size.
 * 
 * \param code The type of mortality causing size to be reduced.
 * \param ndv A pointer to the individual.
 * \param killamt The amount of relative size to remove.
 * 
 * \return TRUE if relative size was reduced. FALSE if killamt
 *         is greater than or equal to the relative size of the
 *         plant (which requires indiv_Kill_Complete()), or if
 *         the individual isn't clonal.
 * 
 * \sideeffect The size of ndv is reduced and it's growth rate is modified.
 * 
 * \sa indiv_Kill_Complete()
 * \sa IndivType
 * 
 * \ingroup INDIVIDUAL
 */
Bool indiv_Kill_Partial( MortalityType code,
                          IndivType *ndv,
                          RealF killamt) {
  SppIndex sp;
  Bool result = FALSE;

  sp = ndv->myspecies;
  if ( GT(ndv->relsize, killamt) && Species[sp]->isclonal) {
    result            = TRUE;
    ndv->killed       = TRUE;
    ndv->relsize     -= killamt;
    ndv->killedby     = code;
    ndv->growthrate   = 0.0;
    ndv->prob_veggrow = Species[sp]->prob_veggrow[code];
  }

  return( result);
}

/**
 * \brief Kill a portion of an individual plant.
 * 
 * Remove individual proportionally and adjust relative size.
 * Also keep up with survivorship data.
 * 
 * \param ndv A pointer to the individual.
 * \param killType The MortalityType code. This parameter is currently unused.
 * \param proportKilled value between 0 and 1. The percent total biomass to remove.
 * 
 * \sideeffect ndv->relsize is adjusted. 
 * 
 * \ingroup INDIVIDUAL
 */
void indiv_proportion_Kill(IndivType *ndv, int killType, RealF proportKilled)
{
    #define xF_DELTA (20*F_DELTA)
	#define xD_DELTA (20*D_DELTA)
	#define ZERO(x) \
		( (sizeof(x) == sizeof(float)) \
				? ((x)>-xF_DELTA && (x)<xF_DELTA) \
						: ((x)>-xD_DELTA && (x)<xD_DELTA) )

	if( ndv->age > Species[ndv->myspecies]->max_age )
	{
	     LogError(logfp, LOGWARN, "%s dies older than max_age (%d > %d). Iter=%d, Year=%d\n",
			                    Species[ndv->myspecies]->name,
			                    ndv->age, Species[ndv->myspecies]->max_age,
			                    Globals->currIter, Globals->currYear);
	}

	//if (!UseGrid)
	//	insertIndivKill(ndv->id, killType);

    //kill indiv Proportionally or adjust their real size irrespective of being annual or perennial, both will have this effect
	// saving killing year real size here that is going to use for calculating next year proportional recovery
	ndv->prv_yr_relsize = ndv->relsize;

	RealF reduction = -(ndv->relsize * proportKilled);
//	printf("inside indiv_proportion_Kill() old indiv rel_size=%f, reduction=%f \n",ndv->relsize,reduction);

	ndv->relsize = ndv->relsize + reduction;
//	printf("inside indiv_proportion_Kill() new indiv rel_size=%f \n",ndv->relsize);

	if (ZERO(ndv->relsize) || LT(ndv->relsize, 0.0))
	{
		ndv->relsize =0.0;
		// increase mortality count only if relsize has become zero due to fire
		species_Update_Kills(ndv->myspecies, ndv->age);
	}

	#undef xF_DELTA
	#undef xD_DELTA
	#undef ZERO
}

/**
 * \brief Reduces biomass of an individual proportionally
 * 
 * Implement grazing for each individual. Also keep up with survivorship data.
 * 
 * \param ndv A pointer to the individual.
 * \param proportionGrazing Value between 0 and 1. The proportion of biomass to remove.
 * 
 * \sideeffect ndv->relsize is adjusted.
 * 
 * \sa Species_Proportion_Grazing()
 * 
 * \ingroup INDIVIDUAL
 */
void indiv_proportion_Grazing( IndivType *ndv, RealF proportionGrazing)
{
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
	( (sizeof(x) == sizeof(float)) \
			? ((x)>-xF_DELTA && (x)<xF_DELTA) \
					: ((x)>-xD_DELTA && (x)<xD_DELTA) )

    RealF grazing_reduce = -(ndv->normal_growth * proportionGrazing);
    //	printf("inside indiv_proportion_Grazing() old indiv rel_size=%f, grazing_reduce=%f \n",ndv->relsize,grazing_reduce);
    ndv->relsize = ndv->relsize + grazing_reduce;
    //	printf("inside indiv_proportion_Grazing() new indiv rel_size=%f \n",ndv->relsize);

	if (ZERO(ndv->relsize) || LT(ndv->relsize, 0.0))
	{
        ndv->relsize = 0.0;
    }
#undef xF_DELTA
#undef xD_DELTA
#undef ZERO
}

/**
 * \brief Recover some biomass proportionally after a disturbance event.
 * 
 * Recover individuals proportionally after fire.
 * Also keep up with survivorship data. 
 * 
 * \param ndv Pointer to the individual.
 * \param killType The MortalityType code of what killed the individual.
 * \param proportionRecovery Value between 0 and 1. The proportion of individual relsize to recover.
 * \param proportionKilled Value between 0 and 1. The proportion of the individual killed
 *                         by the disturbance event.
 * 
 * \sideeffect ndv->relsize is modified.
 * 
 * \sa Species_Proportion_Recovery
 * 
 * \ingroup INDIVIDUAL
 */
void indiv_proportion_Recovery(IndivType *ndv, int killType, RealF proportionRecovery, RealF proportionKilled) 
{
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
	( (sizeof(x) == sizeof(float)) \
			? ((x)>-xF_DELTA && (x)<xF_DELTA) \
					: ((x)>-xD_DELTA && (x)<xD_DELTA) )

    /* Utilize individual relsize saved before killing and proportionKilled 
     * to determine proportional recovery */
    RealF prev_reduction = ndv->prv_yr_relsize * proportionKilled;
    RealF increase = prev_reduction * proportionRecovery;

    //printf("previous year ndv->relsize before = %f\n, Species = %s \n", ndv->prv_yr_relsize, Species[ndv->myspecies]->name);
    //printf("ndv->relsize before = %f\n, Species = %s \n", ndv->relsize, Species[ndv->myspecies]->name);
    //printf("increase = %f\n, Species = %s \n", increase, Species[ndv->myspecies]->name);

    ndv->relsize = ndv->relsize + increase;
    //printf("ndv->relsize after = %f\n,Species = %s \n", ndv->relsize, Species[ndv->myspecies]->name);

    /* This should never happen because proportion recovered should always be 
     * positive or zero */
    if (LT(ndv->relsize, 0.0)) {
        // this should never happen because `increase` should always be positive
        LogError(logfp, LOGWARN, "'indiv_proportion_Recovery': an individual of "\
      "%s reached relsize < 0 (increase = %.3f): for " \
      "killType = %d, proportionKilled = %.2f, proportionRecovery = %.2f",
                Species[ndv->myspecies]->name, increase,
                killType, proportionKilled, proportionRecovery);
    }

    /* If the relsize is still zero or is less than zero, remove the individual */
    if (ZERO(ndv->relsize) || LT(ndv->relsize, 0.0)) {
        _delete(ndv);
    }

#undef xF_DELTA
#undef xD_DELTA
#undef ZERO
}

/**
 * \brief Completely kill an individual.
 * 
 * Kills the individual, deallocates it, and removes it from the linked
 * list of individuals stored in Species.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \param ndv Pointer to the individual to kill.
 * \param killType MortalityType code. This parameter is currently unused.
 * 
 * \sideeffect The Species[ndv->myspecies]->IndivHead linked list is updated.\n
 *             The Species[ndv->myspecies]->est_count is updated.\n
 *             The individual is deallocated.
 * 
 * \ingroup INDIVIDUAL
 */
void indiv_Kill_Complete( IndivType *ndv, int killType) 
{
  if( ndv->age > Species[ndv->myspecies]->max_age ) {
    LogError(logfp, LOGWARN, "%s dies older than max_age (%d > %d). Iter=%d, Year=%d\n",
                    Species[ndv->myspecies]->name,
                    ndv->age, Species[ndv->myspecies]->max_age,
                    Globals->currIter, Globals->currYear);
  }
 // if(!UseGrid)
//	  insertIndivKill(ndv->id,killType);
  species_Update_Kills(ndv->myspecies, ndv->age);

  _delete(ndv); // `_delete` updates the Species[ndv->myspecies]->est_count, i.e., removes one individual

}

/**
 * \brief Deletes an individual.
 * 
 * \param ndv A pointer to the individual.
 * 
 * \sideeffect Species[ndv->myspecies]->est_count is updated.\n
 *             The Species[ndv->myspecies]->IndivHead linked list us updated.
 * 
 * \sa indiv_Kill_Complete() where this function is called.
 * 
 * \ingroup INDIVIDUAL_PRIVATE
 */
void _delete (IndivType *ndv) 
{
  SppIndex sp;
  SpeciesType *s;

  sp = ndv->myspecies;
  s = Species[sp];

  #ifdef DEBUG_MEM
    assert(fValidPointer(ndv, sizeof(IndivType)));
  #endif

  /* Detach indiv's data object from list */
  if (ndv == s->IndvHead) {
    if (ndv->Next == NULL)
      s->IndvHead = NULL;
    else {
      s->IndvHead = ndv->Next;
      s->IndvHead->Prev = NULL;
    }
  } else {
    ndv->Prev->Next = ndv->Next;
    if (ndv->Next != NULL)
      ndv->Next->Prev = ndv->Prev;
  }

  // update Species[ndv->myspecies]->est_count, i.e.,
  // remove one individual from tally
  if( --s->est_count == 0) {
    // if there are no individual left of this species, then remove species
    // from resource group and update `est_count` of the resource group
    rgroup_DropSpecies(sp);
  }

  if ((s->est_count > 0 && s->IndvHead == NULL)
     || (s->est_count == 0 && s->IndvHead != NULL))
     LogError(logfp, LOGFATAL,
              "PGMR: Indiv Count out of sync in _delete()");

  Mem_Free(ndv);
}

/**
 * \brief Sort a list of pointers to individuals according to the
 *             size of the individuals.
 * 
 * This function is irrespective of species, age, etc.
 * 
 * \param sorttype indicates whether ascending or descending
 * \param n number of individuals to be sorted
 * \param list an array of n pointers to individuals to be sorted.
 * 
 * \sideeffect The list is returned sorted. 
 * 
 * \ingroup INDIVIDUAL
 */
void Indiv_SortSize( const byte sorttype,
                     const size_t n, IndivType **list) {
  int (*cmpfunc)(const void*, const void*);

  if ( n < 1) return;  /* shouldn't happen */
  switch (sorttype) {
    case SORT_A: cmpfunc = Indiv_CompSize_A; break;
    case SORT_D: cmpfunc = Indiv_CompSize_D; break;
    default:
      LogError(logfp, LOGFATAL,
             "Invalid sort mode in Indiv_SortSize");
  }

  qsort( list,
         n,
         sizeof(IndivType **),
         cmpfunc
       );
/*
  for(i=0;i<=n;i++)printf("%d:spp=%d, size=%5.4f\n",
                       i,list[i]->myspecies,list[i]->relsize);
*/
}


/**
 * \brief Comparison function for qsort, ascending order
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 8/2/01.
 * 
 * \param key1 Pointer to first IndivType to compare.
 * \param key2 Pointer to first IndivType to compare.
 * 
 * \return -1 if key1 < key2
 * \return 0 if key1 == key2
 * \return 1 if key1 > key2
 * 
 * \ingroup INDIVIDUAL
 */
int Indiv_CompSize_A( const void *key1, const void *key2) {
  int r =0;
  IndivType **p1=((IndivType **)(key1)),
            **p2=((IndivType **)(key2));

  if      ( LT((*p1)->relsize, (*p2)->relsize) ) r = -1;
  else if ( GT((*p1)->relsize, (*p2)->relsize) ) r = 1;

  return r;
}

/**
 * \brief Comparison function for qsort, descending order
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 8/2/01.
 * 
 * \param key1 Pointer to first IndivType to compare.
 * \param key2 Pointer to first IndivType to compare.
 * 
 * \return 1 if key1 < key2
 * \return 0 if key1 == key2
 * \return -1 if key1 > key2
 * 
 * \ingroup INDIVIDUAL
 */
int Indiv_CompSize_D( const void *key1, const void *key2) {
  int r =0;
  IndivType **p1=((IndivType **)(key1)),
            **p2=((IndivType **)(key2));

  if      ( LT((*p1)->relsize, (*p2)->relsize) ) r = 1;
  else if ( GT((*p1)->relsize, (*p2)->relsize) ) r = -1;

  return r;
}
