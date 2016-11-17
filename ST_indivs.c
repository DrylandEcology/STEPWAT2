/********************************************************/
/********************************************************/
/*  Source file: indivs.c
/*  Type: module
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: This module manages the comings and goings of
 *           individual plants, plus a function for sorting
 *           based on size.
/*  History:
/*     (6/15/2000) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdlib.h>
#include <memory.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "filefuncs.h"
#include "myMemory.h"


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

//_delete function also called from ST_species.c
void _delete (IndivType *ndv);

/***********************************************************/
/****************** Begin Function Code ********************/


/***********************************************************/
Bool indiv_New( SppIndex sp) {
/*======================================================*/
/* PURPOSE */
/* Use this routine to when a new plant is established.
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
*/
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/


  IndivType *p;
  static int id=0;

  if (Species[sp]->est_count == MAX_INDIVS_PER_SPP) {
    LogError(logfp, LOGWARN, "Limit reached: %s is about to get %d "
                   "indivs (max=%d)\n", Species[sp]->name,
                   Species[sp]->est_count +1,
                   MAX_INDIVS_PER_SPP);
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

  //sql for inserting new indiv
  if(!UseGrid)
	  insertIndiv(p);
  id++;
  return( TRUE);
}

/**************************************************************/
static IndivType *_create ( void) {
/*======================================================*/
/* PURPOSE */
/* Local routine creates object, initializes to zero and
 * returns a pointer to it.  Returned indiv has no identity
 * until called by Indiv_New().
 *
 * Guarantees the creation of a valid object or it fails.
*/
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

  IndivType *p;

  p = (IndivType *) Mem_Calloc( 1, sizeof(IndivType),
                                "indiv_create");
  return (p);

}

/**************************************************************/
Bool indiv_Kill_Partial( MortalityType code,
                          IndivType *ndv,
                          RealF killamt) {
/*======================================================*/
/* PURPOSE */
/* Clonal plants can be partially killed. If so, the type
 * of mortality must be recorded as it determines when or
 * if the plant can vegetatively propogate.  Amount of
 * damage/killage/shrinkage is subtracted from size and
 * RGroup and Species sizes are updated.
 *
 * code - stored in indiv to control vegetative propogation
 *        (see rgroup_Grow() for prop., see mort_Main() for
 *         vegetative reduction).
 * ndv  - pointer to the individual to be shrunken.
 * killamt - relative amount of the plant to be removed
 *           from the relative size variable.
*/
/* HISTORY */
/*   Chris Bennett @ LTER-CSU 6/15/2000            */
/*   cwb - 2-Dec-02 -- Another bug became apparent while
 *       adding SOILWAT code and this exists in the old
 *       model as well.  The problem is in the list
 *       processing: partial kills don't remove the
 *       indiv object whereas if killamt > relsize the
 *       object was removed from the list though the caller
 *       may not be able to know this.  My solution is
 *       to return TRUE if killamt < relsize (and relsize
 *       etc is updated) or FALSE otherwise which allows
 *       the caller to kill completely and handle the
 *       removal properly.
 *
 *   2/25/03 - Egad! Can't believe I never noticed there was
 *       no call to update the size in the higher echelons.
 *       Added call to Species_Update_Newsize().

/*------------------------------------------------------*/
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
    Species_Update_Newsize(sp, -killamt);
  }

  return( result);
}

void indiv_proportion_Kill(IndivType *ndv, int killType, RealF proportKilled)
{
	/*======================================================*/
	/* PURPOSE */
	/* Remove individual proportionally and adjust relative sizes of the
	 * RGroup and Species downward proportionally by the size of the indiv.
	 * Also keep up with survivorship data.
	 */
	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 6/15/2000
	 *  09/23/15 -AT  -Added proportionKilled
	 *  Nov 4 2015- AT - Modified code for doing proportional kill for annual as well and not deleting species 
	 *  and indi from memory only adjusting their real size  */

	/*------------------------------------------------------*/

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
			                    Globals.currIter, Globals.currYear);
	}

	if (!UseGrid)
		insertIndivKill(ndv->id, killType);

    //kill indiv Proportionally or adjust their real size irrespective of being annual or perennial, both will have this effect
	species_Update_Kills(ndv->myspecies, ndv->age);
	// saving killing year real size here that is going to use for calculating next year proportional recovery
	ndv->prv_yr_relsize = ndv->relsize;

	RealF reduction = -(ndv->relsize * proportKilled);
	printf("inside indiv_proportion_Kill() old rel_size=%f, reduction=%f \n",ndv->relsize,reduction);

	ndv->relsize = ndv->relsize + reduction;
	printf("inside indiv_proportion_Kill() new  rel_size=%f \n",ndv->relsize);
	Species_Update_Newsize(ndv->myspecies, reduction);

	if (ZERO(ndv->relsize) || LT(ndv->relsize, 0.0))
	{
		ndv->relsize =0.0;
	}

	#undef xF_DELTA
	#undef xD_DELTA
	#undef ZERO
}
void indiv_proportion_Grazing( IndivType *ndv, RealF proportionGrazing)
{
	/*======================================================*/
	/* PURPOSE */
	/* Do individual grazing proportionally and adjust relative sizes of the
	 * RGroup and Species downward proportionally by the size of the indiv.
	 * Also keep up with survivorship data.
	 */
	/* HISTORY */
	/* 1st Nov 2015- AT
	/*------------------------------------------------------*/

	#define xF_DELTA (20*F_DELTA)
	#define xD_DELTA (20*D_DELTA)
	#define ZERO(x) \
	( (sizeof(x) == sizeof(float)) \
			? ((x)>-xF_DELTA && (x)<xF_DELTA) \
					: ((x)>-xD_DELTA && (x)<xD_DELTA) )


	RealF grazing_reduce = -(ndv->relsize * proportionGrazing);
	ndv->relsize = ndv->relsize + grazing_reduce;
	Species_Update_Newsize(ndv->myspecies, grazing_reduce);

	if (ZERO(ndv->relsize) || LT(ndv->relsize, 0.0))
	{
		ndv->relsize = 0.0;
	}
	#undef xF_DELTA
	#undef xD_DELTA
	#undef ZERO
}

void indiv_proportion_Recovery( IndivType *ndv, int killType,RealF proportionRecovery,RealF proportionKilled)
{
	/*======================================================*/
	/* PURPOSE */
	/* Recover individual proportionally and adjust relative sizes of the
	 * RGroup and Species upward proportionally by the size of the indiv.
	 * Also keep up with survivorship data.
	 */
	/* HISTORY */
	/* 1st Nov 2015- AT
	/*------------------------------------------------------*/

	#define xF_DELTA (20*F_DELTA)
	#define xD_DELTA (20*D_DELTA)
	#define ZERO(x) \
	( (sizeof(x) == sizeof(float)) \
			? ((x)>-xF_DELTA && (x)<xF_DELTA) \
					: ((x)>-xD_DELTA && (x)<xD_DELTA) )

   // using  individual killing year old real size and reduction for making base for calculating proportional recovery
	RealF prev_reduction = ndv->prv_yr_relsize * proportionKilled;
  	RealF increase = prev_reduction * proportionRecovery;
	ndv->relsize = ndv->relsize + increase;
	Species_Update_Newsize(ndv->myspecies, increase);

	if (ZERO(ndv->relsize) || LT(ndv->relsize, 0.0))
	{
		_delete(ndv);
	}
	#undef xF_DELTA
	#undef xD_DELTA
	#undef ZERO
}

/**************************************************************/
void indiv_Kill_Complete( IndivType *ndv, int killType) {
/*======================================================*/
/* PURPOSE */
/* Remove individual and adjust relative sizes of the
 * RGroup and Species downward by the size of the indiv.
 * Also keep up with survivorship data.
*/
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

  if( ndv->age > Species[ndv->myspecies]->max_age ) {
    LogError(logfp, LOGWARN, "%s dies older than max_age (%d > %d). Iter=%d, Year=%d\n",
                    Species[ndv->myspecies]->name,
                    ndv->age, Species[ndv->myspecies]->max_age,
                    Globals.currIter, Globals.currYear);
  }
  if(!UseGrid)
	  insertIndivKill(ndv->id,killType);
  species_Update_Kills(ndv->myspecies, ndv->age);
  Species_Update_Newsize(ndv->myspecies, -ndv->relsize);
  _delete(ndv);

}

/**************************************************************/
void _delete (IndivType *ndv) {
/*======================================================*/
/* PURPOSE */
/* Local routine to remove the data object of an individual.
 * Called from indiv_Kill_Complete().
*/
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000
 *   a species' list of indivs is kept as a doubly linked list.
 *   although a singly linked list would work, I implemented
 *   the double just in case it might be useful in the future.
 *   as of 12/02 it hasn't been, but who knows?

/*------------------------------------------------------*/
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

  if( --s->est_count == 0)
    rgroup_DropSpecies(sp);


  if ((s->est_count > 0 && s->IndvHead == NULL)
     || (s->est_count == 0 && s->IndvHead != NULL))
     LogError(logfp, LOGFATAL,
              "PGMR: Indiv Count out of sync in _delete()");

  Mem_Free(ndv);
}

/**********************************************************/
void Indiv_SortSize( const byte sorttype,
                     const size_t n, IndivType **list) {
/*======================================================*/
/* Sort a list of pointers to individuals according to the
 * size of the individuals irrespective of species, age, etc.
 *
 * sorttype - indicates whether ascending or descending
 * n - number of individuals to be sorted
 * list[]  - an array of n pointers to individuals to be sorted.
 *           the list is returned sorted.
 *

/* HISTORY */
/* Chris Bennett @ LTER-CSU 12/15/2000            */
/*     8/2/01 - cwb - replaced shell sort with qsort(). */
/*------------------------------------------------------*/
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


/**********************************************************/
int Indiv_CompSize_A( const void *key1, const void *key2) {
/*======================================================*/
/* Comparison function for qsort, ascending order.
 * compare key1->relsize with key2->relsize.
 * if key1 < key2, return -1
 * if key1 == key2, return 0
 * if key1 > key2, return 1
 */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 8/2/01            */

/*------------------------------------------------------*/
  int r =0;
  IndivType **p1=((IndivType **)(key1)),
            **p2=((IndivType **)(key2));

  if      ( LT((*p1)->relsize, (*p2)->relsize) ) r = -1;
  else if ( GT((*p1)->relsize, (*p2)->relsize) ) r = 1;

  return r;
}

/**********************************************************/
int Indiv_CompSize_D( const void *key1, const void *key2) {
/*======================================================*/
/* Comparison function for qsort, descending order.
 * compare key1->relsize with key2->relsize.
 * if key1 < key2, return 1
 * if key1 == key2, return 0
 * if key1 > key2, return -1
 */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 8/2/01            */

/*------------------------------------------------------*/
  int r =0;
  IndivType **p1=((IndivType **)(key1)),
            **p2=((IndivType **)(key2));

  if      ( LT((*p1)->relsize, (*p2)->relsize) ) r = 1;
  else if ( GT((*p1)->relsize, (*p2)->relsize) ) r = -1;

  return r;
}

