/********************************************************/
/********************************************************/
/*  Source file: species.c
/*  Type: module
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: Read / write and otherwise manage the
 *           site specific information.  See also the
 *           Layer module.
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
#include "ST_steppe.h"
#include "ST_globals.h"
#include "myMemory.h"
#include "rands.h"


/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
void rgroup_AddSpecies( GrpIndex rg, SppIndex sp) ;
void rgroup_DropSpecies( SppIndex sp) ;
Bool indiv_New( SppIndex sp);
void indiv_Kill_Complete( IndivType *ndv, int killType);
void indiv_proportion_Kill( IndivType *ndv, int killType,RealF proportionKilled);

void _delete (IndivType *ndv);

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* but have to be declared.                             */
void species_Update_Kills( SppIndex sp, IntS age );
void species_Update_Estabs( SppIndex sp, IntS num );
SppIndex species_New( void);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static SpeciesType *_create ( void);


/****************** Begin Function Code ********************/
/***********************************************************/


IntS Species_NumEstablish( SppIndex sp) {
/*======================================================*/
/* PURPOSE */
/* return a number: 0 or more seedlings that establish */


/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

	//special conditions if we're using the grid and seed dispersal options (as long as its not during the spinup, because we dont use seed dispersal during spinup)
	if(UseGrid && UseSeedDispersal && !DuringSpinup)
		if(Species[sp]->sd_sgerm) {
			if(Species[sp]->max_seed_estab <= 1)
				return 1;
			else
				return (IntS) RandUniRange(1, Species[sp]->max_seed_estab);		
		} else
			return 0;

	//float biomass = Species[sp]->relsize * Species[sp]->mature_biomass; //This line does nothing!
	if ( RGroup[Species[sp]->res_grp]->est_annually ||
			LE( RandUni(), Species[sp]->seedling_estab_prob) || (Species[sp]->sd_sgerm)) {
		if (Species[sp]->max_seed_estab <= 1)
			return 1;
		else
			return (IntS) RandUniRange(1, Species[sp]->max_seed_estab);
		/*    return Species[sp]->max_seed_estab; */
	} else  {
		return 0 ;
	}
}

/**************************************************************/
RealF Species_GetBiomass (SppIndex sp)  {
/*======================================================*/
/* PURPOSE */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/
  return (Species[sp]->relsize * Species[sp]->mature_biomass);
}

/**************************************************************/
void Species_Add_Indiv( SppIndex sp, Int new_indivs) {
/*======================================================*/
/* PURPOSE */
/* Add n=new_indivs individuals to the established list for
 * species=sp.  add the species to the established list for
 * the species group if needed and update the relsizes.   */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

  Int i;
  GrpIndex rg;
  RealF newsize = 0.0; /*accumulate total relsize for new indivs*/

  if (0 == new_indivs) return;

  rg = Species[sp]->res_grp;

  /* add individuals until max indivs */
  for(i=1; i<= new_indivs; i++ ) {
    if ( ! indiv_New(sp))
      LogError(logfp,LOGFATAL, "Unable to add new individual in Species_Add_Indiv()");
    Species[sp]->est_count++;
    newsize += Species[sp]->relseedlingsize;
  }

  /* add species to species group if new*/
  rgroup_AddSpecies( rg, sp);

/* accumulate sizes and resources used/available*/
  Species_Update_Newsize( sp, newsize);
}


/**************************************************************/
void species_Update_Kills( SppIndex sp, IntS age ) {
/*======================================================*/
/* PURPOSE */
/* accumulate frequencies of kills by age for survivorship */
/* 'kills' is a pointer to a dynamically sized array created
 * in params_check_species().
 *
 * Important note: age is base1 so it must be decremented
 * for the base0 arrays.
 */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 5/18/2001            */

/*------------------------------------------------------*/
	
  if(!isnull(Species[sp]->kills)) {
    age--;
    if(age >= Species[sp]->max_age)  //a quick check to keep from writing off the end of the kills arrays (which was happening in really obscure cases)...
    	return;
    Species[sp]->kills[age] ++;
    RGroup[Species[sp]->res_grp]->kills[age] ++;
  }

}

/**************************************************************/
void species_Update_Estabs( SppIndex sp, IntS num ) {
/*======================================================*/
/* PURPOSE */
/* accumulate number of indivs established by species */
/* for fecundity rates.                               */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 5/21/2001            */

/*------------------------------------------------------*/

   Species[sp]->estabs += num;
   RGroup[Species[sp]->res_grp]->estabs += num;

}

/**************************************************************/
void Species_Update_Newsize( SppIndex sp, RealF newsize ) {
/*======================================================*/
/* PURPOSE */
/* This is the point at which any changes in the individuals,
   whether by growth or mortality, is reflected in the overall
   relative size at the Species and RGroup levels.

   There is a gotcha with floating point math (at least in C;
   FORTRAN may hide this problem) in that rounding errors and
   slight inaccuracies in representing rational numbers in the
   least significant decimal places can cause values to be
   nonzero when they should be 0.0.  This is especially
   problematic when individuals grow and die partially, and
   then an individual is killed completely:  due to the
   internal representation, the value may never be exactly
   zero.  The ZRO() macro (defined in "generic.h") tests
   for a value near zero.  Likewise, tests for equality are
   subject to representational error, so there are macros for
   that as well.

*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000
 *
 *    2/26/03 - cwb - The numeric representation problem now
 *            appears to  be worse than it first seemed.  The
 *            ZRO macro uses the F_DELTA macro, defined in
 *            "generic.h".  Because F_DELTA is used in macros
 *            sprinkled everywhere, that value should be as
 *            small as practical.  However, the changes in
 *            RGroup and Species relative sizes occur only
 *            here and it turns out that the small F_DELTA is
 *            not enough to prevent rounding error from
 *            becoming a substantial problem over time.  Thus,
 *            the ZRO macro is redefined to be more inclusive.
 *    -- Also, (and more importantly) the line to correct a
 *            reduction of size more than current size was
 *            incorrect.  It was change from
 *                newsize = Species[sp]->relsize;
 *            to
 *                newsize = -Species[sp]->relsize;
 *
 *    3-Apr-03 - Big change.  Via some form of miscommunication,
 *            the group size was equal to the sum of of all the
 *            individuals' sizes, meaning that any one full-size individual
 *            used the group's full complement of space and resources.
 *            Now the assumption is the group's space and resources
 *            can support the equivalent of one full sized indiv
 *            of each species, eg, RGroup->relsize = 1.0 can contain
 *            the sum of the max biomass of the species in the group.
 *            See new function RGroup-Update_Newsize().
 *
 *
 *    7-Nov-03 (cwb) Adding new algorithms for annuals. Primarily
 *            this means that the mechanism for adding and deleting
 *            indivs is unneeded, so we have to make sure it doens't
 *            get referenced when the species is an annual.  This also
 *            means that est_count has no meaning, since we're only
 *            fiddling with the relative size.

 */

/*------------------------------------------------------*/
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
( (sizeof(x) == sizeof(float)) \
  ? ((x)>-xF_DELTA && (x)<xF_DELTA) \
  : ((x)>-xD_DELTA && (x)<xD_DELTA) )

  GrpIndex rg;

  rg = Species[sp]->res_grp;

  if ( LT(Species[sp]->relsize, 0.0) ) {
    LogError(logfp, LOGWARN, "Species_Update_Newsize: %s relsize < 0.0 (=%.1f)"
                   " year=%d, iter=%d",  Species[sp]->name, Species[sp]->relsize,
                   Globals.currYear, Globals.currIter);
  }
  if ( GT(Species[sp]->relsize, 100.) ) {
    LogError(logfp, LOGNOTE, "Species_Update_Newsize: %s relsize very large (=%.1f)"
                   " year=%d, iter=%d",  Species[sp]->name, Species[sp]->relsize,
                   Globals.currYear, Globals.currIter);
  }

  /* if this cond. true, we're off a bit from zeroing. fix it */
  if ( Species[sp]->est_count == 1 &&
       LT(newsize,-Species[sp]->relsize ) )
    newsize = -Species[sp]->relsize;

  Species[sp]->relsize += newsize;
  RGroup_Update_Newsize(rg);

  if ( Species[sp]->max_age != 1) {
    /* make sure zeros are actually zeroed */
    if (Species[sp]->est_count < 0) Species[sp]->est_count = 0;
  }
  if (ZERO(Species[sp]->relsize))  Species[sp]->relsize = 0.0;


#undef xF_DELTA
#undef xD_DELTA
#undef ZERO
}


/**************************************************************/
SppIndex species_New( void) {
/*======================================================*/
/* PURPOSE */
/* Create a new species object and give it the next
 * consecutive identifier.
 *
 * Return the index of the new object.
 *
 * Initialization is performed in parm_Species_Init
 * but the list of individuals in the species is
 * maintained by a linked list which is of course
 * empty when the species is first created.  See the
 * Indiv module for details.

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/
  SppIndex i = (SppIndex) Globals.sppCount;

  if (++Globals.sppCount > MAX_SPECIES) {
    LogError(logfp, LOGFATAL, "Too many species specified (>%d)!\n"
            "You must adjust MAX_SPECIES and recompile!\n",
            MAX_SPECIES);
  }

  Species[i] = _create();
  Species[i]->IndvHead = NULL;
  return i;
}


/**************************************************************/
static SpeciesType *_create ( void) {
/*======================================================*/
/* PURPOSE */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

  SpeciesType *p;

  p = (SpeciesType *)
      Mem_Calloc( 1, sizeof(SpeciesType), "Species_Create");

  return (p);

}

/**************************************************************/
SppIndex Species_Name2Index (const char *name) {
/*======================================================*/
/* PURPOSE */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/
  Int i, sp = -1;
  ForEachSpecies(i) {
    if (strcmp(name, Species[i]->name) == 0) {
      sp = i;
      break;
    }
  }
  return( (SppIndex)sp);
}

void Species_Proportion_Kill (const SppIndex sp, int killType, RealF proportionKilled ){
	/*======================================================*/
	/* PURPOSE */
	/* Proportion Killed all established individuals in a species.
	 *
	 * Note the special loop construct.  we have to save the
	 * pointer to next prior to killing because the object
	 * is deleted.

	/* HISTORY */
	/* Chris Bennett @ LTER-CSU 11/15/2000            */
	/*   8/3/01 - cwb - added linked list processing.
	 *   09/23/15 -AT  -Added proportionKilled
	 */

	/*------------------------------------------------------*/
#define xF_DELTA (20*F_DELTA)
#define xD_DELTA (20*D_DELTA)
#define ZERO(x) \
( (sizeof(x) == sizeof(float)) \
  ? ((x)>-xF_DELTA && (x)<xF_DELTA) \
  : ((x)>-xD_DELTA && (x)<xD_DELTA) )

	  IndivType *p = Species[sp]->IndvHead,
	            *t;

	  if (Species[sp]->max_age == 1) {
	    Species_Update_Newsize(sp, -Species[sp]->relsize);
	  } else {
	    while(p) {
	      t = p->Next;
	      indiv_proportion_Kill( p, killType,proportionKilled);
	      p = t;
	    }
	  }

	  if ( LT(Species[sp]->relsize, 0.0)  || ZERO(Species[sp]->relsize) )
		{
		//"Warning:ST_species.c Species_Proportion_Kill()  Species[sp]->relsize is either zero or negative so deleting all the individual in species and making species rel_size to zero
		    Species[sp]->relsize = 0.0 ;
			IndivType *p1 = Species[sp]->IndvHead, *t1;

			while (p1)
			{
				t1 = p1->Next;
				_delete(p);
				p1 = t1;
			}
			rgroup_DropSpecies(sp);
		}

	#undef xF_DELTA
	#undef xD_DELTA
	#undef ZERO

}

/**************************************************************/
void Species_Kill (const SppIndex sp, int killType) {
/*======================================================*/
/* PURPOSE */
/* Kill all established individuals in a species.
 *
 * Note the special loop construct.  we have to save the
 * pointer to next prior to killing because the object
 * is deleted.

/* HISTORY */
/* Chris Bennett @ LTER-CSU 11/15/2000            */
/*   8/3/01 - cwb - added linked list processing.
 */

/*------------------------------------------------------*/
  IndivType *p = Species[sp]->IndvHead,
            *t;


  if (Species[sp]->max_age == 1) {
    Species_Update_Newsize(sp, -Species[sp]->relsize);
  } else {
    while(p) {
      t = p->Next;
      indiv_Kill_Complete( p, killType);
      p = t;
    }
  }

  rgroup_DropSpecies(sp);



}



#ifdef DEBUG_MEM
#include "myMemory.h"
/*======================================================*/
void Species_SetMemoryRefs( void) {
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
  SppIndex sp;
  IndivType *p;

  ForEachSpecies(sp) {
    NoteMemoryRef(Species[sp]);
    NoteMemoryRef(Species[sp]->kills);
    p = Species[sp]->IndvHead;
    while (p) {
      NoteMemoryRef(p);
      p = p->Next;
    }
  }

}

#endif


