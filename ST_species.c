/**
 * \file ST_species.c
 * \brief Contains function definitions for all species-specific functions.
 * 
 * This file modifies the \ref Species global variable and all of it's fields
 * including the linked list of [IndivTypes](\ref IndivType).
 * 
 * For the most part these functions will be used by their resgroup counterparts
 * in \ref ST_resgroup.c. There very few times were we would want to call these 
 * functions from \ref main().
 * 
 * \author
 *     Kyle Palmquist\n
 *     Chris Bennett\n
 *     Chandler Haukap\n
 *     Freddy Pierson
 * \date 22 August 2019
 * \ingroup SPECIES
 */

#include <stdlib.h>
#include <string.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sw_src/rands.h"
#include "ST_spinup.h"
#include "ST_seedDispersal.h"

extern
  pcg32_random_t species_rng;

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
void rgroup_AddSpecies(GrpIndex rg, SppIndex sp);
void rgroup_DropSpecies(SppIndex sp);
Bool indiv_New(SppIndex sp);
void indiv_Kill_Complete(IndivType *ndv, int killType);
void indiv_proportion_Kill(IndivType *ndv, int killType, RealF proportionKilled);
void indiv_proportion_Recovery(IndivType *ndv, int killType,
		RealF proportionRecovery, RealF proportionKilled);
void indiv_proportion_Grazing(IndivType *ndv, RealF proportionGrazing);

void _delete(IndivType *ndv);
void save_annual_species_relsize(void);
void copy_individual(const IndivType* src, IndivType* dest);

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* but have to be declared.                             */
void species_Update_Kills(SppIndex sp, IntS age);
void species_Update_Estabs(SppIndex sp, IntS num);
SppIndex species_New(void);
void copy_species(const SpeciesType* src, SpeciesType* dest);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static SpeciesType *_create(void);

/****************** Begin Function Code ********************/
/***********************************************************/

/**
 * \brief Calculates the number of individuals of a species to establish.
 * 
 * This function stochastically determines for each species, the number of 
 * seedlings that will establish this year based on the maximum number of 
 * seedlings that can establish in any year.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000. 
 * Overhauled by Chandler Haukap.
 * 
 * \param sp the index in \ref Species of the species to establish
 * 
 * \return A number of seedlings between 0 and \ref Species[sp]->max_seed_estab
 * 
 * \author Chandler Haukap
 * \date 19 December 2019
 * \ingroup SPECIES
 */
IntS Species_NumEstablish(SppIndex sp)
{
    // If we are using seed dispersal
    if(UseSeedDispersal && Species[sp]->use_dispersal){
        if(Species[sp]->seedsPresent && 
           RandUni(&species_rng) <= Species[sp]->seedling_estab_prob){
			// printf("** %s used dispersal **\n", Species[sp]->name);
            return (IntS) RandUniIntRange(1, Species[sp]->max_seed_estab, 
										  &species_rng);
        } else {
			// printf("%s tried dispersal but was unsuccessfull\n", Species[sp]->name);
            return 0;
        }
    }

    // If we are forcing this species to establish every year
    if(RGroup[Species[sp]->res_grp]->est_annually){
		// printf("%s forced establishment\n", Species[sp]->name);
        return (IntS) RandUniIntRange(1, Species[sp]->max_seed_estab, &species_rng);
    }

    // Otherwise, run normal establishment
    if(RandUni(&species_rng) <= Species[sp]->seedling_estab_prob){
		// printf("%s used traditional establishment\n", Species[sp]->name);
        return (IntS) RandUniIntRange(1, Species[sp]->max_seed_estab, &species_rng);
    }

    // If we didn't get caught by any of the three "if" statements above then
    // This species will not establish.
    return 0;
}

/**
 * \brief Returns the biomass in grams of a given species.
 * 
 * \param sp is the index in \ref Species of the species.
 * 
 * \return The biomass of the species. The value will always be
 *         greater than or equal to 0.
 * 
 * \ingroup SPECIES
 */
RealF Species_GetBiomass(SppIndex sp) {
	if (Species[sp]->est_count == 0) return 0.0;
	return (getSpeciesRelsize(sp) * Species[sp]->mature_biomass);
}

/**
 * \brief Adds a number of new individuals to the list of individuals of
 *        the given species.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \param sp is the index in \ref Species of the species to add individuals to.
 * \param new_indivs is the number of individuals to add.
 * 
 * \sideeffect Adds \ref Species to \ref RGroup if this species had been dropped previously.\n
 *             Modifies the linked list of individuals stored in \ref Species[ \ref sp ]
 * 
 * \ingroup SPECIES
 */
void Species_Add_Indiv(SppIndex sp, Int new_indivs)
{
	Int i;
	GrpIndex rg;
	RealF newsize = 0.0; /*accumulate total relsize for new indivs*/

	if (0 == new_indivs)
		return;

	rg = Species[sp]->res_grp;

	//printf("Inside Species_Add_Indiv() spIndex=%d, new_indivs=%d \n ",sp,  new_indivs);

	/* add individuals until max indivs */
	for (i = 1; i <= new_indivs; i++)
	{
		if (!indiv_New(sp))
		{
			LogError(logfp, LOGFATAL, "Unable to add new individual in Species_Add_Indiv()");
		}

		Species[sp]->est_count++;
		newsize += Species[sp]->relseedlingsize;
		//printf("Loop Index i=%d, Species[sp]->relseedlingsize=%.5f now newsize=%.5f \n ",i,  Species[sp]->relseedlingsize,newsize);
	}

	/* add species to species group if new*/
	rgroup_AddSpecies(rg, sp);

	//printf("Inside Species_Add_Indiv() calculated total newsize=%.5f \n ",newsize);
}

/**
 * \brief Accumulates the frequencies of kills by age for a given species.
 * 
 * Note that age should be entered in base 1.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 5/18/2001.
 * 
 * \param sp is the index in \ref Species of the requested species.
 * \param age is the age for which to update the statistics.
 * 
 * \sideeffect \ref Species[sp]->kills[age] is updated.\n
 *             \ref RGroup[Species[sp]->res_grp]->kills[age] is updated.
 * 
 * \ingroup SPECIES
 */
void species_Update_Kills(SppIndex sp, IntS age)
{
	if (!isnull(Species[sp]->kills))
	{
		age--;
		if (age >= Species[sp]->max_age) //a quick check to keep from writing off the end of the kills arrays (which was happening in really obscure cases)...
			return;
		Species[sp]->kills[age]++;
		RGroup[Species[sp]->res_grp]->kills[age]++;
	}
}

/**
 * \brief Accumulate a number of individuals established in the given species.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 5/21/2001.
 * 
 * \param sp is the index in \ref Species of the species to update.
 * \param num is the number of individuals established since the last time this function
 *            was called.
 * 
 * \sideeffect \ref Species and \ref RGroup estabs are updated.
 * 
 * \ingroup SPECIES
 */
void species_Update_Estabs(SppIndex sp, IntS num)
{
	Species[sp]->estabs += num;
	RGroup[Species[sp]->res_grp]->estabs += num;
}

/** 
 * \brief Returns the relsize of the requested species.
 * 
 * \param sp the index in the \ref Species array that you want to measure.
 * 
 * \return RealF greater than of equal to 0 representing the summed
 *         relsizes of all individuals in \ref Species[sp]
 * 
 * \sa getRGroupRelsize()
 * 
 * \ingroup SPECIES
 */
RealF getSpeciesRelsize(SppIndex sp)
{
	IndivType *p = Species[sp]->IndvHead;
    double sum = 0;

	if(p)
	{
    	while(p)
    	{
        	sum += p->relsize;
			p = p->Next;
    	}
	}

	return (RealF) (sum + Species[sp]->extragrowth);
}

/**
 * \brief Get the height of the tallest individual of this species.
 * 
 * \param sp A pointer to the \ref SpeciesType.
 * 
 * \return A float. The height of the tallest individual of the species in 
 *         centimeters.
 * 
 * \author Chandler Haukap
 * 
 * \ingroup SPECIES
 */
RealF getSpeciesHeight(SpeciesType* sp)
{
    IndivType* indiv;
    RealF maxrelsize = 0;

    // If there are no individuals established.
    if(sp->est_count < 1){
        return 0;
    }
    
    // Find the biggest individual.
    ForEachIndiv(indiv, sp){
        if(indiv->relsize > maxrelsize){
            maxrelsize = indiv->relsize;
        }
    }

	printf("within getSpeciesHeight maxrelsize = %f, Species mature_biomass = %f, maxHeight= %f, heightSlope= %f\n ", maxrelsize, sp->mature_biomass, sp->maxHeight, sp->heightSlope);

    // (maxrelsize * sp->mature_biomass) is the biomass of the individual.
    return sp->maxHeight * (1 - exp(
        -(sp->heightSlope * maxrelsize * sp->mature_biomass)));
}

/**
 * \brief Create a new species and integrate it into \ref Species.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \return index of the new species inside \ref Species.
 * 
 * \sideeffect Allocates memory for the new species and integrates it into
 * 			   \ref Species.
 * 
 * \sa Species
 * \sa _create()
 * 
 * \ingroup SPECIES
 */
SppIndex species_New(void)
{
	SppIndex i = (SppIndex) Globals->sppCount;

	if (++Globals->sppCount > MAX_SPECIES)
	{
		LogError(logfp, LOGFATAL, "Too many species specified (>%d)!\n"
				"You must adjust MAX_SPECIES and recompile!\n",
		MAX_SPECIES);
	}

	Species[i] = _create();
	Species[i]->IndvHead = NULL;
	return i;
}

/**
 * \brief Copy one species to another species.
 * 
 * Copies every field of the source \ref SpeciesType struct to the destination
 * \ref SpeciesType struct. This function also deallocates the destination's
 * \ref IndivType linked list then allocates a new \ref IndivType linked list.
 * 
 * \param src Is a pointer to the \ref SpeciesType struct to copy information 
 *            FROM.
 * \param dest Is a pointer to the \ref SpeciesType struct to copy information
 *             TO.
 * 
 * \sideeffect 
 *     The \ref IndivType linked list of the destination will be reallocated to
 *     the size of the source. And, of course, all fields of the destination 
 *     will be overwritten.
 * 
 * \ingroup SPECIES
 */
void copy_species(const SpeciesType* src, SpeciesType* dest){
	int i;
	IndivType *srcIndv, *destIndv, *next;

	// Error checking. If src == dest we would just end up loosing the linked 
    // list and some other variables.
	if(!src || src == dest){
		return;
	}

	/* ---- Reallocate and copy any arrays ---- */
	// kills: Note that this array is allocated if and only if MortFlags.summary.
	if(src->max_age > 1 && src->use_me && MortFlags.summary){
		Mem_Free(dest->kills);
		dest->kills = (IntUS*) Mem_Calloc(src->max_age, sizeof(IntUS), "copy_species: kills");
		for(i = 0; i < src->max_age; ++i){
			dest->kills[i] = src->kills[i];
		}
	}
	// seedprod
	Mem_Free(dest->seedprod);
	dest->seedprod = (IntUS*) Mem_Calloc(src->viable_yrs, sizeof(IntUS), "copy_species: seedprod");
	for(i = 0; i < src->viable_yrs; ++i){
		dest->seedprod[i] = src->seedprod[i];
	}

	/* ----------- Copy all fields ----------- */
	dest->alpha = src->alpha;
	dest->ann_mort_prob = src->ann_mort_prob;
	dest->beta = src->beta;
	dest->cohort_surv = src->cohort_surv;
	dest->disturbclass = src->disturbclass;
	dest->est_count = src->est_count;
	dest->estabs = src->estabs;
	dest->exp_decay = src->exp_decay;
	dest->extragrowth = src->extragrowth;
	dest->intrin_rate = src->intrin_rate;
	dest->isclonal = src->isclonal;
	dest->lastyear_relsize = src->lastyear_relsize;
	dest->mature_biomass = src->mature_biomass;
	dest->max_age = src->max_age;
	dest->max_rate = src->max_rate;
	dest->max_seed_estab = src->max_seed_estab;
	dest->max_slow = src->max_slow;
	dest->max_vegunits = src->max_vegunits;
	strcpy(dest->name, src->name);
	dest->prob_veggrow[0] = src->prob_veggrow[0];
	dest->prob_veggrow[1] = src->prob_veggrow[1];
	dest->prob_veggrow[2] = src->prob_veggrow[2];
	dest->prob_veggrow[3] = src->prob_veggrow[3];
	dest->pseed = src->pseed;
	dest->received_prob = src->received_prob;
	dest->relseedlingsize = src->relseedlingsize;
	dest->res_grp = src->res_grp;
    dest->maxHeight = src->maxHeight;
    dest->heightSlope = src->heightSlope;
	dest->minReproductiveSize = src->minReproductiveSize;
	dest->seedsPresent = src->seedsPresent;
    dest->maxDispersalProbability = src->maxDispersalProbability;
	dest->seedbank = src->seedbank;
	dest->seedling_biomass = src->seedling_biomass;
	dest->seedling_estab_prob = src->seedling_estab_prob;
	dest->seedling_estab_prob_old = src->seedling_estab_prob_old;
	dest->sp_num = src->sp_num;
	dest->tempclass = src->tempclass;
	dest->use_dispersal = src->use_dispersal;
	dest->use_me = src->use_me;
	dest->use_temp_response = src->use_temp_response;
	dest->var = src->var;
	dest->viable_yrs = src->viable_yrs;

	/* ------------- DEEP Copy linked list ---------------- */
	// Destroy the old linked list
	destIndv = dest->IndvHead;
	while(destIndv){
		next = destIndv->Next;
		Mem_Free(destIndv);
		destIndv = next;
	}

	srcIndv = src->IndvHead;
	// If there is a list at all.
	if(srcIndv){
		// Allocate a new individual
		destIndv = (IndivType*) Mem_Calloc(1, sizeof(IndivType), "copy_species: individual");
		// This individual is the head of the list
		dest->IndvHead = destIndv;
		// Copy the individual information across
		copy_individual(srcIndv, destIndv);
		// Since this is the head there is nothing behind it.
		destIndv->Prev = NULL;
		// While there are more individuals
		while(srcIndv->Next){
			// Move to the next individual in src
			srcIndv = srcIndv->Next;
			// Allocate the next entry in dest.
			destIndv->Next = (IndivType*) Mem_Calloc(1, sizeof(IndivType), "copy_species: individual");
			// Doubly link the list before moving on.
			destIndv->Next->Prev = destIndv;
			// Move to the new entry
			destIndv = destIndv->Next;
			// copy the individual information across
			copy_individual(srcIndv, destIndv);
		}
		// srcIndv->Next is false: we are at the end of the list.
		destIndv->Next = NULL;
	}
}

/**
 * \brief Allocates memory for a new \ref SpeciesType.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \return pointer to the newly allocated struct.
 * 
 * Called from \ref species_New().
 * 
 * \sideeffect A \ref SpeciesType struct is allocated.
 * 
 * \ingroup SPECIES_PRIVATE
 */
static SpeciesType *_create(void)
{
	SpeciesType *p;

	p = (SpeciesType *) Mem_Calloc(1, sizeof(SpeciesType), "Species_Create");
        p->name = Mem_Calloc(SuperGlobals.max_speciesnamelen + 1, sizeof(char), "Species_Create");

	return (p);
}

/**
 * \brief Converts a species name into it's location in \ref Species.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * \param name is the name of the species.
 * 
 * \return The index in \ref Species of the species with this name. 
 * \return -1 if the given name is not the name of any species.
 * 
 * \ingroup SPECIES
 */
SppIndex Species_Name2Index(const char *name)
{
	Int i, sp = -1;
	ForEachSpecies(i)
	{
		if (strcmp(name, Species[i]->name) == 0)
		{
			sp = i;
			break;
		}
	}
	return ((SppIndex) sp);
}

/**
 * \brief Kills all individuals in the given species.
 * 
 * This is meant to be called for annual species only.
 * 
 * \param sp the index in \ref Species of the requested annual species.
 * \param killType is a deprecated parameter that does nothing.
 * 
 * \sideeffect The entire linked list of individuals is deleted.\n
 *             The species is dropped from the RGroup.
 * 
 * \sa rgroup_DropSpecies() which is called from this function.
 * 
 * \ingroup SPECIES
 */
void Species_Annual_Kill(const SppIndex sp, int killType)
{
	IndivType *p1 = Species[sp]->IndvHead, *t1;
	while (p1)
	{
		t1 = p1->Next;
		_delete(p1);
		p1 = t1;
	}
	rgroup_DropSpecies(sp);
}

/**
 * \brief Kill some proportion of the given species.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 6/15/2000.
 * 
 * This proportion is removed from every individual in the species. 
 * 
 * \param sp is the index in \ref Species of the individual.
 * \param killType is deprecated and does nothing.
 * \param proportionKilled is a value between 0 and 1 denoting the proportion of
 *                         relative size to remove.
 * 
 * \sideeffect All individuals in the \ref Species[sp]->IndivHead linked list 
 *             will be reduced in size.
 * 
 * \ingroup SPECIES
 */
void Species_Proportion_Kill(const SppIndex sp, int killType,
		RealF proportionKilled)
{
	IndivType *p = Species[sp]->IndvHead, *t;
	//kill  all the species individuals  proportionally or adjust their real size irrespective of being annual or perennial, both will have this effect
	while (p)
	{
		t = p->Next;
		indiv_proportion_Kill(p, killType, proportionKilled);
		p = t;
	}
}

/**
 * \brief Performs grazing on a given species at a given proportion.
 * 
 * This function will perform grazing on every individual in the given species.
 * It will also graze the extra (superfluous) growth stored at the species level.
 * 
 * Initial programming by AT on 1/11/15.
 * Extra growth grazing added by Chandler Haukap 14/8/18.
 * 
 * \param sp is the index in \ref Species of the species.
 * \param proportionGrazing is the proportion of relative size to remove. Value
 *                          between 0 and 1.
 * 
 * \sideeffect All individuals in the given species' relsizes will be reduced
 *             proportionally.\n
 *             \ref Species[sp]->relsize will be reduced by the same proportion.
 * 
 * \sa indiv_proportion_Grazing() which is the function this function calls to 
 *                                perform grazing on each individual.
 * 
 * \ingroup SPECIES
 */
void Species_Proportion_Grazing(const SppIndex sp, RealF proportionGrazing)
{
	//CH- extra growth is only stored at the species level. This will graze extra
	//    growth for the whole species.
	//    loss represents the proportion of extragrowth that is eaten by livestock
	RealF loss = Species[sp]->extragrowth * proportionGrazing;

	//CH- Remove the loss from Species extra growth.
	Species[sp]->extragrowth -= loss;	// remove the loss from extragrowth

	//Implement grazing on normal growth for all individuals in each species.
	IndivType *t, *p = Species[sp]->IndvHead;
	while (p) //while p points to an individual
	{
		t = p->Next; //must store Next since p might be deleted at any time.
		indiv_proportion_Grazing(p, proportionGrazing);
		p = t; //move to the next plant.
	}
}

/**
 * \brief Proportion Recovery (representing re-sprouting after fire) for all 
 *        established individuals in a species.
 * 
 * \param sp is the index in \ref Species of the species.
 * \param killType is the \ref MortalityType code that killed this individual.
 * \param proportionRecovery The proportion of previous (before the fire) biomass 
 *                           to recover. Value between 0 and 1.
 * \param proportionKilled is the proportion of biomass removed by the disturbance
 *                         event. Value between 0 and 1.
 * 
 * \sideeffect All individual's relative sizes will be increased.
 * 
 * \sa indiv_proportion_Recovery() which is called by this function to perform
 *                                 the recovery on an individual level.
 * 
 * \ingroup SPECIES
 */
void Species_Proportion_Recovery(const SppIndex sp, int killType,
        RealF proportionRecovery, RealF proportionKilled) {
    IndivType *p = Species[sp]->IndvHead, *t;
    //Recover biomass for each perennial species that is established
    while (p) {
        t = p->Next;
        indiv_proportion_Recovery(p, killType, proportionRecovery,
                proportionKilled);
        p = t;
    }
    //printf("'within proportion_recovery after first killing': Species = %s, relsize = %f, est_count = %d\n",Species[sp]->name, Species[sp]->relsize, Species[sp]->est_count);
}

/**
 * \brief Kills all established individuals in a species.
 * 
 * The species will be completely killed and dropped from it's RGroup.
 * 
 * Initial programming by Chris Bennett @ LTER-CSU 11/15/2000.
 * 
 * \param sp is the index in \ref Species of the species.
 * \param killType is the \ref MortalityType code that killed the species.
 * 
 * \sideeffect 
 *    \* Every individual in the linked list will be deleted.
 *    \* The species will be dropped from it's RGroup.
 * 
 * \sa rgroup_DropSpecies() for what happens when a species is dropped.
 *     indiv_Kill_Complete() for how each individual is killed.
 * 
 * \ingroup SPECIES
 */
void Species_Kill(const SppIndex sp, int killType)
{
    IndivType *p = Species[sp]->IndvHead, *t;

    while (p)
    {
        t = p->Next;
        indiv_Kill_Complete(p, killType);
        p = t;
    }
    
    rgroup_DropSpecies(sp);
}

/**
 * \brief Records the relative size of all annual species.
 * 
 * Loops through all species and if the species is an annual it's relative size is
 * saved to Species[sp]->lastyear_relsize. This function is useful for determining
 * what proportion of relative size was killed due to disturbances as well as 
 * determining annual establishment.
 * 
 * \sideeffect the lastyear_relsize field of \ref Species is populated for all
 *             annual species.
 * 
 * \ingroup SPECIES
 */
void save_annual_species_relsize() {
    int sp = 0;

    ForEachSpecies(sp) {
        if (Species[sp]->max_age == 1) {
            Species[sp]->lastyear_relsize = getSpeciesRelsize(sp);
        }
    }
}

#ifdef DEBUG_MEM
#include "sw_src/myMemory.h"
/*======================================================*/
void Species_SetMemoryRefs( void)
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
	SppIndex sp;
	IndivType *p;

	ForEachSpecies(sp)
	{
		NoteMemoryRef(Species[sp]);
		NoteMemoryRef(Species[sp]->kills);
		p = Species[sp]->IndvHead;
		while (p)
		{
			NoteMemoryRef(p);
			p = p->Next;
		}
	}

}

#endif
