/**
 * \file ST_defines.h
 * \brief Generic model definitions such as constants,
 *        enums, and looping contructs.
 * 
 * \author
 *     Kyle Palmquist\n
 *     Chandler Haukap\n
 *     Freddy Pierson\n 
 *     Chris Bennett
 * 
 * \date 23 August 2019
 * 
 * \ingroup STEPPE
 */

#ifndef STEPPE_DEF_H
#define STEPPE_DEF_H

#include "generic.h"

/* see #include "ST_structs.h" below */

/**
 * \brief STEPWAT2's version of true. Can be used with any Bool or int.
 * \ingroup STEPPE
 */
#define TRUE 1
/**
 * \brief STEPWAT2's version of false. Can be used with any Bool or int.
 * \ingroup STEPPE
 */
#define FALSE 0

/***************************************************
 * Basic definitions
 ***************************************************/
/**
 * \brief Macro for the maximum number of species allowed in the simulation.
 * 
 * Replaced by (SuperGlobals.max_spp_per_grp * SuperGlobals.max_rgroups) which are both
 * read from inputs.
 * 
 * \ingroup SPECIES
 */
#define MAX_SPECIES (SuperGlobals.max_spp_per_grp * SuperGlobals.max_rgroups)

/**
 * \brief MAX_OUTFIELDS The maximum number of fields our output functions are 
 *                      capable of outputting.
 * \ingroup STEPPE 
 */
#define MAX_OUTFIELDS ((MAX_SPECIES * 2) + (SuperGlobals.max_rgroups * 6) + 8)

/**
 * \brief Maximum size of an output field.
 * 
 * For numbers this isn't very important, but for header strings, like "sagebrush_bmass",
 * it is important to leave enough room for long names.
 * 
 * This is defined as Globals.max_groupnamelen + 6 currently, but can be increased
 * as needed.
 * 
 * \ingroup STEPPE
 */
#define MAX_FIELDLEN (SuperGlobals.max_groupnamelen + 6)  /* +6 for xtra chars like _RSize, etc */
/**
 * \brief defines the maximum number of cells in gridded mode.
 * \ingroup GRIDDED
 */
#define MAX_CELLS 10000

/**
 * \brief Means "sort this array in ascending order".
 * 
 * \sa Indiv_SortSize() is one function that could take this
 *                      as a parameter.
 * 
 * \ingroup STEPPE
 */
#define SORT_A 'a'
/**
 * \brief Means "sort this array in descending order".
 * 
 * \sa Indiv_SortSize() is one function that could take this
 *                      as a parameter.
 * 
 * \ingroup STEPPE
 */
#define SORT_D 'd'
/**
 * \brief Means "don't sort this array".
 * 
 * \ingroup STEPPE
 */
#define SORT_0 '\0'

/*********** define types for indices to data structures *************/
/**
 * \brief Used to iterate across resource groups stored in the \ref RGroup global variable.
 * 
 * GrpIndex exists to clarify what a variable is meant to access.
 * Any variable of type GrpIndex is meant to access \ref RGroup.
 * 
 * \sa RGroup
 * 
 * \ingroup RGROUP
 */
typedef IntS GrpIndex;
/**
 * \brief Used to iterate across species stored in the \ref Species global variable.
 * 
 * SppIndex exists to clarify what a variable is meant to access.
 * Any variable of type SppIndex is meant to access \ref Species.
 * 
 * \sa Species
 * 
 * \ingroup SPECIES
 */
typedef IntS SppIndex;


/********************** Define enumerators ***************************/
/* FYI, enums start at 0 by default.
 * I put a few "firsts" and "lasts" in here to allow
 * variable length enums where the code doesn't have
 * to know in advance, but I haven't used them much
 * yet.
 */

/**
 * \brief Enumerator for temperature class.
 * 
 * Each species has a temperature class.
 * 
 * \sa species_st which instanciates this enumerator.
 * 
 * \ingroup STEPPE
 */
typedef enum {NoSeason=-1, CoolSeason, WarmSeason}
  TempClass;

/**
 * \brief All types of disturbances.
 * 
 * Used to determine what type of disturbance is on the plot.
 * 
 * \sa plot_st which instanciates this enumerator.
 * 
 * \ingroup MORTALITY
 */
typedef enum {NoDisturb, FecalPat, AntMound, Burrow, LastDisturb}
  DisturbEvent;

/**
 * \brief All sensitivity levels to disturbances.
 * 
 * Used to determine how severely a disturbance will affect a species.
 * 
 * \sa species_st which instanciates this enumerator.
 * 
 * \ingroup MORTALITY
 */
typedef enum {VerySensitive, Sensitive, Insensitive, VeryInsensitive}
  DisturbClass;

/**
 * \brief All possible precipitation classes.
 * 
 * Used to determine whether this year's precipitation is above average, 
 * average, or below average.
 * 
 * \sa environs_st which instantiates this enumerator.
 * 
 * \ingroup STEPPE
 */
typedef enum {Ppt_Wet, Ppt_Norm, Ppt_Dry}
  PPTClass;

/**
 * \brief It is unclear why this was created.
 * 
 * It is used in ST_environs.c to access Succulent arrays
 * and Globals arrays, but only Intcpt and Slope are used.
 * Also, This enumerator is never explicitly instanciated.
 * 
 * \sa _make_disturbance() which is an example of a function
 *                         that uses Slope.
 * 
 * \ingroup STEPPE
 */
typedef enum {Intcpt, Slope, P0=0, P1, P2, P3, P4}
  Params;

/**
 * \brief How deep a group's roots penetrate.
 * 
 * Defined at the rgroup level.
 * 
 * \sa resourcegroup_st which instantiates this enumerator.
 * 
 * \ingroup STEPPE
 */
typedef enum {DepthNonComp, DepthShallow, DepthMedium, DepthDeep, DepthLast}
  DepthClass;

/**
 * \brief Enumerates all of the different input files.
 * 
 * ST_params.c stores the file names according to this convention.
 * 
 * \sa ST_params.c for its usage.
 * 
 * \ingroup STEPPE
 */
typedef enum {F_First, F_Log, F_Model, F_Env, F_Plot, F_RGroup, F_Species,
              F_BMassFlag, F_BMassPre, F_BMassAvg,
              F_MortFlag,  F_MortPre,  F_MortAvg,
              F_SXW, F_MaxRGroupSpecies, F_EXE}
  ST_FileIndex;

/**************************************************************/
/* This section defines some shorthands for loops.  There is a
 * serious tendency in this program toward off-by-one errors
 * because some arrays, eg, the main plant data structures such
 * as RGroup[], RGroup[]->est_spp, etc etc, work better with a
 * base of 1 (ie, no plant #0 exists) and some, such as the
 * temporary arrays for sorting etc, are much more naturally
 * treated as normal C arrays.  This was a design issue made
 * at the start of the process--to think of accessing groups,
 * species, and plants as #1, #2, etc--it has helped much more
 * than it hindered, but it remains to be seen whether that
 * will hold out for someone else.
 *
 * To mitigate the negative, I created these macros to maintain
 * consistency and to make the code somewhat more readable.
 *
 * Well, eventually the limitations became horribly apparent.
 * I spent a couple of months debugging diabolical errors that
 * taught me a lot about memory management.  Number 1: don't
 * mix base0 and base1 arrays.  Off-by-1 random memory overwrites
 * drove me to distraction and forced the recoding of all the
 * arrays to start at 0 and the redefining of the macros below
 * to start at 0.  I'm glad I had all the loops coded as these
 * macros, because that helped quite a bit.
 */
/**************************************************************/

/**
 * \brief Loop through each group.  
 * 
 * Generates a GrpIndex to use when accessing RGroup (usually). 
 * 
 * \param c is a GrpIndex. Use it to access each element in RGroup
 *          like RGroup[c]
 * 
 * \ingroup RGROUP
 */
#define ForEachGroup(c)    for((c)=0; (c)< Globals->grpCount; (c)++)

/**
 * \brief Loop through all species regardless of group.
 * 
 * This function is useful for performing an operation on every
 * species.
 * 
 * \param s will increment every iteration of the loop. Use it to
 *          access every element in Species like Species[s].
 * 
 * \sa Species
 * 
 * \ingroup SPECIES
 */
#define ForEachSpecies(s)    for((s)=0; (s)< Globals->sppCount; (s)++)

/**
 * \brief Traverses a species' linked list of individuals.
 * 
 * \param i is a IndivType* which will point to a different
 *          individual every iteration of the loop.
 * \param s is a SpeciesType* which points to the desired
 *          SpeciesType. Must be set before calling this macro.
 * 
 * \ingroup INDIVIDUAL 
 */
#define ForEachIndiv(i,s)  for((i)=(s)->IndvHead; (i)!=NULL; (i)=(i)->Next)

/**
 * \brief Macro that loops through all established species in
 *        the designated rgroup.
 * 
 * Established species are a bit complicated because we
 * have to pull the species number from the RGroup->est_spp
 * array so we need three parameters.  
 * 
 * \param s is an SppIndex
 * \param g is a GrpIndex
 * \param i is an int. 
 * 
 * Note that s and i are changed.
 * s becomes the species number of each established species
 * but isn't required to be used in the code (and this can
 * give compiler warnings to that effect--ignore them).
 * i must be a declared lvalue but is only used internally to
 * the loop.  it is the index into the list of all species
 * that are actually established.
 * 
 * \ingroup SPECIES
 */
#define ForEachEstSpp(s,g,i) for((i)=0,(s)=RGroup[g]->est_spp[i];\
                                  (i) < RGroup[g]->est_count;    \
                                  (s)=RGroup[g]->est_spp[++(i)])

/**
 * \brief Macro that loops n times where n is the number of
 *        established species.
 * 
 * \param g is a GrpIndex that should be set to the desired
 *          group.
 * \param i is an int that will increment by 1 every loop of
 *          the program.
 * 
 * \ingroup SPECIES
 */
#define ForEachEstSpp2(g,i) for((i)=0; (i) < RGroup[g]->est_count; (i)++)

/**
 * \brief Loop over each possible species in the group established or not.
 * 
 * \param s should be a SppIndex. This is the value that is modified. it
 *          does not have to be predefined.
 * \param g should be a GrpIndex. This should be set to the desired group
 *          BEFORE calling the macro.
 * \param i should be an int. This does not have to be predefined.
 * 
 * Works exactly like ForEachEstSpp() but looks up a group's
 * possible species instead of established species.
 * 
 * \sa ForEachEstSpp() 
 * 
 * \ingroup RGROUP
 */
#define ForEachGroupSpp(s,g,i) for((i)=0,(s)=RGroup[g]->species[i];\
                                   (i) < RGroup[g]->max_spp;    \
                                   (s) = RGroup[g]->species[++(i)])

/**
 * \brief A compatability macro left behind from a time when STEPWAT2
 *        didn't explicitly store max age.
 * 
 * \param s should be a SppIndex
 * 
 * This is not to be used by future developers. Instead use Species[s]->max_age
 * 
 * \sa SppIndex
 * 
 * \ingroup SPECIES
 */
#define SppMaxAge(s)      (Species[s]->max_age)

/**
 * \brief A compatability macro left behind from a time when STEPWAT2
 *        didn't explicitly store max age.
 * 
 * \param g should be a GrpIndex.
 * 
 * This is not to be used by future developers. Instead use RGroup[g]->max_age.
 * 
 * \sa GrpIndex
 * 
 * \ingroup RGROUP
 */
#define GrpMaxAge(g)      (RGroup[g]->max_age)
/**************************************************************/


/* now define the structures (via ST_structs.h)
 * and their monikers
 */
#include "ST_structs.h"
typedef struct indiv_st IndivType;
typedef struct species_st SpeciesType;
typedef struct resourcegroup_st GroupType;
typedef struct succulent_st SucculentType;
typedef struct environs_st EnvType;
typedef struct plot_st PlotType;
typedef struct globals_st ModelType;
typedef struct bmassflags_st BmassFlagsType;
typedef struct mortflags_st MortFlagsType;
typedef struct superglobals_st GlobalType;


#endif
