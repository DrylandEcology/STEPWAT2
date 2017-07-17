/********************************************************/
/********************************************************/
/*  Source file: ST_defines.h
/*  Type: header
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: Generic model definitions such as constants,
 *           enums, and looping contructs.
/*  History:
/*     (6/15/2000) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

#ifndef STEPPE_DEF_H
#define STEPPE_DEF_H

#include "generic.h"

/* see #include "ST_structs.h" below */

/***************************************************
 * Basic definitions
 ***************************************************/
#define MAX_YEARS            1000
#define MAX_INDIVS_PER_SPP   100
#define MAX_SPP_PER_GRP      10
#define MAX_RGROUPS          10
#define MAX_SPECIES (MAX_SPP_PER_GRP * MAX_RGROUPS)
#define MAX_INDIVS (MAX_INDIVS_PER_SPP * MAX_SPP_PER_GRP * MAX_RGROUPS)
#define MAX_INDIVS_PER_GRP (MAX_INDIVS_PER_SPP * MAX_SPP_PER_GRP)
#define MAX_GROUPNAMELEN     15
#define MAX_SPECIESNAMELEN   4
#define MAX_OUTFIELDS (MAX_SPECIES + (MAX_RGROUPS *2) + 5 +1 )
#define MAX_FIELDLEN MAX_GROUPNAMELEN + 6  /* +6 for xtra chars like _RSize, etc */
#define MAX_CELLS 10000 // defines the maximum number of cells in the grid option

/* Constants for flagging whether a sort is
   ascending or descending or none */
#define SORT_A 'a'
#define SORT_D 'd'
#define SORT_0 '\0'

/* define types for indices to data structures.
 * they're defined this way so they might be changed
 * (say to object pointers) in the future
 */
typedef IntS GrpIndex;
typedef IntS SppIndex;

/* FYI, enums start at 0 by default*/
/* I put a few "firsts" and "lasts" in here to allow
 * variable length enums where the code doesn't have
 * to know in advance, but I haven't used them much
 * yet.
 */
typedef enum {NoSeason=-1, CoolSeason, WarmSeason}
  TempClass;

typedef enum {Slow, NoResources, Intrinsic, Disturbance, LastMort}
  MortalityType;

typedef enum {NoDisturb, FecalPat, AntMound, Burrow, LastDisturb}
  DisturbEvent;

typedef enum {VerySensitive, Sensitive, Insensitive, VeryInsensitive}
  DisturbClass;

typedef enum {Ppt_Wet, Ppt_Norm, Ppt_Dry}
  PPTClass;

typedef enum {Intcpt, Slope, P0=0, P1, P2, P3, P4}
  Params;

typedef enum {DepthNonComp, DepthShallow, DepthMedium, DepthDeep, DepthLast}
  DepthClass;

typedef enum {F_First, F_Log, F_Model, F_Env, F_Plot, F_RGroup, F_Disturbance, F_Species,
              F_BMassFlag, F_BMassPre, F_BMassAvg,
              F_MortFlag,  F_MortPre,  F_MortAvg,
              F_SXW, F_EXE}
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

/* Basic loop through each group.  Generates a GrpIndex to use
 * when accessing RGroup (usually). 'c' must already be declared
 * preferably as GrpIndex, but note that it gets changed, so it
 * can only be an lvalue.  */
/* void ForEachGroup(GrpIndex) */
#define ForEachGroup(c)    for((c)=0; (c)< Globals.grpCount; (c)++)

/* Generate an SppIndex to access Species[] to loop over each
 * defined species, irrespective of group.   */
/* void ForEachSpecies(SppIndex) */
#define ForEachSpecies(s)    for((s)=0; (s)< Globals.sppCount; (s)++)

/* Same for individuals within a species. Traverses a species'
 * linked list of IndivType objects.  i is a pointer to
 * IndivType, s is a pointer to SpeciesType.  See the Indiv_*()
 * functions for more details. */
/* void ForEachIndiv(IndivType*, SpeciesType*) */
#define ForEachIndiv(i,s)  for((i)=(s)->IndvHead; (i)!=NULL; (i)=(i)->Next)

/* Established species are a bit more complicated because we
 * have to pull the species number from the RGroup->est_spp
 * array so we need three parameters.  s is SppIndex, g is
 * GrpIndex, i is int. Note that s and i are changed.
 * s becomes the species number of each established species
 * but isn't required to be used in the code (and this can
 * give compiler warnings to that effect--ignore them).
 * i must be a declared lvalue but is only used internally to
 * the loop.  it is the index into the list of all species
 * that is actually established. */
/* void ForEachEstSpp(SppIndex, GrpIndex, int) */
#define ForEachEstSpp(s,g,i) for((i)=0,(s)=RGroup[g]->est_spp[i];\
                                  (i) < RGroup[g]->est_count;    \
                                  (s)=RGroup[g]->est_spp[++(i)])
#define ForEachEstSpp2(g,i) for((i)=0; (i) < RGroup[g]->est_count; (i)++)

/* loop over each possible species in the group established or not.
 * Works exactly like ForEachEstSpp() but looks up a group's
 * possible species instead of established species.  */
/* void ForEachGroupSpp(SppIndex, GrpIndex, IntU) */
#define ForEachGroupSpp(s,g,i) for((i)=0,(s)=RGroup[g]->species[i];\
                                   (i) < RGroup[g]->max_spp;    \
                                   (s) = RGroup[g]->species[++(i)])

/* shorthand moniker.  used to be more complicated but now
 * it's just a convenience.  See parms_Check_Species() to
 * see how the max ages get set.  */
/* IntS SppMaxAge(SppIndex) */
#define SppMaxAge(s)      (Species[s]->max_age)

/* Another convenience mostly for consistency, but also in case
 * it seems better to make it more complicated in some future
 * version.
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


#endif
