/********************************************************/
/********************************************************/
/*  Source file: ST_functions.h
/*  Type: header
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: Declares public functions used throughout the
 *           model, although some are only used in specific
 *           places.
/*  History:
/*     (6/15/2000) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

#ifndef FUNCTION_DEF
#define FUNCTION_DEF

#include "ST_defines.h"

void Env_Generate( void );

/* See steppe_main.c for declarations of the following

   The following declarations are to be used only in the noted
   source files because they only make sense there.
*/
/*  --- main.c ---- */


/* the rest of these functions might reasonably be called
 * from several locations.
 */
char *Parm_name(ST_FileIndex i);

IndivType **RGroup_GetIndivs( GrpIndex rg, const char sort, IntS *num);
GrpIndex RGroup_New( void);
void RGroup_Kill( GrpIndex rg);
GrpIndex RGroup_Name2Index (const char *name) ;
void RGroup_Update_Newsize( GrpIndex rg);
RealF RGroup_GetBiomass( GrpIndex rg) ;

RealF Species_GetBiomass (SppIndex sp);
void Species_Add_Indiv( SppIndex sp, Int new_indivs);
void Species_Update_Newsize( SppIndex sp, RealF newsize );
void save_annual_species_relsize(SppIndex sp, RealF newsize);
SppIndex Species_Name2Index (const char *name);
void Species_Kill (const SppIndex sp, int killType);
void Species_Proportion_Kill (const SppIndex sp, int killType, RealF proportionKilled );
void Species_Proportion_Recovery (const SppIndex sp, int killType, RealF proportionRecovery,RealF proportionKilled);
void Species_Proportion_Grazing(const SppIndex sp,RealF proportionGrazing );
void Species_Annual_Kill(const SppIndex sp, int killType);
IntS Species_NumEstablish( SppIndex sp);

void Indiv_SortSize( const byte sorttype,
                     const size_t n, IndivType **list);
int Indiv_CompSize_A( const void *key1, const void *key2);
int Indiv_CompSize_D( const void *key1, const void *key2);


#ifdef DEBUG_MEM
  void RGroup_SetMemoryRefs(void);
  void Species_SetMemoryRefs(void);
  void Parm_SetMemoryRefs(void);
  void Stat_SetMemoryRefs(void);
#endif

#endif
