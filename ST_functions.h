
/**
 * \file ST_functions.h
 * \brief Declares public functions used throughout the model.
 * 
 * Ideally this file would be removed and all of the functions within would be
 * declared in module-specific header files. However, that is outside of the 
 * scope of my current issue.
 *  
 * \author CWB (initial programing)
 * \date 15 June 2000
 * \author Chandler Haukap (author of this documentation)
 * \date 16 January 2020
 * \ingroup STEPPE
 */

#ifndef FUNCTION_DEF
#define FUNCTION_DEF

#include "ST_defines.h"

void Env_Generate( void );
void copy_environment(const EnvType* src, EnvType* dest);
void copy_plot(const PlotType* src, PlotType* dest);
void copy_succulent(const SucculentType* src, SucculentType* dest);

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
GrpIndex RGroup_Name2Index (const char *name);
RealF getRGroupRelsize(GrpIndex rg);
void RGroup_Update_GrpResProp( GrpIndex rg);
RealF RGroup_GetBiomass( GrpIndex rg) ;

RealF Species_GetBiomass (SppIndex sp);
void Species_Add_Indiv( SppIndex sp, Int new_indivs);
RealF getSpeciesRelsize(SppIndex sp);
RealF getSpeciesHeight(SpeciesType* sp);
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
