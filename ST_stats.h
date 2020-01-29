/**
 * \file ST_stats.h 
 * \brief Defines all structs used in ST_stats.c
 * 
 * It was created to allow the gridded code to instantiate accumulators
 * without including ST_stats.c. 
 * 
 * To read the documentation related to creating this file see issues
 * #259 and #266 on GitHub 
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup STATISTICS
 */

#ifndef STATS_STRUCT_DEF
#define STATS_STRUCT_DEF

#include "ST_defines.h"

/* Basic struct that holds average, sum of differences squared, standard deviation
 * and number of entries. This is enough information to calculate running values without
 * storing raw data. If you would like to create a new accumulator do not use this
 * struct. Instead, use StatType.*/
struct accumulators_st {
  double ave, sum_dif_sqr, sd;
  unsigned long nobs;
};

/* Accumulator along with the RGroup or Species name. */
struct stat_st {
  char *name; /* array of ptrs to names in RGroup & Species */
  struct accumulators_st *s;
} typedef StatType;

/* Struct for wildfire and prescribed fire stats. */
struct fire_st {
  int *wildfire;
  int **prescribedFire;
} typedef FireStatsType;

/*----------------------- Exported Functions ---------------------------------- */
void stat_Collect( Int year ) ;
void stat_Collect_GMort ( void ) ;
void stat_Collect_SMort ( void ) ;
void stat_Output_YrMorts( void ) ;
void stat_Output_AllMorts( void) ;
void stat_Output_AllBmass(void) ;
void stat_free_mem( void );
void stat_Copy_Accumulators(StatType* newDist, StatType* newPpt, StatType* newTemp, StatType* newGrp, 
                            StatType* newGsize, StatType* newGpr, StatType* newGmort, StatType* newGestab, 
                            StatType* newSpp, StatType* newIndv, StatType* newSmort, StatType* newSestab, 
                            StatType* newSrecieved, FireStatsType* newGwf, Bool firstTime);
void make_header( char *buf);
void make_header_with_std( char *buf);

#endif