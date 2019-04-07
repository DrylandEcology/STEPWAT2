/* Chandler Haukap- This file defines all structs used in ST_stats.c.
 * It was created to allow the gridded code to instanciate accumulators
 * without including ST_stats.c. 
 * 
 * To read the documentation related to creating this struct see issues
 * #259 and #266 on GitHub */

#ifndef STATS_STRUCT_DEF
#define STATS_STRUCT_DEF


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

typedef struct  {
  struct accumulators_st *dist, *temp, *ppt, **grp1, **gsize, **gpr2,  **gwf2,  **gpf2, 
  							**gmort, **gestab, **spp, **indv, **smort, **sestab, **sreceived;
} accumulators_grid_st;

/* Struct for wildfire and prescribed fire stats. */
struct fire_st {
  int *wildfire;
  int **prescribedFire;
} typedef FireStatsType;

// Structure for holding sum values of all the grid cells
struct accumulators_grid_cell_st {
  double sum, sum_std;
  unsigned long nobs;
};

struct stat_grid_cell_st {
  char *name; /* array of ptrs to names in RGroup & Species */
  struct accumulators_grid_cell_st *s;    /* array of holding all the years values */
} typedef GridStatsType;

#endif