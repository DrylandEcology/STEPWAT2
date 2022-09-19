/*********************************************************************************
 * ST_seedDispersal.h
 * 
 * Contains declarations for exported functions from the seed dispersal module.
 * It also contains the declaration for the UseSeedDispersal flag and the
 * seed dispersal struct.
 *********************************************************************************/

#ifndef SEEDDISPERSAL_H
#define SEEDDISPERSAL_H

/* Holds seed dispersal information. */
typedef struct _grid_sd_struct
{ //for seed dispersal
	/* TRUE if seeds are present. */
	Bool seeds_present;
	/* TRUE if this cell has recieved any seeds. */
	Bool seeds_received;
	/* dispersalProb[row][col] = the probability that this cell will disperse seeds to cell (row,col). */
	double **dispersalProb;
	/* Last year's precipitation. */
	double lyppt;
} Grid_SD_St;


/* =================================================== */
/*            Externed Global Variables                */
/* --------------------------------------------------- */
extern Bool UseSeedDispersal;
extern pcg32_random_t dispersal_rng;

/* =================================================== */
/*             Global Function Declarations            */
/* --------------------------------------------------- */
void disperseSeeds(void);
void initDispersalParameters(void);

#endif
