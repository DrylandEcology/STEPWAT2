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
struct _grid_sd_struct
{ //for seed dispersal
	/* TRUE if seeds are present. */
	Bool seeds_present;
	/* TRUE if this cell has recieved any seeds. */
	Bool seeds_received;
	/* dispersalProb[row][col] = the probability that this cell will disperse seeds to cell (row,col). */
	double **dispersalProb;
}typedef Grid_SD_St;

/* TRUE if we should run seed dispersal between years during the main simulation. */
Bool UseSeedDispersal;

void disperseSeeds(void);
void initDispersalParameters(void);

#endif