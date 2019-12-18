/**
 * \file ST_seedDispersal.h
 * \brief Declarations for all functions and structures exported from the seed
 *        dispersal module.
 * 
 * The seed dispersal module is intended to run in gridded mode. It offers
 * an alternative to traditional establishment by allowing growth only when 
 * seeds reach a cell from a nearby cell.
 * 
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL
 */

#ifndef SEEDDISPERSAL_H
#define SEEDDISPERSAL_H

/**
 * \brief Holds seed dispersal information.
 * 
 * \ingroup SEED_DISPERSAL
 */
struct _grid_sd_struct
{ //for seed dispersal
	/* TRUE if seeds are present. */
	Bool seeds_present;
	/* TRUE if this cell has recieved any seeds. */
	Bool seeds_received;
	/* dispersalProb[row][col] = the probability that this cell will disperse seeds to cell (row,col). */
	double **dispersalProb;
}typedef Grid_SD_St;

/** 
 * \brief If TRUE then seed dispersal will be used. 
 * \ingroup SEED_DISPERSAL
 */
Bool UseSeedDispersal;

// See ST_seedDispersal.c for documentation of these functions.
void disperseSeeds(void);
void initDispersalParameters(void);

#endif