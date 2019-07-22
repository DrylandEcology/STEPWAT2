/*********************************************************************************
 * ST_seedDispersal.h
 * 
 * Contains declarations for exported functions from the seed dispersal module.
 * It also contains the declaration for the UseSeedDispersal flag.
 *********************************************************************************/

#ifndef SEEDDISPERSAL_H
#define SEEDDISPERSAL_H

/* TRUE if we should run seed dispersal between years during the main simulation. */
Bool UseSeedDispersal;

void disperseSeeds(void);
void initDispersalParameters(void);

#endif