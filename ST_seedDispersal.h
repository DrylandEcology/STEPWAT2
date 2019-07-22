

#ifndef SEEDDISPERSAL_H
#define SEEDDISPERSAL_H

/* TRUE if we should run seed dispersal between years during the main simulation. */
Bool UseSeedDispersal;

void disperseSeeds(void);
void initDispersalParameters(void);

#endif