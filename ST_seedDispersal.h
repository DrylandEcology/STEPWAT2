/**
 * \file ST_seedDispersal.h
 * \brief Declarations for all functions and structures exported from the seed
 *        dispersal module.
 * 
 * The seed dispersal module is intended to run in [gridded mode](\ref GRID).
 * It offers an alternative to traditional establishment by allowing 
 * establishment only when seeds reach a cell from a nearby cell.
 * 
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL
 */

#ifndef SEEDDISPERSAL_H
#define SEEDDISPERSAL_H

#include "ST_defines.h"

/**
 * \brief A struct for a single dispersal event.
 *
 * A linked list of these events can be used to output any statistics you could
 * want about seed dispersal.
 *
 * \author Chandler Haukap
 * \date 28 January 2020
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
typedef struct dispersal_event_st {
    int year;
    int iteration;
    int fromCell;
    int toCell;
    char name[5];
    struct dispersal_event_st* next;
} DispersalEvent;

/* =================================================== */
/*            Externed Global Variables                */
/* --------------------------------------------------- */
extern Bool UseSeedDispersal;
extern Bool recordDispersalEvents;
extern pcg32_random_t dispersal_rng;

/* =================================================== */
/*             Global Function Declarations            */
/* --------------------------------------------------- */
// See ST_seedDispersal.c for documentation of these functions.
void disperseSeeds(int year);
void outputDispersalEvents(char* filePrefix);
void freeDispersalMemory(void);

#endif
