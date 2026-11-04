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
 * @brief Maximum number of bytes allocated for buffered seed availability
 * output records.
 */
#define SEED_AVAIL_BUFFER_BYTES (1024*1024)
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
    int seedN;
    char name[5];
    struct dispersal_event_st* next;
} DispersalEvent;

/**
 * @brief Stores a single seed availability output record.
 *
 * Contains seed availability information for one species in one grid cell
 * during a specific simulation iteration and year. Records are temporarily
 * stored in the seed availability output buffer before being written to the
 * CSV output file.
 */
typedef struct seed_availability_data_st {
    int iteration;          /**< Simulation iteration associated with the record. */
    int year;               /**< Simulation year associated with the record. */
    int cell;               /**< Grid cell identifier. */

    char name[5];           /**< Species name/code, limited to four characters plus
                                 the null terminator. */

    int seedsReceived;      /**< Number of seeds received by the species in the cell. */
    int seedsProduced;      /**< Number of seeds produced by the species in the cell. */
    int eind;               /**< Effective maximum number of individuals allowed to establish. */
    double pestab;          /**< Effective seedling establishment probability. */
} SeedAvailabilityData;

typedef struct seed_availability_st{
    SeedAvailabilityData *data; /**< Buffered seed availability records. */
    FILE *file;                 /**< Seed availability CSV output file. */
    size_t position;            /**< Next available position in the buffer. */
    size_t capacity;            /**< Maximum number of records in the buffer. */
    Bool headerSet;             /**< Whether the CSV header has been written. */
}SeedAvailability;

/* =================================================== */
/*            Externed Global Variables                */
/* --------------------------------------------------- */
extern Bool UseSeedDispersal;
extern Bool recordDispersalEvents;
extern sw_random_t dispersal_rng;
extern Bool outputSDData;
extern Bool outputSeedAvailability;

/* =================================================== */
/*             Global Function Declarations            */
/* --------------------------------------------------- */
// See ST_seedDispersal.c for documentation of these functions.
void disperseSeeds(int year);
void outputDispersalEvents(char* filePrefix);
void freeDispersalMemory(void);
void initDispersalMemory(void);

#endif
