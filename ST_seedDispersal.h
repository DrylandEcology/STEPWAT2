/**
 * \file ST_seedDispersal.h
 * \brief Declarations for all functions and structures exported from the seed
 *        dispersal module.
 * 
 * The seed dispersal module is intended to run in gridded mode. It offers
 * an alternative to traditional establishment by allowing establishment only when 
 * seeds reach a cell from a nearby cell.
 * 
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL
 */

#ifndef SEEDDISPERSAL_H
#define SEEDDISPERSAL_H

/** 
 * \brief If TRUE then seed dispersal will be used. 
 * \ingroup SEED_DISPERSAL
 */
Bool UseSeedDispersal;

// See ST_seedDispersal.c for documentation of these functions.
void disperseSeeds(void);

#endif