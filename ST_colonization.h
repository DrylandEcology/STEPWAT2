/**
 * \file ST_colonization.h
 * 
 * Structs and functions exported from the [colonization](\ref COLONIZATION) 
 * module
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION
 */

#ifndef COLONIZATION_H
#define COLONIZATION_H

#include "ST_grid.h"

/* ------------------------------- Functions ------------------------------- */
// See ST_colonization.c for a description of these functions.
Bool colonize(int year);
void initColonization(char* fileName);
void freeColonizationMemory(void);

/* -------------------------------- structs -------------------------------- */
/**
 * \brief A struct to represent a single colonization event.
 * 
 * \ingroup COLONIZATION
 */
struct colonizationEvent_st {
    /** \brief Number of the first [cell](\ref CellType) to be colonized. */
    int fromCell;
    /** \brief Number of the last [cell](\ref CellType) to be colonized. */
    int toCell;
    /** \brief Year to begin colonization. */
    int startYear;
    /** \brief How many years colonization lasts. */
    int duration;
    /** \brief The index in \ref Species of the colonizing species. */
    SppIndex species;
} typedef ColonizationEvent;

#endif