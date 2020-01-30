/** 
 * \file ST_spinup.h
 * \brief Defines the exported functions from the [spinup](\ref SPINUP) module.
 * 
 * Spinup is a way to allow species to establish and grow before the actual 
 * simulation begins. Basically, it performs the same tasks as \ref runGrid
 * except without [seed dispersal](\ref SEED_DISPERSAL), 
 * [grazing](\ref grazing_EndOfYear) or [fire](\ref mort_EndOfYear).
 * 
 * Spinup is currently only available for [gridded mode](\ref GRID).
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup SPINUP
 */

#ifndef SPINUP_H
#define SPINUP_H

/* This module only functions in the context of gridded mode. */
#include"ST_grid.h"

/********** Exported functions defined in ST_spinup.c ***********/

void runSpinup(void);
void loadSpinupConditions(void);
void freeSpinupMemory(void);

/************************ Exported variables ****************************/

/**
 * \brief Stores the state of the cells following spinup.
 * 
 * This information must be stored so it can be loaded in every iteration of
 * the spinup.
 * 
 * \ingroup SPINUP
 */
CellType** spinupCells;

/**
 * \brief TRUE if spinup should be run.
 * \ingroup SPINUP
 */
Bool shouldSpinup;

/**
 * \brief TRUE if the program is currently in spinup.
 * \ingroup SPINUP
 */
Bool DuringSpinup;

#endif