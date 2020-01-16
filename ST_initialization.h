/** 
 * \file ST_initialization.h
 * \brief Defines the exported functions from the 
 *        [initialization](\ref INITIALIZATION) module.
 * 
 * Initialization is currently only available for [gridded mode](\ref GRID).
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup INITIALIZATION
 */

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

/* This module only functions in the context of gridded mode. */
#include"ST_grid.h"

/**************** Enumerator for initialization types ******************/

/* Possible methods of initialization. */
typedef enum
{
	INIT_WITH_SPINUP,
	INIT_WITH_SEEDS,
	INIT_WITH_NOTHING
} InitializationMethod;

/********** Exported functions defined in ST_initialization.c ***********/

void runInitialization(void);
void loadInitializationConditions(void);
void freeInitializationMemory(void);

/************************ Exported variables ****************************/

/* Stores the state of the cells following spinup. */
CellType** initializationCells;
/* The method of initialization specified in inputs. */
InitializationMethod initializationMethod;
/* TRUE if the program is currently in initialization. */
Bool DuringInitialization;

#endif