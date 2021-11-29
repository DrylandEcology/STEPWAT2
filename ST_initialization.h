/***********************************************************************/
/* ST_initialization.h
    This header defines the functions exported from ST_initialization.c.
    Initialization is currently only availible when running gridded
    mode.

    Initial programming by Chandler Haukap in August of 2019. */
/***********************************************************************/

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

/* This module only functions in the context of gridded mode. */
#include "ST_grid.h"

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


/* =================================================== */
/*            Externed Global Variables                */
/* --------------------------------------------------- */
extern CellType** initializationCells;
extern InitializationMethod initializationMethod;
extern Bool DuringInitialization;

#endif
