/**
 * \file ST_progressBar.h
 * \brief Defines all exported objects from ST_progressBar.c. 
 *
 * TO ADD A PROGRESS BAR TO NEW CODE:
 * 
 * Incorporating a progress bar into new code is easy. Start by adding your 
 * entry to the \ref Status enum. Then, define what message your status should
 * print in the switch statement in \ref logProgress. Finally, if your Status
 * requires a progress bar, create a function for calculating progress in 
 * \ref ST_progressBar.c (see \ref _calculateSimulationProgress for an example
 * of such a function) then add your function to the logic in 
 * \ref _calculateProgress().
 *
 * \author Chandler Haukap 
 * \date August 2019
 * \ingroup PROGRESS_BAR
 */

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

/*********************** Enumerator(s) ****************************/

/* States of the program that the progress bar recognizes. */
typedef enum 
{
	SPINUP,
	SIMULATION,
	OUTPUT,
	DONE
} Status;

/******************** Exported Function(s) ************************/

void logProgress(int iteration, int year, Status status);

#endif