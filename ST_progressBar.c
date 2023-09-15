/**
 * \file ST_progressBar.c
 * \brief Function definitions for a progress bar printed to the terminal.
 * 
 * For a description of how to integrate the progress bar into new code,
 * see \ref ST_progressBar.h.
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup PROGRESS_BAR
 */

#include "ST_progressBar.h"
#include "ST_globals.h"
#include <string.h>

/*************** Local Function(s). Treat these as private. ***************/

double _calculateProgress(int innerLoopIteration, int outerLoopIteration, Status status);
double _calculateSpinupProgress(int year);
double _calculateSimulationProgress(int year, int iteration);

/*********************** Function Definitions *****************************/

/**
 * \brief Log the program's progress using a progress bar.
 * 
 * This function expects to be called inside a nested iterations and years
 * loop, such as the one used in \ref main().
 * 
 * You may use this progress bar in a single for loop by setting either 
 * iteration or year to 0. 
 * 
 * Each \ref Status parameter must define their own algorithm for processing
 * iteration and year. For example, see \ref _calculateSpinupProgress
 * or \ref _calculateSimulationProgress. For a detailed description of how to 
 * implement a new \ref Status, see \ref ST_progressBar.h.
 * 
 * \param iteration integer greater than 0. Input 0 if and only if the program
 *                  is not currently in an iteration loop.
 * \param year integer greater than 0. Input 0 if and only if the program is
 *             not currently in a years loop.
 * \param status Use the \ref Status enum to choose a value.
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup PROGRESS_BAR
 */
void logProgress(int iteration, int year, Status status){
	static char progressString[256];
	int index = 0;					// Where we are in progressString
	Bool needsProgressBar = FALSE;	// By default we do not need a progress bar
	progressString[0] = '\0';		// Empty the string
	iteration--;					// iteration loops are 1 indexed, but we need 0 indexing.

	switch(status){
		case SPINUP:
			strcpy(progressString, "Spinning up  |");
			index += 14;	// We have copied over 14 characters
            needsProgressBar = TRUE;
			break;
		case SIMULATION:
			strcpy(progressString, "Simulating   |");
			index += 14;	// We have copied over 14 characters
			needsProgressBar = TRUE;
			break;
		case OUTPUT:
			strcpy(progressString, "Writing files...                   ");
			index += 13;	// We have copied over 13 characters
			break;
		case DONE:
			strcpy(progressString, "Done");
			// We need to pad the string with spaces to make sure we overwrite any progress bars.
			for(index = 4; index < 40; ++index) progressString[index] = ' ';
			break;
		default:
			break;
	}

	// If our program is currently in a looping state we can make a progress bar using year and iteration.
	if(needsProgressBar){
		int numberOfCharacters = 0;
		double percentComplete = _calculateProgress(year, iteration, status);
		// Add '=' characters for every 5% complete.
		while(percentComplete > 0){
			progressString[index] = '=';
			index++;
			numberOfCharacters++;
			percentComplete -= 5;
		}
		// Add spaces until we hit 20 characters
		while(numberOfCharacters < 20){
			progressString[index] = ' ';
			numberOfCharacters++;
			index++;
		}
		progressString[index++] = '|';
		progressString[index++] = '\0';
	}

	printf("\r%s", progressString);	// print the string we generated
	fflush(stdout);					// Explicitly print the output.

	// If we are done we want to print a newline character so the terminal isn't appended to our "Done" string.
	if(status == DONE){
		printf("\n");
	}
}

/**
 * \brief Calculates the program's progress inside of a looping structure.
 *
 * \param innerLoopIteration is the iteration of the inner loop.
 * \param outerLoopIteration is the iteration of the outer loop.
 * \param status Use the "status" enumerator. Valid options are SPINUP or 
 *               SIMULATION. 
 * 
 * \return A double between 0 and 100 representing how close the program is to
 *         completing the given loop(s).
 * 
 * The \ref Status you choose determines how innerLoopIteration and 
 * outerLoopIteration will be used to calculate an overall percentage.
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup PROGRESS_BAR_PRIVATE
 */
double _calculateProgress(int innerLoopIteration, int outerLoopIteration, Status status){
	double percentComplete;
	if(status == SPINUP){
		percentComplete = _calculateSpinupProgress(innerLoopIteration);
	} else if(status == SIMULATION) {
		percentComplete = _calculateSimulationProgress(innerLoopIteration, outerLoopIteration);
	} else {
		return 0;	// No other Status has defined how to calculate progress.
	}
	return percentComplete;
}

/**
 * \brief Algorithm for calculating how far along spinup is.
 * 
 * This function is intended to be called by \ref _calculateProgress when
 * the SPINUP \ref Status is used. Note that while 
 * \ref _calculateProgress takes two loop parameters this function only needs
 * one, because [spinup](\ref SPINUP) always runs for 1
 * iteration.
 * 
 * \param year The current year in the "years" loop that spinup is
 *             running.
 * 
 * \return A double between 0 and 100 where 100 means "Spinup 
 *         complete". 
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup PROGRESS_BAR_PRIVATE
 */
double _calculateSpinupProgress(int year){
    return (year / (double) SuperGlobals.runSpinupYears) * 100;
}

/**
 * \brief Algorithm for calculating how far along the main simulation is.
 * 
 * The main simulation is expected to use 2 nested for loops: one for the
 * iteration and one for the year. This function therefore takes two parameters
 * and uses [SuperGlobals.runModelYears](\ref SuperGlobals) and 
 * [SuperGlobals.runModelIterations](\ref SuperGlobals) to determine a 
 * percentage.
 * 
 * \param year The current year of the inner loop
 * \param iteration The current iteration of the outer loop.
 * 
 * \return A double between 0 and 100 where 100 means "Simulation complete". 
 * 
 * \author Chandler Haukap
 * \date August 2019
 * \ingroup PROGRESS_BAR_PRIVATE
 */
double _calculateSimulationProgress(int year, int iteration){
    double prog = ((iteration * SuperGlobals.runModelYears) + year) 
						  / (double) (SuperGlobals.runModelIterations * SuperGlobals.runModelYears);

    return prog * 100;
}
