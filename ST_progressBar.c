/**************************************************************************/
/* ST_progressBar.c
    Function definitions for a progress bar printed to the terminal.
    See ST_progressBar.h for a description of how to add a new Status.
 */
/**************************************************************************/

#include "ST_progressBar.h"
#include "ST_defines.h"
#include "ST_globals.h"
#include<string.h>

/*************** Local Function(s). Treat these as private. ***************/

double _calculateProgress(int year, int iteration, Status status);
double _calculateInitializationProgress(int year);
double _calculateSimulationProgress(int year, int iteration);

/*********************** Function Definitions *****************************/

/* Log the program's progress using a progress bar.
	Param iteration: integer greater than 0. Input 0 iff the program is not currently in an iteration loop.
	Param year: integer greater than 0. Input 0 iff the program is not currently in a years loop.
	Param status: Use the "status" enum to choose a value.  */
void logProgress(int iteration, int year, Status status){
	static char progressString[256];
	int index = 0;					// Where we are in progressString
	Bool needsProgressBar = FALSE;	// By default we donot need a progress bar
	progressString[0] = '\0';		// Empty the string
	iteration--;					// iteration loops are 1 indexed, but we need 0 indexing.

	switch(status){
		case INITIALIZATION:
			strcpy(progressString, "Initializing |");
			index += 14;	// We have copied over 16 characters
            needsProgressBar = TRUE;
			break;
		case SIMULATION:
			strcpy(progressString, "Simulating   |");
			index += 14;	// We have copied over 12 characters
			needsProgressBar = TRUE;
			break;
		case OUTPUT:
			strcpy(progressString, "Writing files");
			index += 13;	// We have copied over 12 characters
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

/* Returns a double between 0 and 100 representing how close the program is to completing a given loop.
     Param year: Current year in the loop.
	 Param iteration: Current iteration in the loop.
	 Param status: Use the "status" enumerator. Valid options are SPINUP or SIMULATION. */
double _calculateProgress(int year, int iteration, Status status){
	double percentComplete;
	if(status == INITIALIZATION){
		percentComplete = _calculateInitializationProgress(year);
	} else if(status == SIMULATION) {
		percentComplete = _calculateSimulationProgress(year, iteration);
	} else {
		return 0;	// No other Status has defined how to calculate progress.
	}
	return percentComplete;
}

/* Algorithm for calculating how far along initialization is.
   Returns a percentage between 0 and 100. */
double _calculateInitializationProgress(int year){
    return (year / (double) SuperGlobals.runInitializationYears) * 100;
}

/* Algorithm for calculating how far along the simulation is.
   Returns a percentage between 0 and 100. */
double _calculateSimulationProgress(int year, int iteration){
    double prog = ((iteration * SuperGlobals.runModelYears) + year) 
						  / (double) (SuperGlobals.runModelIterations * SuperGlobals.runModelYears);

    return prog * 100;
}