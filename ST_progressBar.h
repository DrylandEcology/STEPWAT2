/******************************************************************/
/* ST_progressBar.h
    Defines all exported objects from ST_progressBar.c. 

    TO ADD STATUSES:
        Adding a Status is easy. Start by adding your entry to the
        Status enum. Then, define what message your status should 
        print in the switch statement in logProgress() defined in
        ST_progressBar.c. Finally, if your Status requires a
        progress bar, create a function for calculating progress
        in ST_progressBar.c then add your function to 
        _calculateProgress().
*/
/******************************************************************/

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

/*********************** Enumerator(s) ****************************/

/* States of the program that the progress bar recognises. */
typedef enum 
{
	INITIALIZATION,
	SIMULATION,
	OUTPUT,
	DONE
} Status;

/******************** Exported Function(s) ************************/

void logProgress(int iteration, int year, Status status);

#endif