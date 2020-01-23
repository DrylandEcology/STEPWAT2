/**
 * \file ST_colonization.c
 * \brief Function definitions and private variables for the 
 *        [colonization](\ref COLONIZATION) module.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION_PRIVATE
 */

#include "ST_colonization.h"
#include "myMemory.h"
#include "filefuncs.h"
#include "ST_functions.h"
#include "ST_globals.h"

#define MAX_COLONIZATION_EVENTS 50

/* --------------------------- Private Functions --------------------------- */
void _copyEvent(ColonizationEvent* dest, ColonizationEvent* src);

/**
 * \brief All of the events this module will simulate, sorted by their start
 *        year.
 * 
 * \ingroup COLONIZATION_PRIVATE
 */
ColonizationEvent* _allEvents = 0;

/**
 * \brief The number of [events](\ref ColonizationEvent) in the 
 *        \ref _allEvents array.
 * 
 * \ingroup COLONIZATION_PRIVATE
 */
int _numberOfEvents = 0;

/**
 * \brief The main function of the colonization module.
 * 
 * This function will perform all colonization events, meaning if an event is
 * supposed to occur in the given year it will provides seeds for the requested
 * species in the requested cell.
 * 
 * \return TRUE if colonization ocurred in any cell.
 * \return FALSE if nothing colonized in any cell.
 * 
 * \sideeffect
 *     If the species was turned off this function will also turn it on.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION
 */
Bool colonize(void) {
  Bool somethingColonized = FALSE;
  int i = 0;

  // While the startYear is earlier than the current year.
  while(_numberOfEvents > i && _allEvents[i].startYear <= Globals->currYear) {
    // If the event is occuring this year
    if(_allEvents[i].startYear + _allEvents[i].duration > Globals->currYear) {
      CellType* cell = &gridCells[_allEvents[i].row][_allEvents[i].column];
      // Setting both of these things to TRUE should ensure the species has a
      // chance to establish.
      cell->mySpecies[_allEvents[i].species]->seedsPresent = TRUE;
      cell->mySpecies[_allEvents[i].species]->use_me = TRUE;
      somethingColonized = TRUE;
    }
    ++i;
  }

  return somethingColonized;
}

/**
 * \brief Initialize the colonization module from a file.
 * 
 * This function will allocate memory for every 
 * [colonization event](\ref ColonizationEvent) the user specified in the input
 * file.
 * 
 * It will then populate the \ref _allEvents array with the events, sorted from
 * earliest startYear to latest.
 * 
 * \param fileName The name of the file to open. This can be a relative path
 *                 or an absolute path.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION
 */
void initColonization(char* fileName) {
  char inbuf[128];
  char name[5];
  int valuesRead;
  int numOfEvents = 0;
  int cell;
  int startYear;
  int duration;
  int i, j;

  // Open the file.
  FILE* file = fopen(fileName, "r");
  if(!file) {
    LogError(logfp, LOGERROR, "Error in Colonization module:"
             "\n\tCannot open input file %s", fileName);
  }

  // A temporary variable to store colonization events
  ColonizationEvent* tempEvents = Mem_Calloc(MAX_COLONIZATION_EVENTS,
                                             sizeof(ColonizationEvent),
                                             "initColonization");

  // Throw away the header
  GetALine(file, inbuf);

  // Read the entire file.
  while(GetALine(file, inbuf)) {
    valuesRead = sscanf(inbuf, "%d,%s,%d,%d", &cell, name, &startYear, &duration);
        
    /* ----------------------- Input Validation ------------------------ */
    if(numOfEvents >= MAX_COLONIZATION_EVENTS) {
      LogError(logfp, LOGWARN, "Error reading colonization file:"
               "A maximum of %d colonization events can be specified. "
               "The rest of the events will be ignored.", 
               MAX_COLONIZATION_EVENTS);
      break;
    }
    if(valuesRead != 4) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tIncorrect number of input arguments.");
      return;
    }
    if(cell < 0 || cell > (grid_Rows * grid_Cols)) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tInvalid cell (%d) specified.", cell);
      return;
    }
    if(startYear < 1 || startYear > SuperGlobals.runModelYears) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tInvalid start year (%d) specified.", startYear);
      return;
    }
    if(Species_Name2Index(name) == -1) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tUnrecognized species name (%s).", name);
      return;
    }

    /* ------------- Add the event to our temporary array -------------- */
    tempEvents[numOfEvents].row = cell / grid_Cols;
	  tempEvents[numOfEvents].column = cell % grid_Cols;
    tempEvents[numOfEvents].duration = duration;
    tempEvents[numOfEvents].species = Species_Name2Index(name);
    tempEvents[numOfEvents].startYear = startYear;
    numOfEvents++;
  }

  // Ensure we never leak memory
  if(_allEvents) {
    Mem_Free(_allEvents);
    _allEvents = 0;
  }

  /* ------------ Allocate and assign the file level variables ----------- */
  // This algorithm sorts the events by startYear as is copies them over.
  // Sorting the entries now will save time later.
  int lowestYear, lowestYearIndex;
  Mem_Calloc(numOfEvents, sizeof(ColonizationEvent), "initColonization");
  for(i = 0; i < numOfEvents; ++i) {
    // I set the default to 32,767 meaning if the simulation is run for
    // 32,768 years there could be a bug. However, I doubt that would ever
    // happen.
    lowestYear = 32767;
    lowestYearIndex = 0;

    // Find the lowest start year.
    for(j = 0; j < numOfEvents; ++j){
      if(tempEvents[j].startYear < lowestYear) {
        lowestYear = tempEvents[j].startYear;
        lowestYearIndex = j;
      }
    }

    _copyEvent(&_allEvents[i], &tempEvents[lowestYearIndex]);
    // This ensures we never pick that index again.
    tempEvents[lowestYearIndex].startYear = 32767;
  }
  _numberOfEvents = numOfEvents;

  // Free up the memory we used.
  Mem_Free(tempEvents);
  fclose(file);
}

/**
 * \brief Free this module's memory
 * 
 * This function is safe to call more than once, but it only needs to be called
 * once per simulation.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION
 */
void freeColonizationMemory(void) {
  if(_allEvents){
    Mem_Free(_allEvents);
    _allEvents = 0;
  }
}

/**
 * \brief Copy one [event](\ref ColonizationEvent)'s information to an other.
 * 
 * Note that both events MUST be allocated prior to calling this function.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION_PRIVATE
 */
void _copyEvent(ColonizationEvent* dest, ColonizationEvent* src) {
  dest->row = src->row;
  dest->column = src->column;
  dest->duration = src->duration;
  dest->species = src->species;
  dest->startYear = src->startYear;
}