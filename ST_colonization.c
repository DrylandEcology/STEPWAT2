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
#include "sw_src/myMemory.h"
#include "sw_src/filefuncs.h"
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
 * \param year is the current year. Year is input as a parameter, rather than 
 *             infered from a cell, to make sure that there is no descrepency 
 *             between the current year of any [cell](\ref CellType).
 * 
 * \return TRUE if colonization occurred in any cell.
 * \return FALSE if nothing colonized in any cell.
 * 
 * \sideeffect
 *     If the species was turned off this function will also turn it on.
 * 
 * \author Chandler Haukap
 * \date January 2020
 * \ingroup COLONIZATION
 */
Bool colonize(int year) {
  Bool somethingColonized = FALSE;
  ColonizationEvent* event;
  int i = 0, cell;
  
  // While the startYear is earlier than the current year.
  while(i < _numberOfEvents && _allEvents[i].startYear <= year) {
    // If the event is occurring this year
    if(_allEvents[i].startYear + _allEvents[i].duration > year) {
      event = &_allEvents[i];

      for(cell = event->fromCell; cell <= event->toCell; ++cell){
          load_cell(cell / grid_Cols, cell % grid_Cols);
          // Setting both of these things to TRUE should ensure the species has a
          // chance to establish.
          Species[event->species]->seedsPresent = TRUE;
          Species[event->species]->use_me = TRUE;
          Species[event->species]->use_dispersal = TRUE;

          unload_cell();
      }
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
  int fromCell, toCell;
  int startYear;
  int duration;
  int i, j;
  _numberOfEvents = 0;

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
    // Assume the event includes multiple cells.
    valuesRead = sscanf(inbuf, "%d-%d,%d,%d,%s", &fromCell, &toCell,&startYear,
                        &duration, name);
    if(valuesRead == 5) {
      // if fromCell > toCell, we will swap them. For example "8-4" --> "4-8".
      if(toCell < fromCell) {
        int temp = toCell;
        toCell = fromCell;
        fromCell = temp;
      }
    } else {
      // If the event only occurs in a single cell.
      valuesRead = sscanf(inbuf, "%d,%d,%d,%s", &fromCell, &startYear, 
                          &duration, name);
      // -1 means all cells.
      if(fromCell == -1){
        fromCell = 0;
        toCell = (grid_Rows * grid_Cols) - 1;

      // Otherwise this event must occur in only 1 cell.
      } else {
        toCell = fromCell;
      }
    }
        
    /* ----------------------- Input Validation ------------------------ */
    if(_numberOfEvents >= MAX_COLONIZATION_EVENTS) {
      LogError(logfp, LOGWARN, "Error reading colonization file:"
               "A maximum of %d colonization events can be specified. "
               "The rest of the events will be ignored.", 
               MAX_COLONIZATION_EVENTS);
      break;
    }
    if(valuesRead != 4 && valuesRead != 5) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:\n\tIncorrect"
               " number of input arguments.\n\tIs it possible you put a space"
               " somewhere? Remember this file cannot contain spaces.");
      return;
    }
    if(fromCell < 0 || toCell > ((grid_Rows * grid_Cols) - 1)) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tInvalid cell range ( %d - %d ) specified.", fromCell, toCell);
      return;
    }
    if(startYear < 1 || startYear > SuperGlobals.runModelYears) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tInvalid start year (%d) specified.", startYear);
      return;
    }
    // For the species info we need to have a cell loaded.
    load_cell(toCell / grid_Cols, toCell % grid_Cols);
    if(Species_Name2Index(name) == -1) {
      Mem_Free(tempEvents);
      LogError(logfp, LOGFATAL, "Error reading colonization file:"
               "\n\tUnrecognized species name (%s).", name);
      return;
    }

    /* ------------- Add the event to our temporary array -------------- */
    tempEvents[_numberOfEvents].fromCell = fromCell;
    tempEvents[_numberOfEvents].toCell = toCell;
    tempEvents[_numberOfEvents].duration = duration;
    tempEvents[_numberOfEvents].species = Species_Name2Index(name);
    tempEvents[_numberOfEvents].startYear = startYear;
    _numberOfEvents++;
    unload_cell();
  }

  // Ensure we never leak memory
  if(_allEvents) {
    Mem_Free(_allEvents);
    _allEvents = 0;
  }

  // If there were no events specified we are done.
  if(_numberOfEvents == 0) {
    Mem_Free(tempEvents);
    fclose(file);
    return;
  }

  /* ------------ Allocate and assign the file level variables ----------- */
  // This algorithm sorts the events by startYear as is copies them over.
  // Sorting the entries now will save time later.
  int lowestYear, lowestYearIndex;
  _allEvents = Mem_Calloc(_numberOfEvents, sizeof(ColonizationEvent), 
                          "initColonization");
  for(i = 0; i < _numberOfEvents; ++i) {
    lowestYear = 0;
    lowestYearIndex = -1;

    // Find the lowest start year.
    for(j = 0; j < _numberOfEvents; ++j){
      // (If this is the lowest start year we've seen) && it hasn't been used.
      if((tempEvents[j].startYear < lowestYear || lowestYearIndex == -1) &&
         tempEvents[j].startYear != -1) {
        lowestYear = tempEvents[j].startYear;
        lowestYearIndex = j;
      }
    }
    _copyEvent(&_allEvents[i], &tempEvents[lowestYearIndex]);

    // Mark this event as processed.
    tempEvents[lowestYearIndex].startYear = -1;
  }

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
  dest->fromCell = src->fromCell;
  dest->toCell = src->toCell;
  dest->duration = src->duration;
  dest->species = src->species;
  dest->startYear = src->startYear;
}
