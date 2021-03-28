/**
 * \file ST_seedDispersal.c
 * \brief Function definitions for all seed dispersal specific functions.
 * 
 * Note that this module uses the underscore prefix to denote "private"
 * functions and variables.
 *
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL_PRIVATE
 */

#include "ST_seedDispersal.h"
#include "ST_globals.h"
#include "ST_grid.h"
#include "sw_src/myMemory.h"
#include "sw_src/rands.h"

float _distance(int x1, int y1, int x2, int y2, float cellWidth);
Bool _shouldProduceSeeds(SppIndex sp);
float _rateOfDispersal(float PMD, float maxHeight, float maxDistance);
float _probabilityOfDispersal(float rate, float height, float distance);
float _maxDispersalDistance(float height);
void _recordDispersalEvent(int year, int iteration, int fromCell, int toCell, 
                           const char* name);

/**
 * \brief A struct for a single dispersal event. 
 * 
 * A linked list of these events can be used to output any statistics you could
 * want about seed dispersal.
 * 
 * \author Chandler Haukap
 * \date 28 January 2020
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
struct dispersal_event_st {
    int year;
    int iteration;
    int fromCell;
    int toCell;
    char name[5];
    struct dispersal_event_st* next;
} typedef DispersalEvent;

/**
 * \brief A module level variable pointing to the first \ref DispersalEvent.
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
DispersalEvent* _firstEvent = NULL;

/**
 * \brief A module level variable pointing to the last \ref DispersalEvent.
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
DispersalEvent* _lastEvent = NULL;

/**
 * \brief The random number generator for the seed dispersal module.
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
pcg32_random_t dispersal_rng;

/**
 * \brief TRUE if \ref dispersal_rng has already been seeded.
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
Bool isRNGSeeded = FALSE;

/**
 * \brief Disperse seeds between cells.
 *
 * Iterates through all senders and recipients and determines which cells
 * received seeds. If a cell does receive seeds for a given species,
 * Species[sp]->seedsPresent will be set to TRUE.
 * 
 * \param year the year of the simulation. This is input as a parameter because
 *             disperseSeeds() is typically called before the gridded mode
 *             updates the year for each cell.
 *
 * \sideeffect
 *    For all cells and all species Species[sp]->seedsPresent will be set to
 *    TRUE if a seed reached the cell and FALSE if no seeds reached the cell.
 *
 * \author Chandler Haukap
 * \date 18 December 2019
 * \ingroup SEED_DISPERSAL
 */
void disperseSeeds(int year) {
  SppIndex sp;
  CellType *receiverCell;
  int row, col;
  int receiverRow, receiverCol;
  // The probability of dispersal
  double Pd;
  // The rate of dispersal
  double rate;
  // The height of the tallest individual
  double height;
  // The distance between a potential sender and recipient
  double distance;

  if (!isRNGSeeded) {
    RandSeed(SuperGlobals.randseed, &dispersal_rng);
    isRNGSeeded = TRUE;
  }

  // Before we do anything we need to reset seedsPresent.
  for (row = 0; row < grid_Rows; ++row) {
    for (col = 0; col < grid_Cols; ++col) {
      load_cell(row, col);
      ForEachSpecies(sp) { Species[sp]->seedsPresent = FALSE; }
      unload_cell();
    }
  }

  for (row = 0; row < grid_Rows; ++row) {
    for (col = 0; col < grid_Cols; ++col) {
      load_cell(row, col);

      // This loop refers to the Species array of the SENDER.
      ForEachSpecies(sp) {
        // Running this algorithm on Species that didn't request dispersal
        // wouldn't hurt, but it would be a waste of time.
        if (!Species[sp]->use_dispersal)
          continue;

        // If there are no individuals of this species that are of reproductive
        // age continue.
        if (!_shouldProduceSeeds(sp))
          continue;

        // These variables are independent of recipient.
        height = getSpeciesHeight(Species[sp]);
		printf("height after getSpeciesHeight called Species = %s, height = %f\n ", Species[sp]->name, height);

        rate = _rateOfDispersal(Species[sp]->maxDispersalProbability,
                                Species[sp]->maxHeight,
                                _maxDispersalDistance(height));

		printf("rate after _rateOfDispersal called Species = %s, rate = %f\n ", Species[sp]->name, rate);

        // Iterate through all possible recipients of seeds.
        for (receiverRow = 0; receiverRow < grid_Rows; ++receiverRow) {
          for (receiverCol = 0; receiverCol < grid_Cols; ++receiverCol) {
            receiverCell = &gridCells[receiverRow][receiverCol];

            // This algorithm wouldn't hurt anything but it would waste time.
            if (!receiverCell->mySpecies[sp]->use_dispersal) {
              continue;
            }

            // If this cell already has seeds there is no point in continuing
            //if (receiverCell->mySpecies[sp]->seedsPresent) {
            //  continue;
            //}

            // These variables depend on the recipient.
            distance = _distance(col, row, receiverCol, receiverRow,
                                 Globals->plotsize);
            Pd = _probabilityOfDispersal(rate, height, distance);

            // Stochastically determine if seeds reached the recipient.
            if (RandUni(&dispersal_rng) < Pd) {
              // Remember that Species[sp] refers to the sender, but in this
              // case we are refering to the receiver.
              receiverCell->mySpecies[sp]->seedsPresent = TRUE;

              // If the user requested statistics.
              if(recordDispersalEvents) {
                _recordDispersalEvent(year, Globals->currIter, 
                                      (row * grid_Cols) + col, (receiverRow *
                                      grid_Cols) + receiverCol, 
                                      Species[sp]->name);
              }
            }
          } // END for each receiverCol
        }   // END for each receiverRow
      }     // END ForEachSpecies(sp)
      unload_cell();
    } // END for each col
  }   // END for each row
}

/**
 * \brief Output a summary of every [dispersal event](\ref DispersalEvent) that
 *        has occurred.
 * 
 * This function will output a file for every \ref gridCell. The files will 
 * contain one entry for every time the associated cell received seeds from 
 * another [cell](\ref CellType).
 * 
 * \param filePrefix is the name all of the files should have. The actual names
 *                   of the files will be "<filePrefix><N>.csv" where N is the
 *                   number of the [cell](\ref CellType) associated with the 
 *                   file.
 * 
 * \author Chandler Haukap
 * \date January 28 2020
 * \ingroup SEED_DISPERSAL
 */
void outputDispersalEvents(char* filePrefix) {
    char fileName[1024];
    int i;
    DispersalEvent* thisEvent = _firstEvent;
    FILE** files = Mem_Calloc(grid_Rows * grid_Cols, sizeof(FILE*), 
                              "outputDispersalEvents");

    for(i = 0; i < grid_Rows * grid_Cols; ++i) {
        sprintf(fileName, "%s%d.csv", filePrefix, i);
        files[i] = fopen(fileName, "w");
        fprintf(files[i], "Iteration,Year,From Cell,Species,To Cell\n");
    }

    while(thisEvent) {
        fprintf(files[thisEvent->toCell], "%d,%d,%d,%s,%d\n", 
                thisEvent->iteration, thisEvent->year, thisEvent->fromCell,
                thisEvent->name, thisEvent->toCell);
        
        thisEvent = thisEvent->next;
    }

    for(i = 0; i < grid_Rows * grid_Cols; ++i) {
        fclose(files[i]);
    }
    Mem_Free(files);
}

/**
 * \brief Free the memory allocated in the 
 *        [seed dispersal module](\ref SEED_DISPERSAL)
 * 
 * This function should be called after running the colonization module to
 * ensure that no memory leaks occur. It is safe to call multiple times.
 * 
 * \author Chandler Haukap
 * \date January 28 2020
 * \ingroup SEED_DISPERSAL
 */
void freeDispersalMemory(void) {
    DispersalEvent* thisEvent = _firstEvent;
    DispersalEvent* nextEvent;

    while(thisEvent != NULL){
        nextEvent = thisEvent->next;
        Mem_Free(thisEvent);
        thisEvent = nextEvent;
    }

    _firstEvent = NULL;
    _lastEvent = NULL;
}

/**
 * \brief Calculates the distance between two 2-dimensional points.
 *
 * \param x1 The column of the first cell, i.e. it's x coordinate.
 * \param y1 The row of the first cell, i.e. its y coordinate.
 * \param x2 The column of the second cell, i.e. its x coordinate.
 * \param y2 The row of the second cell, i.e. its y coordinate.
 * \param cellWidth The length of the square cells.
 *
 * \return A Double. The distance between the cells.
 *
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _distance(int x1, int y1, int x2, int y2, float cellWidth) {
  double rowDist = abs((x1 - x2)) * cellWidth;
  double colDist = abs((y1 - y2)) * cellWidth;

  // returns the distance between the two grid cells
  if (rowDist == 0) {
    // If rowdist == 0 these cells are in the same row. Might as well save some
    // computational power.
    return colDist;
  } else if (colDist == 0) {
    // If colDist == 0 these cells are in the same row. Might as well save some
    // computational power.
    return rowDist;
  } else {
    // Pythagorean theorem:
    return sqrt(pow(colDist, 2.0) + pow(rowDist, 2.0));
  }
}

/**
 * \brief Determines if a given species in the [loaded cell](\ref load_cell)
 *        is capable of producing seeds.
 *
 * Note that a cell must be loaded for this function to work.
 *
 * \param sp the index in the \ref Species array of the species to test.
 *
 * \return TRUE if there is a sexually mature individual of the given species.\n
 *         FALSE if there is not.
 *
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
Bool _shouldProduceSeeds(SppIndex sp) {
  IndivType *thisIndiv;
  SpeciesType *thisSpecies = Species[sp];

  ForEachIndiv(thisIndiv, thisSpecies) {
    if (thisIndiv->relsize >= thisSpecies->minReproductiveSize) {
      return TRUE;
    }
  }

  return FALSE;
}

/**
 * \brief Returns the rate of dispersal.
 *
 * \param PMD is the probability of maximum dispersal.
 * \param height is the average height of an individual of the given
 *                   species.
 * \param maxDistance is the maximum distance an individual of this species can
 *                    disperse seeds.
 *
 * \return A float.
 *
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _rateOfDispersal(float PMD, float maxHeight, float maxDistance) {
  return log((PMD) * maxHeight) / maxDistance;
}

/**
 * \brief Returns the probability that seeds will disperse a given distance.
 *
 * \param rate is the rate of seed dispersal.
 * \param height is the height of the tallest individual of the species.
 * \param distance is the distance the seeds must travel.
 *
 * \return A float.
 *
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _probabilityOfDispersal(float rate, float height, float distance) {
  return exp((rate * distance) / height);
}

/**
 * \brief Returns the maximum dispersal distance for a given 
 *        [individual](\ref IndivType).
 * 
 * This function is stochastic, meaning it will return a different value even if
 * it is given the same input parameter.
 * 
 * \param height The height, in cm, of the individual.
 * 
 * \return A float greater than or equal to 0.
 * 
 * \sa getSpeciesHeight which will give you the height of the tallest
 *     individual in a species.
 * 
 * \author Chandler Haukap
 * \date 27 February 2020
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _maxDispersalDistance(float height) {
    float constant = 8.337;
    return RandUniFloatRange(0, 2 * constant * (height/100), &dispersal_rng);
}

/**
 * \brief Add a [dispersal event](\ref DispersalEvent) to the linked list of
 *        dispersal events.
 * 
 * \param year is the year of the simulation when this event occurred.
 * \param fromCell is the origin of the seeds.
 * \param toCell is the recipient of the seeds.
 * \param name is the name of the species.
 * 
 * \sideeffect
 *     This will allocate memory for a new DispersalEvent.
 * 
 * \author Chandler Haukap
 * \date 28 January 2020
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
void _recordDispersalEvent(int year, int iteration, int fromCell, int toCell,
                           const char* name) {
  // Allocate a new event.
  DispersalEvent* newEvent = Mem_Calloc(1, sizeof(DispersalEvent), 
                                        "_recordDispersalEvent");

  // Populate the struct that we just allocated.
  newEvent->next = NULL;
  newEvent->year = year;
  newEvent->iteration = iteration;
  newEvent->fromCell = fromCell;
  newEvent->toCell = toCell;
  newEvent->name[0] = name[0];
  newEvent->name[1] = name[1];
  newEvent->name[2] = name[2];
  newEvent->name[3] = name[3];
  newEvent->name[4] = '\0';
  
  // Add the event to the linked list.
  if(!_firstEvent){
    _firstEvent = newEvent;
    _lastEvent = newEvent;
  } else {
    _lastEvent->next = newEvent;
    _lastEvent = newEvent;
  }
}
