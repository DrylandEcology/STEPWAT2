/**
 * \file ST_seedDispersal.c
 * \brief Function definitions for all seed dispersal specific functions.
 * 
 * 
 * \author Chandler Haukap
 * \date 17 December 2019
 * 
 * \ingroup SEED_DISPERSAL_PRIVATE 
 */

#include "ST_globals.h"
#include "ST_defines.h"
#include "ST_grid.h"
#include "ST_seedDispersal.h"
#include "rands.h"
#include "myMemory.h"

float _distance(int row1, int row2, int col1, int col2, float cellLen);
Bool _shouldProduceSeeds(SppIndex sp);
float _rateOfDispersal(float PMD, float meanHeight, float maxHeight);
float _probabilityOfDispersal(float rate, float height, float distance);

/* RNG unique to seed dispersal. */
pcg32_random_t dispersal_rng;

/**
 * \brief Allocates any memory necessary for gridded mode.
 * 
 * \ingroup SEED_DISPERSAL
 */
void initDispersalParameters(void) {
  int senderRow, senderCol, receiverRow, receiverCol, MAXDP, maxCells;
  SppIndex sp;
  double maxRate;   /* Dispersability of the seeds */
  double plotWidth; /* width of the plots (all plots are equal and square) */
  double
      distanceBetweenPlots; /* distance between the sender and the receiver */
  CellType *sender;

  RandSeed(SuperGlobals.randseed, &dispersal_rng);

  /* sender denotes that these loops refer to the cell distributing seeds */
  for (senderRow = 0; senderCol < grid_Rows; ++senderRow) {
    for (senderCol = 0; senderCol < grid_Cols; ++senderCol) {
      /* Cell is loaded to ensure the global Species pointer points to a valid
         SpeciesType so the ForEachSpecies loop is able to iterate. */
      load_cell(senderRow, senderCol);
      sender = &gridCells[senderRow][senderCol];

      /* Allocate seed dispersal information. */
      sender->mySeedDispersal = Mem_Calloc(MAX_SPECIES, sizeof(Grid_SD_St),
                                           "initDispersalParameters");

      ForEachSpecies(sp) {
        if (!(Species[sp]->use_me && Species[sp]->use_dispersal)) {
          continue;
        }

        /* These are the three values we need to calculate the probability of
         * dispersal according to EQ 5 in Coffin & Lauenroth 1989. */
        //
        // MAXD = ((Species[sp]->sd_H * Species[sp]->sd_VW) /
        // Species[sp]->sd_VT) / 100.0; // divided by 100 to convert from cm to
        // m.
        maxRate = -(log(Species[sp]->maxDispersalProbability) /
                    Species[sp]->maxDispersalDistance);
        plotWidth = sqrt(Globals->plotsize);
        MAXDP = (int)ceil(Species[sp]->maxDispersalDistance /
                          plotWidth); // MAXD in terms of plots... rounds up to
                                      // the nearest integer
        maxCells = (int)pow((MAXDP * 2) + 1.0, 2.0);

        /* Allocate the dispersalProb 2d array */
        sender->mySeedDispersal[sp].dispersalProb =
            Mem_Calloc(grid_Rows, sizeof(double *),
                       "initDispersalParameters: dispersalProb");
        for (receiverRow = 0; receiverRow < grid_Rows; ++receiverRow) {
          sender->mySeedDispersal[sp].dispersalProb[receiverRow] =
              Mem_Calloc(grid_Cols, sizeof(double),
                         "initDispersalParameters: dispersalProb[i]");
        }

        /* Loop through all possible recipients of seeds. */
        for (receiverRow = 0; receiverRow < grid_Rows; ++receiverRow) {
          for (receiverCol = 0; receiverCol < grid_Cols; ++receiverCol) {
            if (senderRow == receiverRow && senderCol == receiverCol) {
              continue; // No need to calculate a probability for dispersal to
                        // itself
            }

            distanceBetweenPlots = _distance(senderRow, receiverRow, senderCol,
                                              receiverCol, plotWidth);

            /* The value that we are after should be saved to the sender cell.
             * this equation comes directly from equation 4 in Coffin and
             * Lauenroth 1989. */
            sender->mySeedDispersal[sp]
                .dispersalProb[receiverRow][receiverCol] =
                (distanceBetweenPlots > Species[sp]->maxDispersalDistance)
                    ? (0.0)
                    : (exp(-maxRate * distanceBetweenPlots));
          }
        }
      }
    }
  }
  unload_cell();
}

/**
 * \brief Disperse seeds between cells.
 * 
 * Iterates through all senders and recipients and deterimes which plots 
 * received seeds. If a cell does receive seeds for a given species,
 * Species[sp]->seedsPresent will be set to TRUE.
 * 
 * \sideeffect 
 *    For all cells and all species Species[sp]->seedsPresent will be set to 
 *    TRUE if a seed reached the cell and FALSE if no seeds reached the cell.
 *
 * \author Chandler Haukap
 * \date 18 December 2019
 *
 * \ingroup SEED_DISPERSAL
 */
void disperseSeeds(void) {
  SppIndex sp;
  CellType *receiverCell;
  int row, col;
  int receiverRow, receiverCol;
  // The probability of dispersal
  double Pd;
  double rate;
  double height;
  double distance;

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
        // If this species can't produce we are done.
        if (!_shouldProduceSeeds(sp))
          continue;
        // Running this algorithm on Species that didn't request dispersal
        // wouldn't hurt, but it would be a waste of time.
        if (!Species[sp]->use_dispersal){
          continue;
        }

        // These variables are independent of recipient.
        height = getSpeciesHeight(Species[sp]);
        rate = _rateOfDispersal(Species[sp]->maxDispersalProbability,
                                Species[sp]->meanHeight, 
                                Species[sp]->maxHeight);

        // Iterate through all possible recipients of seeds.
        for (receiverRow = 0; receiverRow < grid_Rows; ++receiverRow) {
          for (receiverCol = 0; receiverCol < grid_Cols; ++receiverCol) {
            receiverCell = &gridCells[receiverRow][receiverCol];
            // This algorithm wouldn't hurt anything but it would waste time.
            if (!receiverCell->mySpecies[sp]->use_dispersal){
              continue;
            }

            // If this cell already has seeds theres no point continuing.
            if(receiverCell->mySpecies[sp]->seedsPresent){
              continue;
            }

            // These variables depend on the recipient.
            distance = _distance(row, receiverRow, col, receiverCol,
                                 Globals->plotsize);
            Pd = _probabilityOfDispersal(rate, height, distance);

            // Stochastically determine if seeds reached the recipient.
            if (RandUni(&dispersal_rng) < Pd) {
              // Remember that Species[sp] refers to the sender, but in this
              // case we are refering to the receiver.
              receiverCell->mySpecies[sp]->seedsPresent = TRUE;
            }
          } // END for each receiverCol
        } // END for each receiverRow
      } // END ForEachSpecies(sp)
      unload_cell();
    } // END for each col
  } // END for each row
}

/**
 * \brief Calculates the distance betwen two 2-dimensional points.
 *
 * \param row1 The row of the first cell, i.e. it's x coordinate.
 * \param row2 The row of the second cell, i.e. it's x coordinate.
 * \param col1 The column of the first cell, i.e. it's y coordinate.
 * \param col1 The column of the first cell, i.e. it's y coordinate.
 * \param cellLen The length of the square cells.
 *
 * \return A Double. The distance between the cells.
 *
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _distance(int row1, int row2, int col1, int col2, float cellLen) {
  double rowDist = row1 - row2;
  double colDist = col1 - col2;

  // returns the distance between the two grid cells
  if (row1 == row2) {
    return (abs(colDist) * cellLen);
  } else if (col1 == col2) {
    return (abs(rowDist) * cellLen);
  } else { // row1 != row2 && col1 != col2
    // Pythagorean theorem:
    return sqrt(pow(abs(colDist) * cellLen, 2.0) +
                pow(abs(rowDist) * cellLen, 2.0));
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
 * \param meanHeight is the average height of an individual of the given
 *                   species.
 * \param maxHeight is the maximum height of an individual of the given
 *                  species.
 *
 * \return A float.
 *
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _rateOfDispersal(float PMD, float meanHeight, float maxHeight) {
  return log(PMD) * meanHeight / maxHeight;
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
  return exp(rate / sqrt(height)) * distance;
}