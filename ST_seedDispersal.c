/**
 * \file ST_seedDispersal.c
 * \brief Function definitions for all seed dispersal specific functions.
 *
 * \author Chandler Haukap
 * \date 17 December 2019
 *
 * \ingroup SEED_DISPERSAL_PRIVATE
 */

#include "ST_seedDispersal.h"
#include "ST_globals.h"
#include "ST_grid.h"
#include "myMemory.h"
#include "rands.h"

float _distance(int row1, int row2, int col1, int col2, float cellLen);
Bool _shouldProduceSeeds(SppIndex sp);
float _rateOfDispersal(float PMD, float meanHeight, float maxHeight);
float _probabilityOfDispersal(float rate, float height, float distance);

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
 * \sideeffect
 *    For all cells and all species Species[sp]->seedsPresent will be set to
 *    TRUE if a seed reached the cell and FALSE if no seeds reached the cell.
 *
 * \author Chandler Haukap
 * \date 18 December 2019
 * \ingroup SEED_DISPERSAL
 */
void disperseSeeds(void) {
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
        rate = _rateOfDispersal(Species[sp]->maxDispersalProbability,
                                Species[sp]->meanHeight,
                                Species[sp]->maxDispersalDistance);

        // Iterate through all possible recipients of seeds.
        for (receiverRow = 0; receiverRow < grid_Rows; ++receiverRow) {
          for (receiverCol = 0; receiverCol < grid_Cols; ++receiverCol) {
            receiverCell = &gridCells[receiverRow][receiverCol];

            // This algorithm wouldn't hurt anything but it would waste time.
            if (!receiverCell->mySpecies[sp]->use_dispersal) {
              continue;
            }

            // If this cell already has seeds there is no point in continuing
            if (receiverCell->mySpecies[sp]->seedsPresent) {
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
        }   // END for each receiverRow
      }     // END ForEachSpecies(sp)
      unload_cell();
    } // END for each col
  }   // END for each row
}

/**
 * \brief Seed dispersal initialization.
 *
 * This function allows species that requested seed dispersal establishment to
 * establish without requiring a mature specimen to distribute seeds.
 * Basically, calling this function simulates a seed bank for every species in
 * every cell that requested seed dispersal.
 *
 * This function should be called in place of \ref disperseSeeds when you want
 * to give every species that requested seed dispersal a chance to establish.
 *
 * \author Chandler Haukap
 * \date 16 January 2020
 * \ingroup SEED_DISPERSAL
 */
void initSeedDispersal(void) {
  int row, col;
  SppIndex sp;

  for (row = 0; row < grid_Rows; ++row) {
    for (col = 0; col < grid_Cols; ++col) {
      load_cell(row, col);
      ForEachSpecies(sp) {
        if (Species[sp]->use_dispersal) {
          Species[sp]->seedsPresent = TRUE;
        }
      } // END For Each Species
    }   // END For Each Column
  }     // END For Each Row
  unload_cell();
}

/**
 * \brief Calculates the distance between two 2-dimensional points.
 *
 * \param row1 The row of the first cell, i.e. its x coordinate.
 * \param row2 The row of the second cell, i.e. its x coordinate.
 * \param col1 The column of the first cell, i.e. its y coordinate.
 * \param col1 The column of the first cell, i.e. its y coordinate.
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
 * \param maxDistance is the maximum distance an individual of this species can
 *                    disperse seeds.
 *
 * \return A float.
 *
 * \author Chandler Haukap
 * \date 17 December 2019
 * \ingroup SEED_DISPERSAL_PRIVATE
 */
float _rateOfDispersal(float PMD, float meanHeight, float maxDistance) {
  return log(PMD) * meanHeight / maxDistance;
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