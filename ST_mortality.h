/**
 * \file ST_mortality.h
 * \brief Functions and variables exported from the \ref MORTALITY module.
 * 
 * \author Chandler Haukap
 * \date February 12 2020
 * \ingroup MORTALITY
 */

#ifndef MORTALITY_H
#define MORTALITY_H

#include "generic.h"
#include "sw_src/pcg/pcg_basic.h"

/* --------------------------- Exported Structs ---------------------------- */

/**
 * \brief Information used when simulating the cheatgrass-wildfire loop.
 * 
 * The biggest determinant in cheatgrass driven wildfire is precipitation,
 * specifically precipitation in Spring and Winter. This struct stores the 
 * values from previous Spring and Winter precipitation as well as the running
 * averages of both.
 * 
 * Note that "year" in this context refers to the water year, which runs from
 * October to September.
 * 
 * \sa _updateCheatgrassPrecip, where this struct is updated each year.
 * \author Chandler Haukap
 * \date 13 January 2020
 * \ingroup MORTALITY
 */
struct CheatgrassPrecip_st {
  /** \brief The Spring precipitation in the previous 3 years.
   * The array is indexed from newest to oldest, meaning prevSpring[0] is the 
   * most recent value. */
  double prevSprings[3];
  /** \brief The precipitation in the last Winter, meaning the Oct-Dec
   * precipitation from 2 years ago and the Jan-Mar precipitation from last 
   * year. */
  double lastWinter;
  /** \brief The current year's Spring precipitation. */
  double currentSpring;

  /** \brief The running average of Spring precipitation. */
  double springMean;
  /** \brief The running average of Winter precipitation. */
  double winterMean;

  /** \brief The sum of October - December precipitation information from last
   *         year.
   * 
   * The variable we use to get precipitation information, \ref SXW, stores 
   * values for the current year, but we need the values from 2 years ago to 
   * calcualte last year's winter precipitation. We therefore store it here,
   * even though it does make the struct a little more confusing. 
   */ 
  double lastOctThruDec;
  /** \brief The sum of October - December precipitation information from the
   *         current year.
   * 
   * The variable we use to get precipitation information, \ref SXW, stores 
   * values for the current year, but we need the values from 2 years ago to 
   * calcualte last year's winter precipitation. We therefore store it here,
   * even though it does make the struct a little more confusing. 
   */ 
  double thisOctThruDec;
  /** \brief The sum of January - March precipitation information from last
   *         year.
   * 
   * The variable we use to get precipitation information, \ref SXW, stores 
   * values for the current year, but we need the values from 1 year ago to 
   * calculate last year's winter precipitation. We therefore store it here,
   * even though it does make the struct a little more confusing. 
   */ 
  double thisJanThruMar;
} typedef CheatgrassPrecip;

/* -------------------------- Exported Functions --------------------------- */
// See ST_mortality.c for definitions and documentation of these functions.

void mort_Main(Bool *killed);
void mort_EndOfYear(void);

void freeMortalityMemory(void);
void proportion_Recovery(void);
void grazing_EndOfYear(void);
void killAnnuals(void);
void killMaxage(void);
void killExtraGrowth(void);

void initCheatgrassPrecip(void);
void setCheatgrassPrecip(CheatgrassPrecip* newCheatgrassPrecip);
CheatgrassPrecip* getCheatgrassPrecip(void);

/* ----------------------------- Exported RNG ------------------------------ */

/**
 * \brief The random number generator specific to the 
 *        [mortality](\ref MORTALITY) module.
 * 
 * This RGN is declared in the header file for other modules can seed it. Other
 * modules should NOT use this RNG to generate random numbers.
 * 
 * \ingroup MORTALITY
 */
pcg32_random_t mortality_rng;

/* ---------------------------- Exported Flags ----------------------------- */

/**
 * \brief A flag for turning cheatgrass-driven wildfire on and off.
 * \ingroup MORTALITY
 */
Bool UseCheatgrassWildfire;

/* ---------------------------- Exported Enums ----------------------------- */

/**
 * \brief All types of mortality.
 * 
 * Used to record what killed an individual.
 * 
 * \sa indiv_st which instantiates this enumerator.
 * 
 * \ingroup MORTALITY
 */
typedef enum {
    Slow, 
    NoResources, 
    Intrinsic, 
    Disturbance, 
    LastMort
} MortalityType;

#endif