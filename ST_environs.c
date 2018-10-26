/********************************************************/
/********************************************************/
/*  Source file: environs.c
 *  Type: module
 *  Application: STEPPE - plant community dynamics simulator
 *  Purpose: Controls all environmental phenomenon from
 *           creating ppt and temp to creating the
 *           disturbances. */
/*  History:
 *     (6/15/2000) -- INITIAL CODING - cwb */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <math.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/pcg/pcg_basic.h"


#include "rands.h"

#ifdef STEPWAT
  #include "sxw_funcs.h"
  #include "sw_src/filefuncs.h"
  #include "sw_src/Times.h"
  extern Bool UseSoilwat;
  extern SXW_t SXW;
#endif

/*********** Locally Used Function Declarations ************/
/***********************************************************/

static void _make_ppt( void);
static void _make_temp( void);
static void _set_ppt_reduction( void);
static void _set_temp_reduction( void);
static void _make_disturbance( void);


/**************************************************************/
/**************************************************************/

extern
  pcg32_random_t environs_rng;


void Env_Generate( void) {
/*======================================================*/
/* PURPOSE */
/* Wrapper to generate a new set of environmental factors,
   usually for the current year.  Any new environmental
   generators should be called from this subroutine.
*/
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*------------------------------------------------------*/

  Int rg;

  switch (UseSoilwat) {
    case FALSE:
         /* clear last year's leftover resources */
         ForEachGroup(rg) RGroup[rg]->res_avail = 0.0;
         break;
    case TRUE:
         SXW_Run_SOILWAT();
         break;
  }

  _make_ppt();
  _make_temp();
  _set_ppt_reduction();
  _set_temp_reduction();
  _make_disturbance();
}

/**************************************************************/
static void _make_ppt( void) {
/*======================================================*/
/* If not running SOILWAT,take a random number from normal distribution with*/
/* mean, stddev that is between min & max from */
/* the Globals.ppt structure.*/
/* Also set the growing season precip. */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/* cwb - 6-Dec-02 -- added code to interface with STEPWAT.
 *       The ppt and gsppt are set in _sxw_set_environs()
 *       but we still pass through this code to set the
 *       Dry/Wet/Normal state.
 * KAP 1/26/2017 The above note by CB is not correct. Env.ppt
 *      is set in _sxw_set_environs, but gsppt is not, it is
 *      set below. When using SOILWAT or not, the gsspt, ppt.dry,
 *      and ppt.wet is currently fixed each year and read from env.in
 *      We should consider calculating gsspt in the _sxw_set_environs
 *      function when running SOILWAT (so it is not fixed), and allowing
 *      what constitutes a wet and dry year to vary across sites. */

/*------------------------------------------------------*/

  IntS r=0, i;

#ifdef DEBUG_ENVCONST
  r=320;
#endif

  if (UseSoilwat)
  { // Run with SOILWAT2: we have monthly PPT and temperature to calculate
    // growing season precipitation as sum of monthly precipitation of those
    // months when mean air temperature exceeds a threshold `GROWING_BASE_TEMP`
    Env.gsppt = 0; // gsppt is defined as IntS and units are millimeters

    for (i = 0; i < MAX_MONTHS; i++)
    {
      Env.gsppt += GE(SXW.temp_monthly[i], GROWING_BASE_TEMP) ?
        (IntS) (SXW.ppt_monthly[i] * 10. + 0.5) : 0;
    }

    if (Env.gsppt <= 0)
    {
      LogError(logfp, LOGWARN, "Zero growing season precipitation in "\
        "year = %d of iteration = %d", Globals.currYear, Globals.currIter);
      Env.gsppt = 0;
    }

    /*
    printf("_make_ppt (year = %d of iteration = %d): "\
      "Env.gsppt(SOILWAT2) = %d would be gsppt(fixed proportion) = %.1f\n",
      Globals.currYear, Globals.currIter, Env.gsppt, Globals.gsppt_prop * Env.ppt);
    */

  } else {
    // run as STEPPE without SOILWAT2:
    while ( r < Globals.ppt.min || r > Globals.ppt.max )
      r = (IntS)(RandNorm( Globals.ppt.avg,  Globals.ppt.std, &environs_rng) +.5);
    if (Env.ppt > 0) {
      Env.lyppt = Env.ppt;
      Env.ppt = r;
    } else {
      Env.lyppt = Env.ppt = r;
    }

    Env.gsppt = (IntS) (Globals.gsppt_prop * Env.ppt);
  }

  if ( Env.ppt <= Globals.ppt.dry )
    Env.wet_dry = Ppt_Dry;
  else if (Env.ppt >= Globals.ppt.wet)
    Env.wet_dry = Ppt_Wet;
  else
    Env.wet_dry = Ppt_Norm;
}

/**************************************************************/
static void _make_temp( void) {
/*======================================================*/
/* take a random number from normal distribution with*/
/* mean, stddev, that is between min & max from */
/* the Globals.temp structure.*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/* cwb - 6-Dec-02 -- added interface to STEPWAT.  Temp is
 *       set in _sxw_set_environs(). */

/*------------------------------------------------------*/

  RealF r=0.;

#ifdef DEBUG_ENVCONST
  r=25.;
#endif

  if (!UseSoilwat) {
    while ( r < Globals.temp.min || r > Globals.temp.max )
      r = RandNorm(Globals.temp.avg, Globals.temp.std, &environs_rng);
    Env.temp = r;
  }
}

/**************************************************************/
static void _set_ppt_reduction( void) {
/*======================================================*/


/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

/* EQN 10*/
  Succulent.reduction = fabs(Succulent.growth[Slope]
                           * Env.gsppt
                           + Succulent.growth[Intcpt]);
/* EQN 16*/
  Succulent.prob_death = (Succulent.mort[Slope] * Env.gsppt
                          + Succulent.mort[Intcpt]) / 100.0;
}

/**************************************************************/
static void _set_temp_reduction( void) {
/*======================================================*/
/* This routine implements EQNs 12 and 13*/

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/*
 *  4-Mar-03 - ran into a bug where temperature extremes cause
 *             negative temp_reductions, so I truncate to 0.
 *  7-Mar-03 - found a logic bomb in C&L90!.  Eqns 12 & 13
 *             are reported to be based on daily temps but
 *             the model necessarily uses MAT.  I'm adding
 *             a hack to compensate somewhat until a better
 *             method comes along.  The current fix is to
 *             account for the fact that the MAT is actually
 *             about 1/3 of optimal max yearly temp. */

/*------------------------------------------------------*/
  int i;
  RealF tp[4]; /* parms for temp growth modifier eqn.*/

  for ( i=CoolSeason; i <= WarmSeason; i ++ ){
    tp[1] = Globals.tempparm[i][0];
    tp[2] = Globals.tempparm[i][1];
    tp[3] = Globals.tempparm[i][2];
    tp[0] = Env.temp + tp[1];
    Env.temp_reduction[i] = tp[2]*tp[0] + tp[3] * (tp[0]*tp[0]);
    Env.temp_reduction[i] = max(0., Env.temp_reduction[i]);
  }

#ifdef STEPWAT
  if (UseSoilwat) {
    if (Env.temp < 9.5 ) {
      Env.temp_reduction[CoolSeason] = .9;
      Env.temp_reduction[WarmSeason] = .6;
    } else {
      Env.temp_reduction[CoolSeason] = .6;
      Env.temp_reduction[WarmSeason] = .9;
    }
  }
#endif

}

/**************************************************************/
static void _make_disturbance( void) {
/*======================================================*/
/* PURPOSE */
/* Generate disturbances, if any, for this year. */
/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000 */

/*------------------------------------------------------*/
  RealF pc; /* probability of colonization if current */
            /* disturbance is fecalpat*/
 /*new disturbance event generated, if any*/
  DisturbEvent event;

/*------------------------------------------------------*/
  /* Can't have simultaneous disturbances*/
  if (Plot.disturbance != NoDisturb) {
    switch( Plot.disturbance) {
      case FecalPat:
           if (Plot.pat_removed) {
             Plot.disturbed = 0;
             Plot.pat_removed = FALSE;
             Plot.disturbance = NoDisturb;
           } else {
             pc = Globals.pat.recol[Slope] * Plot.disturbed
                  + Globals.pat.recol[Intcpt];
             if (RandUni(&environs_rng) <= pc) {
               Plot.pat_removed = TRUE;
               /* slight effects for one year*/
               Plot.disturbed = 1;
             } else {
               Plot.pat_removed = FALSE;
               Plot.disturbed++;
             }
           }
           break;
      case NoDisturb: // Does nothing but prevent a compiler warning
      case LastDisturb:
      default:
           Plot.disturbed = (Plot.disturbed) ? Plot.disturbed -1 : 0;
    }
    if (Plot.disturbed == 0)
      Plot.disturbance = NoDisturb;

  }

  /* if the disturbance was expired above, */
  /* we can generate a new one immediately */
  if (Plot.disturbance == NoDisturb) {

    /* pick some type of disturbance (other than none)*/
    event = (DisturbEvent) RandUniIntRange(1, LastDisturb -1, &environs_rng);

    /* make sure this is off unless needed  */
    Plot.pat_removed = FALSE;
    switch( event) {
      case FecalPat:
       if (!Globals.pat.use) {event=NoDisturb; break;}
         event = (RandUni(&environs_rng) <= Globals.pat.occur)
               ? event : NoDisturb;
         if (event == NoDisturb) break;
         Plot.pat_removed = (RandUni(&environs_rng) <= Globals.pat.removal)
                           ? TRUE : FALSE;
         Plot.disturbed = 0;
         break;
      case AntMound:
       if (!Globals.mound.use) {event=NoDisturb; break;}
         event = (RandUni(&environs_rng) <= Globals.mound.occur)
               ? event :NoDisturb;
         if (event == NoDisturb) break;

         Plot.disturbed = RandUniIntRange(Globals.mound.minyr,
                                       Globals.mound.maxyr,
                                       &environs_rng);
         break;
      case Burrow:
       if (!Globals.burrow.use) {event=NoDisturb; break;}
         event = (RandUni(&environs_rng) <= Globals.burrow.occur)
               ? event :NoDisturb;
         if (event == NoDisturb) break;

         Plot.disturbed = (Globals.burrow.minyr > 0)
                        ? RandUniIntRange(1, Globals.burrow.minyr, &environs_rng)
                        : 0;
         break;
     }
     Plot.disturbance = event;

   }

}
