/** 
 * \file ST_environs.c
 * \brief Controls all environmental phenomenon from creating ppt and temp to 
 *        creating the disturbances.
 * 
 * \author CWB (initial coding)
 * \date 15 June 2000
 * \ingroup ENVIRONMENT
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/pcg/pcg_basic.h"
#include "sw_src/rands.h"
#include "sxw_funcs.h"
#include "sw_src/filefuncs.h"
extern SXW_t* SXW;

/*********** Locally Used Function Declarations ************/
/***********************************************************/

static void _make_ppt( void);
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

  SXW_Run_SOILWAT();

  _make_ppt();
  _set_ppt_reduction();
  _set_temp_reduction();
  _make_disturbance();
}

/* Deep copy one EnvType's information to another. Note that both must be initialized
   BEFORE calling this function. */
void copy_environment(const EnvType* src, EnvType* dest){
  if(!src){
    return;
  }

  dest->lyppt = src->lyppt;
  dest->ppt = src->ppt;
  dest->temp = src->temp;
  dest->gsppt = src->gsppt;
  dest->temp_reduction[0] = src->temp_reduction[0];
  dest->temp_reduction[1] = src->temp_reduction[1];
  dest->wet_dry = src->wet_dry;
}

/* Deep copy one PlotType's information to another. Note that both must be initialized
   BEFORE calling this function. */
void copy_plot(const PlotType* src, PlotType* dest){
  if(!src){
    return;
  }

  dest->disturbance = src->disturbance;
  dest->disturbed = src->disturbed;
  dest->pat_removed = src->pat_removed;
}

/* Deep copy one SucculentType's information to another. Note that both must be initialized
   BEFORE calling this function. */
void copy_succulent(const SucculentType* src, SucculentType* dest){
  if(!src){
    return;
  }

  dest->growth[0] = src->growth[0];
  dest->growth[1] = src->growth[1];
  dest->mort[0] = src->mort[0];
  dest->mort[1] = src->mort[1];
  dest->prob_death = src->prob_death;
  dest->reduction = src->reduction;
}

/**************************************************************/
static void _make_ppt( void) {
/*======================================================*/
/* If not running SOILWAT,take a random number from normal distribution with*/
/* mean, stddev that is between min & max from */
/* the Globals->ppt structure.*/
/* Also determine growing season precipitation. */

/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */
/* cwb - 6-Dec-02 -- added code to interface with STEPWAT.
 *       The ppt and gsppt are set in _sxw_set_environs()
 *       but we still pass through this code to set the
 *       Dry/Wet/Normal state. */
/*------------------------------------------------------*/

  IntS i;

#ifdef DEBUG_ENVCONST
  IntS r=320;
#endif

  // Run with SOILWAT2: we have monthly PPT and temperature to calculate
  // growing season precipitation as sum of monthly precipitation of those
  // months when mean air temperature exceeds a threshold that allows for plant growth 'Globals.temp.gstemp'
  Env->gsppt = 0; // gsppt is defined as IntS and units are millimeters

  for (i = 0; i < MAX_MONTHS; i++)
  {
    Env->gsppt += GE(SXW->temp_monthly[i], Globals->temp.gstemp) ?
      (IntS) (SXW->ppt_monthly[i] * 10. + 0.5) : 0;
  }

  if (Env->gsppt <= 0)
  {
    LogError(logfp, LOGWARN, "Zero growing season precipitation in "\
      "year = %d of iteration = %d", Globals->currYear, Globals->currIter);
    Env->gsppt = 0;
  }

  if ( Env->ppt <= Globals->ppt.dry )
    Env->wet_dry = Ppt_Dry;
  else if (Env->ppt >= Globals->ppt.wet)
    Env->wet_dry = Ppt_Wet;
  else
    Env->wet_dry = Ppt_Norm;
}

/**************************************************************/
static void _set_ppt_reduction( void) {
/*======================================================*/


/* HISTORY */
/* Chris Bennett @ LTER-CSU 6/15/2000            */

/*------------------------------------------------------*/

/* EQN 10*/
  Succulent->reduction = fabs(Succulent->growth[Slope]
                           * Env->gsppt
                           + Succulent->growth[Intcpt]);
/* EQN 16*/
  Succulent->prob_death = (Succulent->mort[Slope] * Env->gsppt
                          + Succulent->mort[Intcpt]) / 100.0;
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
    tp[1] = Globals->tempparm[i][0];
    tp[2] = Globals->tempparm[i][1];
    tp[3] = Globals->tempparm[i][2];
    tp[0] = Env->temp + tp[1];
    Env->temp_reduction[i] = tp[2]*tp[0] + tp[3] * (tp[0]*tp[0]);
    Env->temp_reduction[i] = max(0., Env->temp_reduction[i]);
  }

  if (Env->temp < 9.5 ) {
    Env->temp_reduction[CoolSeason] = .9;
    Env->temp_reduction[WarmSeason] = .6;
  } else {
    Env->temp_reduction[CoolSeason] = .6;
    Env->temp_reduction[WarmSeason] = .9;
  }
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
  if (Plot->disturbance != NoDisturb) {
    switch( Plot->disturbance) {
      case FecalPat:
           if (Plot->pat_removed) {
             Plot->disturbed = 0;
             Plot->pat_removed = FALSE;
             Plot->disturbance = NoDisturb;
           } else {
             pc = Globals->pat.recol[Slope] * Plot->disturbed
                  + Globals->pat.recol[Intcpt];
             if (RandUni(&environs_rng) <= pc) {
               Plot->pat_removed = TRUE;
               /* slight effects for one year*/
               Plot->disturbed = 1;
             } else {
               Plot->pat_removed = FALSE;
               Plot->disturbed++;
             }
           }
           break;
      case NoDisturb: // Does nothing but prevent a compiler warning
      case LastDisturb:
      default:
           Plot->disturbed = (Plot->disturbed) ? (Plot->disturbed - 1) : 0;
           break;
    }
    if (Plot->disturbed == 0)
      Plot->disturbance = NoDisturb;

  }

  /* if the disturbance was expired above, */
  /* we can generate a new one immediately */
  if (Plot->disturbance == NoDisturb) {

    /* pick some type of disturbance (other than none)*/
    event = (DisturbEvent) RandUniIntRange(1, LastDisturb -1, &environs_rng);

    /* make sure this is off unless needed  */
    Plot->pat_removed = FALSE;
    switch( event) {
      case FecalPat:
       if (!Globals->pat.use) {event=NoDisturb; break;}
         event = (RandUni(&environs_rng) <= Globals->pat.occur)
               ? event : NoDisturb;
         if (event == NoDisturb) break;
         Plot->pat_removed = (RandUni(&environs_rng) <= Globals->pat.removal);
         Plot->disturbed = 0;
         break;
      case AntMound:
       if (!Globals->mound.use) {event=NoDisturb; break;}
         event = (RandUni(&environs_rng) <= Globals->mound.occur)
               ? event :NoDisturb;
         if (event == NoDisturb) break;

         Plot->disturbed = RandUniIntRange(Globals->mound.minyr,
                                       Globals->mound.maxyr,
                                       &environs_rng);
         break;
      case Burrow:
       if (!Globals->burrow.use) {event=NoDisturb; break;}
         event = (RandUni(&environs_rng) <= Globals->burrow.occur)
               ? event :NoDisturb;
         if (event == NoDisturb) break;

         Plot->disturbed = (Globals->burrow.minyr > 0)
                        ? RandUniIntRange(1, Globals->burrow.minyr, &environs_rng)
                        : 0;
         break;
      case NoDisturb:
        break;
      case LastDisturb:
        break;
      default:
        break;
     }
     Plot->disturbance = event;

   }

}
