/**
 * \file ST_stats.c
 * \brief Record keeping for the entire model.
 * 
 * Statistics are kept for all metrics of a plant's lifecycle.
 * These statistics are kept for [resource groups](\ref RGROUP),
 * [species](\ref SPECIES), [individuals](\ref INDIVIDUAL), and 
 * [mortality events](\ref MORTALITY).
 * 
 *  History:
 *    (6/15/2000) -- INITIAL CODING - cwb
 *    1/9/01 - revised to make extensive use of malloc()
 *	  5/28/2013 (DLM) - added module level variable accumulators (grid_Stat) for the grid and functions.
 *    07/30/2016 (AKT) Fixed bug at std_dev calculation
 *    8/30/2019 (Chandler Haukap) - Removed grid_Stat accumulators and 
 *           functionality. Gridded mode accumulators are now stored in the
 *           Celltype struct in ST_grid.h and swapped into this module using 
 *           a call to stat_Copy_Accumulators().
 * 
 * \author Chris Bennett (initial coding)
 * \date 15 June 2000
 * 
 * \author DLM (added gridded mode variables and functions)
 * \date 28 May 2016
 * 
 * \author AKT (fixed bug in std_dev calculation)
 * \date 30 July 2016
 * 
 * \author Chandler Haukap (overhauled all changes made by DLM)
 * \author Kyle Palmquist
 * \author Freddy Pierson
 * \date 23 August 2019
 * 
 * \ingroup STATISTICS
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <string.h>
#include "ST_steppe.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "ST_stats.h" // Contains most of the function declarations.
#include "ST_seedDispersal.h"
#include "ST_globals.h"

/* ----------------- Local Variables --------------------- */
StatType *_Dist, *_Ppt, *_Temp,
  *_Grp, *_Gsize, *_Gpr, *_Gmort, *_Gestab,
  *_Spp, *_Indv, *_Smort, *_Sestab, *_Sreceived;

FireStatsType *_Gwf;

/*************** Local Function Declarations ***************/
static void _init( void);
static RealF _get_avg( struct accumulators_st *p);
static RealF _get_std( struct accumulators_st *p);

/** \brief A macro for collecting statistics.
 * 
 * \param p is a pointer to the \ref accumulators_st which is collecting the data.
 * \param v is a double which will be collected.
 * 
 * Note that the syntax checker is obviated, so make sure
 * you follow the this prototype:
 * static void _collect_add(struct accumulators_st *p, double v).
 * 
 * \ingroup STATISTICS_PRIVATE
 */
#define _collect_add(p, v) {					\
  (p)->nobs++;							\
  RealF old_ave = (p)->ave;					\
  (p)->ave = get_running_mean((p)->nobs, (p)->ave, (v));	\
  (p)->sum_dif_sqr += get_running_sqr(old_ave, (p)->ave, (v));	\
  (p)->sd = final_running_sd((p)->nobs, (p)->sum_dif_sqr);	\
}

/** \brief A macro that copies the data of p into v.
 * 
 * \param p is a pointer to the \ref accumulators_st to copy from.
 * \param v is a pointer to the \ref accumulators_st to copy to.
 * 
 * The correct usage is _copy_over(struct accumulators_st *p, struct accumulators_st *v).
 * 
 * \ingroup STATISTICS_PRIVATE
 */
#define _copy_over(p, v) { \
	(p)->ave = (v)->ave; \
	(p)->sum_dif_sqr = (v)->sum_dif_sqr; \
	(p)->sd = (v)->sd;   \
	(p)->nobs = (v)->nobs; \
}

static Bool firsttime = TRUE;


/***********************************************************/
/*              BEGIN FUNCTIONS                            */
/***********************************************************/

/**
 * \brief Collects the requested statistics.
 * 
 * \param year is the year you are collecting. year should be base 1.
 * 
 * Statistics are collected for every metric specified in bmassflags.in
 * and mortflags.in.
 * 
 * \sa Stat_Output() which is where the collected statistics are printed
 *     to CSV files.
 * 
 * \ingroup STATISTICS
 */
void stat_Collect( Int year ) {
/* fill data structures with samples to be
   computed later in Stat_Output().

   enter with year base1, but subtract 1 for index (base0)
*/

  SppIndex sp;
  GrpIndex rg;
  double bmass;

  if (firsttime) {
    firsttime = FALSE;
    _init();
  }

  year--;
  if (BmassFlags.dist && Plot->disturbed)
    _Dist->s[year].nobs++;

  if (BmassFlags.ppt)
    _collect_add( &_Ppt->s[year], Env->ppt);

  if (BmassFlags.tmp)
    _collect_add( &_Temp->s[year], Env->temp);

  if (BmassFlags.grpb) {
    if (BmassFlags.wildfire) {
        _Gwf->wildfire[year] += (RGroup[0]->wildfire) ? 1 : 0;
    }
    ForEachGroup(rg) {
      bmass = (double) RGroup_GetBiomass(rg);
      if ( LT(bmass, 0.0) ) {
        LogError(logfp, LOGWARN, "Grp %s biomass(%.4f) < 0 in stat_Collect()",
                        RGroup[rg]->name, bmass);
        bmass = 0.0;
      }
      _collect_add( &_Grp[rg].s[year], bmass);

      if (BmassFlags.size)
        _collect_add( &_Gsize[rg].s[year],
                          getRGroupRelsize(rg));
      if (BmassFlags.pr)
        _collect_add( &_Gpr[rg].s[year],
                          RGroup[rg]->pr);
      if (BmassFlags.prescribedfire){
        _Gwf->prescribedFire[rg][year] += (RGroup[rg]->prescribedfire) ? 1 : 0;
      }
    }
  }

  if (BmassFlags.sppb) {
    ForEachSpecies(sp) {
      bmass = (double) Species_GetBiomass(sp);
      if ( LT(bmass, 0.0) ) {
        LogError(logfp, LOGWARN, "Spp %s biomass(%.4f) < 0 in stat_Collect()",
                       Species[sp]->name, bmass);
        bmass = 0.0;
      }
      _collect_add( &_Spp[sp].s[year], bmass);

      if (BmassFlags.indv)
        _collect_add( &_Indv[sp].s[year],
            (double) Species[sp]->est_count);
    }
  }

  if(UseSeedDispersal && UseGrid) {
  	ForEachSpecies(sp)
		_collect_add( &_Sreceived[sp].s[year], (double) Species[sp]->received_prob);
  }
}


/**
 * \brief initialize the statistics accumulators.
 * 
 * The function will allocate the accumulators. Note that this function only has to be
 * called once per simulation.
 * 
 * \sideeffect Accumulators are allocated memory based on which statistics are requested
 *             in bmassflags.in and mortflags.in.
 * 
 * \ingroup STATISTICS_PRIVATE
 */
static void _init( void) {
/* must be called after model is initialized */
  SppIndex sp;
  GrpIndex rg;

  if (BmassFlags.dist) {
    _Dist = (StatType*) Mem_Calloc(1, sizeof(StatType), "_stat_init(Dist)");
    _Dist->s = (struct accumulators_st *)
               Mem_Calloc( SuperGlobals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Dist)");
  }
  if (BmassFlags.ppt) {
    _Ppt = (StatType*) Mem_Calloc(1, sizeof(StatType), "_stat_init(PPT");
    _Ppt->s  = (struct accumulators_st *)
               Mem_Calloc( SuperGlobals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(PPT)");
  }
  if (BmassFlags.tmp) {
    _Temp = (StatType*) Mem_Calloc(1, sizeof(StatType), "_stat_init(Temp)");
    _Temp->s = (struct accumulators_st *)
               Mem_Calloc( SuperGlobals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Temp)");
  }
  if (BmassFlags.grpb) {
    _Grp = (struct stat_st *)
           Mem_Calloc( Globals->grpCount,
                       sizeof(struct stat_st),
                      "_stat_init(Grp)");
    ForEachGroup(rg)
      _Grp[rg].s = (struct accumulators_st *)
             Mem_Calloc( SuperGlobals.runModelYears,
                         sizeof(struct accumulators_st),
                        "_stat_init(Grp[rg].s)");

    if (BmassFlags.size) {
      _Gsize = (struct stat_st *)
             Mem_Calloc( Globals->grpCount,
                         sizeof(struct stat_st),
                        "_stat_init(GSize)");
      ForEachGroup(rg)
          _Gsize[rg].s = (struct accumulators_st *)
             Mem_Calloc( SuperGlobals.runModelYears,
                         sizeof(struct accumulators_st),
                        "_stat_init(GSize[rg].s)");
    }
    if (BmassFlags.pr) {
      _Gpr = (struct stat_st *)
             Mem_Calloc( Globals->grpCount,
                         sizeof(struct stat_st),
                        "_stat_init(Gpr)");
      ForEachGroup(rg)
          _Gpr[rg].s = (struct accumulators_st *)
             Mem_Calloc( SuperGlobals.runModelYears,
                         sizeof(struct accumulators_st),
                        "_stat_init(Gpr[rg].s)");
    }

    if (BmassFlags.wildfire || BmassFlags.prescribedfire) {
      _Gwf = (struct fire_st *)
             Mem_Calloc( 1,
                         sizeof(struct fire_st),
                        "_stat_init(Gwf)");

      _Gwf->wildfire = (int *)
          Mem_Calloc( 1,
                      sizeof(int) * SuperGlobals.runModelYears,
                      "_stat_init(Gwf->wildfire)");
      
      _Gwf->prescribedFire = (int **)
          Mem_Calloc( 1,
                      sizeof(int **) * SuperGlobals.max_rgroups,
                      "_stat_init(Gwf->prescribedfire");

      ForEachGroup(rg){
        _Gwf->prescribedFire[rg] = (int *)
          Mem_Calloc( SuperGlobals.runModelYears,
                      sizeof(int) * SuperGlobals.runModelYears,
                      "_stat_init(Gwf->prescribedFire)");
      }
    }
  }

  if (MortFlags.group) {

    _Gestab = (struct stat_st *)
             Mem_Calloc( Globals->grpCount,
                         sizeof(struct stat_st),
                         "_stat_init(Gestab)");
    ForEachGroup(rg)
      _Gestab[rg].s = (struct accumulators_st *)
                     Mem_Calloc( 1, sizeof(struct accumulators_st),
                                "_stat_init(Gestab[rg].s)");

    _Gmort = (struct stat_st *)
           Mem_Calloc( Globals->grpCount,
                       sizeof(struct stat_st),
                      "_stat_init(Gmort)");
    ForEachGroup(rg)
        _Gmort[rg].s = (struct accumulators_st *)
           Mem_Calloc( GrpMaxAge(rg),
                       sizeof(struct accumulators_st),
                      "_stat_init(Gmort[rg].s)");
  }

  if (BmassFlags.sppb) {
      _Spp = (struct stat_st *)
               Mem_Calloc( Globals->sppCount,
                           sizeof(struct stat_st),
                          "_stat_init(Spp)");
      ForEachSpecies(sp)
        _Spp[sp].s = (struct accumulators_st *)
               Mem_Calloc( SuperGlobals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Spp[sp].s)");

      if (BmassFlags.indv) {
        _Indv = (struct stat_st *)
               Mem_Calloc( Globals->sppCount,
                           sizeof(struct stat_st),
                          "_stat_init(Indv)");
        ForEachSpecies(sp)
          _Indv[sp].s = (struct accumulators_st *)
               Mem_Calloc( SuperGlobals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Indv[sp].s)");
    }
  }
  if (MortFlags.species) {
    _Sestab = (struct stat_st *)
           Mem_Calloc( Globals->sppCount,
                       sizeof(struct stat_st),
                      "_stat_init(Sestab)");
    ForEachSpecies(sp)
      _Sestab[sp].s = (struct accumulators_st *)
                    Mem_Calloc( 1, sizeof(struct accumulators_st),
                                "_stat_init(Sestab[sp].s)");

    _Smort = (struct stat_st *)
           Mem_Calloc( Globals->sppCount,
                       sizeof(struct stat_st),
                      "_stat_init(Smort)");
    ForEachSpecies(sp)
      _Smort[sp].s = (struct accumulators_st *)
                    Mem_Calloc( SppMaxAge(sp),
                                sizeof(struct accumulators_st),
                                "_stat_init(Smort[sp].s)");
  }

  if (UseSeedDispersal && UseGrid) {
	  _Sreceived = Mem_Calloc( Globals->sppCount, sizeof(struct stat_st), "_stat_init(Sreceived)");
	  ForEachSpecies(sp) {
		  _Sreceived[sp].s = (struct accumulators_st *)Mem_Calloc( SuperGlobals.runModelYears, sizeof(struct accumulators_st), "_stat_init(Sreceived[sp].s)");
		  _Sreceived[sp].name = &Species[sp]->name[0];
	  }
  }

  /* "appoint" names of columns*/
  if (BmassFlags.grpb) {
    ForEachGroup(rg)
      _Grp[rg].name = &RGroup[rg]->name[0];
  }
  if (MortFlags.group) {
    ForEachGroup(rg)
      _Gmort[rg].name = &RGroup[rg]->name[0];
  }
  if (BmassFlags.sppb) {
    ForEachSpecies(sp)
      _Spp[sp].name = &Species[sp]->name[0];
  }
  if (MortFlags.species) {
    ForEachSpecies(sp)
      _Smort[sp].name = &Species[sp]->name[0];
  }
}

/* Shallow copies StatType and FireStatsType pointers to the local pointers.
   This is intended to be used with gridded mode to load in a given cell */
void stat_Copy_Accumulators(StatType* newDist, StatType* newPpt, StatType* newTemp, StatType* newGrp, StatType* newGsize, 
                            StatType* newGpr, StatType* newGmort, StatType* newGestab, StatType* newSpp, StatType* newIndv,
                            StatType* newSmort, StatType* newSestab, StatType* newSrecieved, FireStatsType* newGwf, Bool firstTime){

  /* Move the local pointers to the location of the given pointers */
  _Dist = newDist;
  _Ppt = newPpt;
  _Temp = newTemp;
  _Grp = newGrp;
  _Gsize = newGsize;
  _Gpr = newGpr;
  _Gmort = newGmort;
  _Gestab = newGestab;
  _Spp = newSpp;
  _Indv = newIndv;
  _Smort = newSmort;
  _Sestab = newSestab;
  _Sreceived = newSrecieved;
  _Gwf = newGwf;
  firsttime = firstTime;
}

/***********************************************************/
void stat_free_mem( void ) {
	//frees memory allocated in this module
	GrpIndex gp;
	SppIndex sp;

  	if(BmassFlags.grpb)
  		ForEachGroup(gp) {
  			Mem_Free(_Grp[gp].s);
  			if (BmassFlags.size) Mem_Free(_Gsize[gp].s);
  			if (BmassFlags.pr) Mem_Free(_Gpr[gp].s);
  		}
  	if(BmassFlags.sppb)
  		ForEachSpecies(sp) {
  			Mem_Free(_Spp[sp].s);
  			if(BmassFlags.indv) Mem_Free(_Indv[sp].s);
  		}

    if (BmassFlags.wildfire || BmassFlags.prescribedfire){
      Mem_Free(_Gwf->wildfire);
      ForEachGroup(gp) {
        Mem_Free(_Gwf->prescribedFire[gp]);
  		}
      Mem_Free(_Gwf->prescribedFire);
    }

  	if (BmassFlags.dist) {
      Mem_Free(_Dist->s);
      Mem_Free(_Dist);
    }
  	if (BmassFlags.ppt) {
      Mem_Free(_Ppt->s);
      Mem_Free(_Ppt);
    }
  	if (BmassFlags.tmp) {
      Mem_Free(_Temp->s);
      Mem_Free(_Temp);
    }

  	if(BmassFlags.grpb) {
  		Mem_Free(_Grp);
  		if (BmassFlags.size) Mem_Free(_Gsize);
  		if (BmassFlags.size) Mem_Free(_Gpr);
  	}
  	if (MortFlags.group) {
  		ForEachGroup(gp) {
  			Mem_Free(_Gmort[gp].s);
  			Mem_Free(_Gestab[gp].s);
  		}
  		Mem_Free(_Gmort);
  		Mem_Free(_Gestab);
  	}
  	if(BmassFlags.sppb) {
  		Mem_Free(_Spp);
  		if(BmassFlags.indv) Mem_Free(_Indv);
  	}
  	if (MortFlags.species) {
  		ForEachSpecies(sp) {
  			Mem_Free(_Smort[sp].s);
  			Mem_Free(_Sestab[sp].s);
  		}
  		Mem_Free(_Smort);
  		Mem_Free(_Sestab);
  	}


	if (UseSeedDispersal && UseGrid) {
		ForEachSpecies(sp)
			Mem_Free(_Sreceived[sp].s);
		Mem_Free(_Sreceived);
	}


}

/**
 * \brief Collects mortality statistics across iterations for all entries in \ref RGroup.
 * 
 * Mortality statistics are accumulated in species_Update_Kills(). stat_Collect_GMort should 
 * be called after every iteration to add the iteration to the simulation statistics.
 * 
 * \sideeffect \ref _Gmort will be modified according to the last iteration's mortality stats.
 * 
 * \sa species_Update_Kills().
 * 
 * \ingroup STATISTICS
 */
void stat_Collect_GMort ( void ) {
    IntS rg, age;

    ForEachGroup(rg) {
      if (!RGroup[rg]->use_me) continue;
      _collect_add( _Gestab[rg].s, RGroup[rg]->estabs);
      for (age=0; age < GrpMaxAge(rg); age++)
        _collect_add( &_Gmort[rg].s[age],
              (double) RGroup[rg]->kills[age]);

    }

}

/**
 * \brief Collects mortality statistics across iterations for all entries in \ref Species.
 * 
 * Mortality statistics are accumulated in species_Update_Kills(). stat_Collect_SMort should 
 * be called after every iteration to add the iteration to the simulation statistics.
 * 
 * \sideeffect \ref _Smort will be modified according to the last iteration's mortality stats.
 * 
 * \sa species_Update_Kills().
 * 
 * \ingroup STATISTICS
 */
void stat_Collect_SMort ( void ) {
   SppIndex sp;
   IntS age;

  ForEachSpecies(sp) {
    if ( !Species[sp]->use_me) continue;
    _collect_add( _Sestab[sp].s, Species[sp]->estabs);
    for (age=0; age < SppMaxAge(sp); age++)
        _collect_add( &_Smort[sp].s[age],
              (double) Species[sp]->kills[age]);

  }

}

/**
 * \brief Prints mortality statistics to the file specified in Globals.mort.fp_year.
 * 
 * This function Will create the header and all entries in the yearly mortality output
 * file. The statistics output are those specified in mortflags.in.
 * 
 * \ingroup STATISTICS
 */
void stat_Output_YrMorts( void ) {

  FILE *f = Globals->mort.fp_year;
  IntS age;
  GrpIndex rg;
  SppIndex sp;
  char sep = MortFlags.sep;

  if (!MortFlags.yearly) return;

  fprintf(f,"Age");
  if (MortFlags.group) {
    ForEachGroup(rg) fprintf(f,"%c%s", sep, RGroup[rg]->name);
  }
  if (MortFlags.species) {
    ForEachSpecies(sp)
      fprintf(f,"%c%s", sep, Species[sp]->name);
  }
  fprintf(f,"\n");
  /* end of the first line */

  fprintf(f,"Estabs");
  if (MortFlags.group) {
    ForEachGroup(rg)
      fprintf(f,"%c%d", sep, RGroup[rg]->estabs);
  }
  if (MortFlags.species) {
    ForEachSpecies(sp)
      fprintf(f,"%c%d", sep, Species[sp]->estabs);
  }
  fprintf(f,"\n");

  /* print one line of kill frequencies per age */
  for(age=0; age < Globals->Max_Age; age++) {
    fprintf(f,"%d", age+1);
    if (MortFlags.group) {
      ForEachGroup(rg){
        if( age < GrpMaxAge(rg) )
          fprintf(f,"%c%d", sep, RGroup[rg]->kills[age]);
        else
          fprintf(f,"%c", sep);
      }
    }
    if (MortFlags.species) {
      ForEachSpecies(sp) {
        if (age < SppMaxAge(sp))
          fprintf(f,"%c%d", sep, Species[sp]->kills[age]);
        else
          fprintf(f,"%c", sep);
      }
    }
    fprintf(f,"\n");
  }

  CloseFile(&f);
}

/**
 * \brief Outputs all mortality statistics.
 * 
 * The file they are printed to is denoted by \ref Parm_name().
 * The statistics printed are those denoted in the mortflags.in file.
 * 
 * \ingroup STATISTICS
 */
void stat_Output_AllMorts( void) {
  FILE *f;
  IntS age;
  GrpIndex rg;
  SppIndex sp;
  char sep = MortFlags.sep;

  if (!MortFlags.summary) return;

  f = OpenFile( Parm_name(F_MortAvg), "w");

  fprintf(f,"Age");
  if (MortFlags.group) {
    ForEachGroup(rg) fprintf(f,"%c%s", sep, RGroup[rg]->name);
  }
  if (MortFlags.species) {
    ForEachSpecies(sp)
      fprintf(f,"%c%s", sep, Species[sp]->name);
  }
  fprintf(f,"\n");
  /* end of first line */

  /* print one line of establishments */
  fprintf(f,"Estabs");
  if (MortFlags.group) {
    ForEachGroup(rg)
      fprintf(f,"%c%5.1f", sep, _get_avg(_Gestab[rg].s));
  }
  if (MortFlags.species) {
    ForEachSpecies(sp)
      fprintf(f,"%c%5.1f", sep, _get_avg( _Sestab[sp].s));
  }
  fprintf(f,"\n");

  /* print one line of kill frequencies per age */
  for(age=0; age < Globals->Max_Age; age++) {
  fprintf(f,"%d", age+1);
  if (MortFlags.group) {
      ForEachGroup(rg)
        fprintf(f,"%c%5.1f", sep, ( age < GrpMaxAge(rg) )
                                  ? _get_avg(&_Gmort[rg].s[age])
                                  : 0.);
    }
    if (MortFlags.species) {
      ForEachSpecies(sp) {
      fprintf(f,"%c%5.1f", sep, ( age < SppMaxAge(sp))
                                ? _get_avg(&_Smort[sp].s[age])
                                : 0.);
    }
    }
  fprintf(f,"\n");
  }

  CloseFile(&f);
}

/***********************************************************/
void stat_Output_AllBmass(void) {

  char buf[2048], tbuf[80], sep = BmassFlags.sep;
  IntS yr;
  GrpIndex rg;
  SppIndex sp;
  FILE *f;

  if (!BmassFlags.summary) return;

  f = OpenFile( Parm_name( F_BMassAvg), "w");

  buf[0]='\0';

  if (BmassFlags.header) {
	make_header_with_std(buf);
    fprintf(f, "%s", buf);
  }

  for( yr=1; yr<= SuperGlobals.runModelYears; yr++) {
    *buf = '\0';
    if (BmassFlags.yr)
      sprintf(buf, "%d%c", yr, sep);

    if (BmassFlags.dist) {
      sprintf(tbuf, "%ld%c", _Dist->s[yr-1].nobs,
              sep);
      strcat(buf, tbuf);
    }

    if (BmassFlags.ppt) {
      sprintf(tbuf, "%f%c%f%c",
              _get_avg(&_Ppt->s[yr-1]), sep,
              _get_std(&_Ppt->s[yr-1]), sep);
      strcat( buf, tbuf);
    }

    if (BmassFlags.pclass) {
      sprintf(tbuf, "\"NA\"%c", sep);
      strcat( buf, tbuf);
    }

    if (BmassFlags.tmp) {
      sprintf(tbuf, "%f%c%f%c",
              _get_avg(&_Temp->s[yr-1]), sep,
              _get_std(&_Temp->s[yr-1]), sep);
      strcat( buf, tbuf);
    }

    if (BmassFlags.grpb) {
      if (BmassFlags.wildfire) {
          sprintf(tbuf, "%d%c",
                  _Gwf->wildfire[yr-1], sep);
          strcat( buf, tbuf);
      }
      ForEachGroup(rg) {
        sprintf(tbuf, "%f%c%f%c",
                _get_avg(&_Grp[rg].s[yr-1]), sep,
                _get_std(&_Grp[rg].s[yr-1]), sep);
        strcat( buf, tbuf);

        if (BmassFlags.size) {
          sprintf(tbuf, "%f%c",
                  _get_avg( &_Gsize[rg].s[yr-1]), sep);
          strcat( buf, tbuf);
        }

        if (BmassFlags.pr) {
          sprintf(tbuf, "%f%c%f%c",
                  _get_avg( &_Gpr[rg].s[yr-1]), sep,
                  _get_std( &_Gpr[rg].s[yr-1]), sep);
          strcat( buf, tbuf);
        }
        /* Output the sum of all the wild and prescribed fire numbers across all iterations;
         ATTENTION: the other output index are the average values across all iterations*/
        if (BmassFlags.prescribedfire) {
          sprintf(tbuf, "%d%c",
                  _Gwf->prescribedFire[rg][yr-1], sep);
          strcat( buf, tbuf);
        }
      }
    }

		if (BmassFlags.sppb)
		{
			for ((sp) = 0; (sp) < Globals->sppCount - 1; (sp)++)
			{
				sprintf(tbuf, "%f%c", _get_avg(&_Spp[sp].s[yr - 1]), sep);
				strcat(buf, tbuf);

				if (BmassFlags.indv)
				{
					sprintf(tbuf, "%f%c", _get_avg(&_Indv[sp].s[yr - 1]), sep);
					strcat(buf, tbuf);
				}
			}

			if (BmassFlags.indv)
			{
				sprintf(tbuf, "%f%c", _get_avg(&_Spp[sp].s[yr - 1]), sep);
				strcat(buf, tbuf);

				sprintf(tbuf, "%f", _get_avg(&_Indv[sp].s[yr - 1]));
				strcat(buf, tbuf);
			}
			else
			{
				sprintf(tbuf, "%f", _get_avg(&_Spp[sp].s[yr - 1]));
				strcat(buf, tbuf);
			}

		}

    fprintf( f, "%s\n", buf);
  }  /* end of foreach year */
  CloseFile(&f);

}

/**
 * \brief returns the average value of an accumulator.
 * 
 * This function works, but is deprecated. You can reference the
 * average directly with p->ave.
 * 
 * \param p is a pointer to the \ref accumulators_st.
 * 
 * \ingroup STATISTICS_PRIVATE
 */
static RealF _get_avg( struct accumulators_st *p) 
{
	return p->ave;
}

/**
 * \brief returns the standard deviation of an accumulator.
 * 
 * This function works, but is deprecated. You can reference the
 * standard deviation directly with p->sd.
 * 
 * \param p is a pointer to the \ref accumulators_st.
 * 
 * \ingroup STATISTICS_PRIVATE
 */
static RealF _get_std(struct accumulators_st *p)
{
	return p->sd;
}

/**
 * \brief prints the header for biomass statistics with standard deviations.
 * 
 * This function is called when a header is requested in bmassflags.in.
 * stat_Output_AllBmass() takes care of calling this function when
 * requested.
 * 
 * \sa _make_header()
 * 
 * \ingroup STATISTICS_PRIVATE
 */
void make_header_with_std( char *buf) {

  char **fields;
  char tbuf[80];
  GrpIndex rg;
  SppIndex sp;
  Int i, fc=0;
  
  fields = (char **)Mem_Calloc(MAX_OUTFIELDS * 2, sizeof(char *), "make_header_with_std");
  
  for (i = 0; i < MAX_OUTFIELDS * 2; i++) {
      fields[i] = (char *)Mem_Calloc(MAX_FIELDLEN + 1, sizeof(char), "make_header_with_std");
  }

  /* Set up headers */
  if (BmassFlags.yr)
    strcpy(fields[fc++], "Year");
  if (BmassFlags.dist)
    strcpy(fields[fc++], "Disturbs");
  if (BmassFlags.ppt) {
    strcpy(fields[fc++], "PPT");
    strcpy(fields[fc++], "StdDev");
  }
  if (BmassFlags.pclass)
    strcpy(fields[fc++], "PPTClass");
  if (BmassFlags.tmp) {
    strcpy(fields[fc++], "Temp");
    strcpy(fields[fc++], "StdDev");
  }
  if (BmassFlags.grpb) {
    if (BmassFlags.wildfire) {
        strcpy(fields[fc++], "WildFire");
    }

    ForEachGroup(rg) {
      strcpy(fields[fc++], RGroup[rg]->name);

      strcpy(fields[fc], RGroup[rg]->name);
      strcat(fields[fc++], "_std");

      if (BmassFlags.size) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "_RSize");
      }
      if (BmassFlags.pr) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++],"_PR");
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "_PRstd");
      }
      if (BmassFlags.prescribedfire) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "_PFire");
      }
    }
  }

  if (BmassFlags.sppb) {
    ForEachSpecies(sp) {
      strcpy(fields[fc++], Species[sp]->name);
      if (BmassFlags.indv) {
        strcpy(fields[fc], Species[sp]->name);
        strcat(fields[fc++], "_Indivs");
      }
    }
  }



  /* Put header line in global variable */
    for (i=0; i< fc-1; i++) {
      sprintf(tbuf,"%s%c", fields[i], BmassFlags.sep);
      strcat(buf, tbuf);
    }
    sprintf(tbuf,"%s\n", fields[i]);
    strcat(buf, tbuf);

    for (i = 0; i < MAX_OUTFIELDS * 2; i++) {
        Mem_Free(fields[i]);
    }
    
    Mem_Free(fields);
}

/**
 * \brief prints the header for biomass statistics without standard deviations.
 * 
 * This function is called when a header is requested in bmassflags.in.
 * stat_Output_AllBmass() takes care of calling this function when
 * requested.
 * 
 * \sa _make_header_with_std()
 * 
 * \ingroup STATISTICS_PRIVATE
 */
void make_header( char *buf) {

  char **fields;
  char tbuf[80];
  GrpIndex rg;
  SppIndex sp;
  Int i, fc=0;

  fields = (char **)Mem_Calloc(MAX_OUTFIELDS * 2, sizeof(char *), "make_header");
  
  for (i = 0; i < MAX_OUTFIELDS * 2; i++) {
      fields[i] = (char *)Mem_Calloc(MAX_FIELDLEN + 1, sizeof(char), "make_header");
  }
  
  /* Set up headers */
  if (BmassFlags.yr)
    strcpy(fields[fc++], "Year");
  if (BmassFlags.dist)
    strcpy(fields[fc++], "Disturbs");
  if (BmassFlags.ppt) {
    strcpy(fields[fc++], "PPT");
    strcpy(fields[fc++], "StdDev");
  }
  if (BmassFlags.pclass)
    strcpy(fields[fc++], "PPTClass");
  if (BmassFlags.tmp) {
    strcpy(fields[fc++], "Temp");
    strcpy(fields[fc++], "StdDev");
  }
  if (BmassFlags.grpb) {
    if (BmassFlags.wildfire) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "WildFire");
    }
    ForEachGroup(rg) {
      strcpy(fields[fc++], RGroup[rg]->name);
      if (BmassFlags.size) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "_RSize");
      }
      if (BmassFlags.pr) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++],"_PR");
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "_PRstd");
      }
      if (BmassFlags.prescribedfire) {
        strcpy(fields[fc], RGroup[rg]->name);
        strcat(fields[fc++], "_PrescribedFire");
      }
    }
  }

  if (BmassFlags.sppb) {
    ForEachSpecies(sp) {
      strcpy(fields[fc++], Species[sp]->name);
      if (BmassFlags.indv) {
        strcpy(fields[fc], Species[sp]->name);
        strcat(fields[fc++], "_Indivs");
      }
    }
  }

  /* Put header line in global variable */
    for (i=0; i< fc-1; i++) {
      sprintf(tbuf,"%s%c", fields[i], BmassFlags.sep);
      strcat(buf, tbuf);
    }
    sprintf(tbuf,"%s\n", fields[i]);
    strcat(buf, tbuf);

    for (i = 0; i < MAX_OUTFIELDS * 2; i++) {
        Mem_Free(fields[i]);
    }
    
    Mem_Free(fields);
}

#ifdef DEBUG_MEM
void Stat_SetMemoryRefs(void) {
/* when debugging memory problems, use the bookkeeping
   code in myMemory.c
 This routine sets the known memory refs in this module
 so they can be  checked for leaks, etc.  All refs will
 have been cleared by a call to ClearMemoryRefs() before
 this, and will be checked via CheckMemoryRefs() after
 this, most likely in the main() function.

 EVERY dynamic allocation must be noted here or the
 check will fail (which is the point, to catch unknown
 or missing pointers to memory).

*/

// DLM - 6/6/2013 : NOTE - The dynamically allocated grid_Stat variable is not accounted for here at the moment, it might need adding for this function to work correctly.  I've been using the valgrind program to debug memory errors, so I haven't bothered adding them since I haven't been using it.

  SppIndex sp;
  GrpIndex rg;

  if (BmassFlags.dist)
    NoteMemoryRef(_Dist.s);

  if (BmassFlags.ppt)
    NoteMemoryRef(_Ppt.s);

  if (BmassFlags.tmp)
    NoteMemoryRef(_Temp.s);

  if (BmassFlags.grpb) {
    NoteMemoryRef(_Grp);

    ForEachGroup(rg)
      NoteMemoryRef(_Grp[rg].s);

    if (BmassFlags.size) {
      NoteMemoryRef(_Gsize);

      ForEachGroup(rg)
        NoteMemoryRef(_Gsize[rg].s);
    }
  }

  if (MortFlags.group) {
    NoteMemoryRef(_Gestab);
    ForEachGroup(rg)
      NoteMemoryRef(_Gestab[rg].s);

    NoteMemoryRef(_Gmort);

    ForEachGroup(rg)
      NoteMemoryRef(_Gmort[rg].s);
  }

  if (BmassFlags.sppb) {
    NoteMemoryRef(_Spp);
    ForEachSpecies(sp)
      NoteMemoryRef(_Spp[sp].s);

    if (BmassFlags.indv) {
      NoteMemoryRef(_Indv);
      ForEachSpecies(sp)
        NoteMemoryRef(_Indv[sp].s);
    }
  }
  if (MortFlags.species) {
    NoteMemoryRef(_Sestab);
    ForEachSpecies(sp)
      NoteMemoryRef(_Sestab[sp].s);

    NoteMemoryRef(_Smort);
    ForEachSpecies(sp)
      NoteMemoryRef(_Smort[sp].s);
  }
}

#endif
