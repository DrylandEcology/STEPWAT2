/********************************************************/
/********************************************************/
//  Source file: stats.c
//  Type: module
//  Application: STEPPE - plant community dynamics simulator
//  Purpose: This is where all of the statistics are kept
//           as the model runs.
//  History:
//     (6/15/2000) -- INITIAL CODING - cwb
//   1/9/01 - revised to make extensive use of malloc() */
//	5/28/2013 (DLM) - added module level variable accumulators (grid_Stat) for the grid and functions to deal with them (stat_Load_Accumulators(), stat_Save_Accumulators() stat_Free_Accumulators(), and stat_Init_Accumulators()).  These functions are called from ST_grid.c and manage the output accumulators so that the gridded version can output correctly.  The accumulators are dynamically allocated, so be careful with them.
// 07/30/2016 (AKT) Fixed bug at std_dev calculation
//
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ST_steppe.h"
#include "filefuncs.h"
#include "myMemory.h"
#include "ST_structs.h"

/************ External Variable Declarations ***************/
/***********************************************************/
#include "ST_globals.h"

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/


/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
  void stat_Collect( Int year ) ;
  void stat_Collect_GMort ( void ) ;
  void stat_Collect_SMort ( void ) ;
  void stat_Output_YrMorts( void ) ;
  void stat_Output_AllMorts( void) ;
  void stat_Output_AllBmass(void) ;
  //Adding below two functions for creating grid cells avg values output file
  void stat_Output_AllBmassAvg(void) ;
  void stat_Output_AllCellAvgBmass(const char * filename);
  void stat_Output_Seed_Dispersal(const char * filename, const char sep, Bool makeHeader); 
  void stat_free_mem( void ) ;
  
  void stat_Load_Accumulators( int cell, int year ); //these accumulators were added to use in the gridded option... there overall purpose is to save/load data to allow steppe to output correctly when running multiple grid cells
  void stat_Save_Accumulators( int cell, int year );
  void stat_Free_Accumulators( void );
  void stat_Init_Accumulators( void );

/************************ Local Structure Defs *************/
/***********************************************************/
struct accumulators_st {
  double sum, sum_sq;
  unsigned long nobs;
};

struct stat_st {
  char *name; /* array of ptrs to names in RGroup & Species */
  struct accumulators_st *s;
} _Dist, _Ppt, _Temp,
  *_Grp, *_Gsize, *_Gpr, *_Gmort, *_Gestab,
  *_Spp, *_Indv, *_Smort, *_Sestab, *_Sreceived;

typedef struct  {
  struct accumulators_st *dist, *temp, *ppt, **grp1, **gsize, **gpr2, 
  							**gmort, **gestab, **spp, **indv, **smort, **sestab, **sreceived;
} accumulators_grid_st;
  
accumulators_grid_st *grid_Stat;


// Local Structure for holding sum values of all the grid cells

struct accumulators_grid_cell_st {
  double sum, sum_std;
  unsigned long nobs;
};

struct stat_grid_cell_st {
  char *name; /* array of ptrs to names in RGroup & Species */
  struct accumulators_grid_cell_st *s;    /* array of holding all the years values */
} _Dist_grid_cell, _Ppt_grid_cell, _Temp_grid_cell,
  *_Grp_grid_cell, *_Gsize_grid_cell, *_Gpr_grid_cell, *_Gmort_grid_cell, *_Gestab_grid_cell,
  *_Spp_grid_cell, *_Indv_grid_cell, *_Smort_grid_cell, *_Sestab_grid_cell, *_Sreceived_grid_cell;




/*************** Local Function Declarations ***************/
/***********************************************************/
static void _init( void);
static RealF _get_avg( struct accumulators_st *p);
static RealF _get_std( struct accumulators_st *p);
//Adding below three functions for copying grid cell values,calculating avg and SD,these values will be used in grid cells avg output file
static void copyStruct(RealF val,RealF std_val,struct accumulators_grid_cell_st *p );
static RealF _get_gridcell_avg( struct accumulators_grid_cell_st *p);
static RealF _get_gridcell_std( struct accumulators_grid_cell_st *p);
static void _make_header( char *buf);
static void _make_header_with_std( char *buf);

/* I'm making this a macro because it gets called a lot, but
/* note that the syntax checker is obviated, so make sure
/* you follow the this prototype:
/* static void _collect_add(struct accumulators_st *p, double v) */
#define _collect_add(p, v) { \
   (p)->sum += (v);              \
   (p)->sum_sq += (v)*(v);       \
   (p)->nobs++;                  \
}

// quick macro to make life easier in the load/save accumulators functions... it just copies the data of p into v
// static void _copy_over(struct accumulators_st *p, struct accumulators_st *v) 
#define _copy_over(p, v) { \
	(p)->sum = (v)->sum; \
	(p)->sum_sq = (v)->sum_sq; \
	(p)->nobs = (v)->nobs; \
}

static Bool firsttime = TRUE;


/***********************************************************/
/*              BEGIN FUNCTIONS                            */
/***********************************************************/

/***********************************************************/
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
  if (BmassFlags.dist && Plot.disturbed)
    _Dist.s[year].nobs++;

  if (BmassFlags.ppt)
    _collect_add( &_Ppt.s[year], Env.ppt);

  if (BmassFlags.tmp)
    _collect_add( &_Temp.s[year], Env.temp);

  if (BmassFlags.grpb) {
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
                          RGroup[rg]->relsize);
      if (BmassFlags.pr)
        _collect_add( &_Gpr[rg].s[year],
                          RGroup[rg]->pr);
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


/***********************************************************/
static void _init( void) {
/* must be called after model is initialized */
  SppIndex sp;
  GrpIndex rg;

// Memory allocation for local structures that will hold grid cells values for calculating avg
  if (UseGrid)
   {

	  if (BmassFlags.dist)
	      _Dist_grid_cell.s = (struct accumulators_grid_cell_st *)
	                 Mem_Calloc( Globals.runModelYears,
	                             sizeof(struct accumulators_grid_cell_st),
	                            "_stat_init(Dist)");
	    if (BmassFlags.ppt)
	      _Ppt_grid_cell.s = (struct accumulators_grid_cell_st *)
	                 Mem_Calloc( Globals.runModelYears,
	                             sizeof(struct accumulators_grid_cell_st),
	                            "_stat_init(PPT)");
	    if (BmassFlags.tmp)
	      _Temp_grid_cell.s = (struct accumulators_grid_cell_st *)
	                 Mem_Calloc( Globals.runModelYears,
	                             sizeof(struct accumulators_grid_cell_st),
	                            "_stat_init(Temp)");

		if (BmassFlags.grpb)
		{
			_Grp_grid_cell = (struct stat_grid_cell_st *)
					 Mem_Calloc( Globals.grpCount, sizeof(struct stat_grid_cell_st), "_stat_init(Grp)");
			ForEachGroup(rg)
				_Grp_grid_cell[rg].s = (struct accumulators_grid_cell_st *)
					 Mem_Calloc( Globals.runModelYears, sizeof(struct accumulators_grid_cell_st), "_stat_init(Grp[rg].s)");

			if (BmassFlags.size)
			{
				_Gsize_grid_cell = (struct stat_grid_cell_st *)
						Mem_Calloc( Globals.grpCount, sizeof(struct stat_grid_cell_st), "_stat_init(GSize)");
				ForEachGroup(rg)
					_Gsize_grid_cell[rg].s = (struct accumulators_grid_cell_st *)
					    Mem_Calloc( Globals.runModelYears, sizeof(struct accumulators_grid_cell_st), "_stat_init(GSize[rg].s)");
			}
			if (BmassFlags.pr)
			{
				_Gpr_grid_cell = (struct stat_grid_cell_st *)
						Mem_Calloc( Globals.grpCount, sizeof(struct stat_grid_cell_st), "_stat_init(Gpr)");
				ForEachGroup(rg)
					_Gpr_grid_cell[rg].s = (struct accumulators_grid_cell_st *)
					    Mem_Calloc( Globals.runModelYears, sizeof(struct accumulators_grid_cell_st), "_stat_init(Gpr[rg].s)");
			}
		}

		 if (BmassFlags.sppb)
		{
			_Spp_grid_cell = (struct stat_grid_cell_st *)
					Mem_Calloc( Globals.sppCount, sizeof(struct stat_grid_cell_st), "_stat_init(Spp)");
			ForEachSpecies(sp)
				_Spp_grid_cell[sp].s = (struct accumulators_grid_cell_st *)
				    Mem_Calloc( Globals.runModelYears, sizeof(struct accumulators_grid_cell_st), "_stat_init(Spp[sp].s)");

			if (BmassFlags.indv)
			{
				_Indv_grid_cell = (struct stat_grid_cell_st *)
						Mem_Calloc( Globals.sppCount, sizeof(struct stat_grid_cell_st), "_stat_init(Indv)");
				ForEachSpecies(sp)
					_Indv_grid_cell[sp].s = (struct accumulators_grid_cell_st *)
					    Mem_Calloc( Globals.runModelYears, sizeof(struct accumulators_grid_cell_st), "_stat_init(Indv[sp].s)");
			}
		}

  }

  if (BmassFlags.dist)
    _Dist.s = (struct accumulators_st *)
               Mem_Calloc( Globals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Dist)");
  if (BmassFlags.ppt)
    _Ppt.s  = (struct accumulators_st *)
               Mem_Calloc( Globals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(PPT)");
  if (BmassFlags.tmp)
    _Temp.s = (struct accumulators_st *)
               Mem_Calloc( Globals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Temp)");
  if (BmassFlags.grpb) {
    _Grp = (struct stat_st *)
           Mem_Calloc( Globals.grpCount,
                       sizeof(struct stat_st),
                      "_stat_init(Grp)");
    ForEachGroup(rg)
      _Grp[rg].s = (struct accumulators_st *)
             Mem_Calloc( Globals.runModelYears,
                         sizeof(struct accumulators_st),
                        "_stat_init(Grp[rg].s)");

    if (BmassFlags.size) {
      _Gsize = (struct stat_st *)
             Mem_Calloc( Globals.grpCount,
                         sizeof(struct stat_st),
                        "_stat_init(GSize)");
      ForEachGroup(rg)
          _Gsize[rg].s = (struct accumulators_st *)
             Mem_Calloc( Globals.runModelYears,
                         sizeof(struct accumulators_st),
                        "_stat_init(GSize[rg].s)");
    }
    if (BmassFlags.pr) {
      _Gpr = (struct stat_st *)
             Mem_Calloc( Globals.grpCount,
                         sizeof(struct stat_st),
                        "_stat_init(Gpr)");
      ForEachGroup(rg)
          _Gpr[rg].s = (struct accumulators_st *)
             Mem_Calloc( Globals.runModelYears,
                         sizeof(struct accumulators_st),
                        "_stat_init(Gpr[rg].s)");
    }
  }

  if (MortFlags.group) {

    _Gestab = (struct stat_st *)
             Mem_Calloc( Globals.grpCount,
                         sizeof(struct stat_st),
                         "_stat_init(Gestab)");
    ForEachGroup(rg)
      _Gestab[rg].s = (struct accumulators_st *)
                     Mem_Calloc( 1, sizeof(struct accumulators_st),
                                "_stat_init(Gestab[rg].s)");

    _Gmort = (struct stat_st *)
           Mem_Calloc( Globals.grpCount,
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
               Mem_Calloc( Globals.sppCount,
                           sizeof(struct stat_st),
                          "_stat_init(Spp)");
      ForEachSpecies(sp)
        _Spp[sp].s = (struct accumulators_st *)
               Mem_Calloc( Globals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Spp[sp].s)");

      if (BmassFlags.indv) {
        _Indv = (struct stat_st *)
               Mem_Calloc( Globals.sppCount,
                           sizeof(struct stat_st),
                          "_stat_init(Indv)");
        ForEachSpecies(sp)
          _Indv[sp].s = (struct accumulators_st *)
               Mem_Calloc( Globals.runModelYears,
                           sizeof(struct accumulators_st),
                          "_stat_init(Indv[sp].s)");
    }
  }
  if (MortFlags.species) {
    _Sestab = (struct stat_st *)
           Mem_Calloc( Globals.sppCount,
                       sizeof(struct stat_st),
                      "_stat_init(Sestab)");
    ForEachSpecies(sp)
      _Sestab[sp].s = (struct accumulators_st *)
                    Mem_Calloc( 1, sizeof(struct accumulators_st),
                                "_stat_init(Sestab[sp].s)");

    _Smort = (struct stat_st *)
           Mem_Calloc( Globals.sppCount,
                       sizeof(struct stat_st),
                      "_stat_init(Smort)");
    ForEachSpecies(sp)
      _Smort[sp].s = (struct accumulators_st *)
                    Mem_Calloc( SppMaxAge(sp),
                                sizeof(struct accumulators_st),
                                "_stat_init(Smort[sp].s)");
  }

  if (UseSeedDispersal && UseGrid) {
	  _Sreceived = Mem_Calloc( Globals.sppCount, sizeof(struct stat_st), "_stat_init(Sreceived)");
	  ForEachSpecies(sp) {
		  _Sreceived[sp].s = (struct accumulators_st *)Mem_Calloc( Globals.runModelYears, sizeof(struct accumulators_st), "_stat_init(Sreceived[sp].s)");
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

/***********************************************************/
void stat_Init_Accumulators( void ) {
	//allocates memory for all of the grid accumulators
	grid_Stat = Mem_Calloc(Globals.nCells, sizeof(accumulators_grid_st), "stat_Init_Accumulators()");
	
  	int i, j;
  	for( i = 0; i < Globals.nCells; i++) {
  		if (BmassFlags.dist) grid_Stat[i].dist = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  		if (BmassFlags.ppt) grid_Stat[i].ppt = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  		if (BmassFlags.tmp) grid_Stat[i].temp = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  		
  		if (BmassFlags.grpb) {
  			grid_Stat[i].grp1 = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st*), "stat_Init_Accumulators()"); // gave grp and gpr numbers attached to them so I wouldn't mix them up lol... bad (confusing) variable names on part of the original creator.
  			if (BmassFlags.size) grid_Stat[i].gsize = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  			if (BmassFlags.pr) grid_Stat[i].gpr2 = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
  		}
  		if (MortFlags.group) {
  			grid_Stat[i].gmort = Mem_Calloc(Globals.grpCount, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
  			grid_Stat[i].gestab = Mem_Calloc(Globals.grpCount, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
  			GrpIndex gp;
  			ForEachGroup(gp) {
  				grid_Stat[i].gmort[gp] = Mem_Calloc(GrpMaxAge(gp), sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  				grid_Stat[i].gestab[gp] = Mem_Calloc(1, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  			}
  		}
  		if(BmassFlags.sppb) {
  			grid_Stat[i].spp = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
  			if(BmassFlags.indv) grid_Stat[i].indv = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
  		}
  		if (MortFlags.species) {
  			grid_Stat[i].smort = Mem_Calloc(Globals.sppCount, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
  			grid_Stat[i].sestab = Mem_Calloc(Globals.sppCount, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
			SppIndex sp;
			ForEachSpecies(sp) {
				grid_Stat[i].smort[sp] = Mem_Calloc(SppMaxAge(sp), sizeof(struct accumulators_st), "stat_Init_Accumulators()");
				grid_Stat[i].sestab[sp] = Mem_Calloc(1, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
			}
		}

		if(UseSeedDispersal && UseGrid) 
			grid_Stat[i].sreceived = Mem_Calloc(Globals.runModelYears, sizeof(struct accumulators_st*), "stat_Init_Accumulators()");
		  		
  		
  		for( j = 0; j < Globals.runModelYears; j++) {
  			if (BmassFlags.grpb) {
  				grid_Stat[i].grp1[j] = Mem_Calloc(Globals.grpCount, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  				if (BmassFlags.size) grid_Stat[i].gsize[j] = Mem_Calloc(Globals.grpCount, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  				if (BmassFlags.pr) grid_Stat[i].gpr2[j] = Mem_Calloc(Globals.grpCount, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  			}
  			if(BmassFlags.sppb) {
  				grid_Stat[i].spp[j] = Mem_Calloc(Globals.sppCount, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  				if(BmassFlags.indv) grid_Stat[i].indv[j] = Mem_Calloc(Globals.sppCount, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  			}
			if(UseSeedDispersal && UseGrid)
				grid_Stat[i].sreceived[j] = Mem_Calloc(Globals.sppCount, sizeof(struct accumulators_st), "stat_Init_Accumulators()");
  		}
  	}
}

/***********************************************************/
void stat_Load_Accumulators(int cell, int year) {
	//loads the accumulators for the cell at the given year

	if (firsttime) {
		firsttime = FALSE;
		_init();
	}
	IntS age;
	int yr;
	yr = year - 1;
	
	if(MortFlags.species) {
		SppIndex sp;

		ForEachSpecies(sp) {
			if ( !Species[sp]->use_me) continue;
			_copy_over(&_Sestab[sp].s[0], &grid_Stat[cell].sestab[sp][0]);
			IntS age;
			for (age=0; age < SppMaxAge(sp); age++)
				_copy_over(&_Smort[sp].s[age], &grid_Stat[cell].smort[sp][age]);
		}
  	}
  	if(MortFlags.group) {
  		GrpIndex rg;

    	ForEachGroup(rg) {
      		if (!RGroup[rg]->use_me) continue;
      		_copy_over(&_Gestab[rg].s[0], &grid_Stat[cell].gestab[rg][0]);
      		for (age=0; age < GrpMaxAge(rg); age++)
        		_copy_over(&_Gmort[rg].s[age], &grid_Stat[cell].gmort[rg][age]);
    	}
  	}

	if (BmassFlags.tmp) _copy_over(&_Temp.s[yr], &grid_Stat[cell].temp[yr]);
	if (BmassFlags.ppt) _copy_over(&_Ppt.s[yr], &grid_Stat[cell].ppt[yr]);
	if (BmassFlags.dist) _copy_over(&_Dist.s[yr], &grid_Stat[cell].dist[yr]);	
	
	if(BmassFlags.grpb) {
		GrpIndex c;
		ForEachGroup(c) {
			_copy_over(&_Grp[c].s[yr], &grid_Stat[cell].grp1[yr][c]);
			if (BmassFlags.size) _copy_over(&_Gsize[c].s[yr], &grid_Stat[cell].gsize[yr][c]);
			if (BmassFlags.pr)	_copy_over(&_Gpr[c].s[yr], &grid_Stat[cell].gpr2[yr][c]);
		}
	}
  
  	if(BmassFlags.sppb) {
  		SppIndex s;
  		ForEachSpecies(s) {
  			_copy_over(&_Spp[s].s[yr], &grid_Stat[cell].spp[yr][s]);
  			if (BmassFlags.indv) _copy_over(&_Indv[s].s[yr], &grid_Stat[cell].indv[yr][s]);
  		}
  	}

	if(UseGrid && UseSeedDispersal) {
		SppIndex s;
		ForEachSpecies(s)
			_copy_over(&_Sreceived[s].s[yr], &grid_Stat[cell].sreceived[yr][s]);
	}
}

/***********************************************************/
void stat_Save_Accumulators(int cell, int year) {
	//saves the accumulators for the cell at the given year
	
	if (firsttime) {
		firsttime = FALSE;
		_init();
	}
	IntS age;
  	int yr;
  	yr = year - 1;
  	
	if(MortFlags.species) {
		SppIndex sp;
		
		ForEachSpecies(sp) {
			if ( !Species[sp]->use_me) continue;
			_copy_over(&grid_Stat[cell].sestab[sp][0], &_Sestab[sp].s[0]);
			for (age=0; age < SppMaxAge(sp); age++)
				_copy_over(&grid_Stat[cell].smort[sp][age], &_Smort[sp].s[age]);
		}
  	}
	if(MortFlags.group) {
		GrpIndex rg;

		ForEachGroup(rg) {
			if (!RGroup[rg]->use_me) continue;
			_copy_over(&grid_Stat[cell].gestab[rg][0], &_Gestab[rg].s[0]);
			for (age=0; age < GrpMaxAge(rg); age++)
				_copy_over(&grid_Stat[cell].gmort[rg][age], &_Gmort[rg].s[age]);
		}
	}
	
	if (BmassFlags.tmp) _copy_over(&grid_Stat[cell].temp[yr], &_Temp.s[yr]);
	if (BmassFlags.ppt) _copy_over(&grid_Stat[cell].ppt[yr], &_Ppt.s[yr]);
	if (BmassFlags.dist) _copy_over(&grid_Stat[cell].dist[yr], &_Dist.s[yr]);

	if(BmassFlags.grpb) {
		GrpIndex c;
		ForEachGroup(c) {
			_copy_over(&grid_Stat[cell].grp1[yr][c], &_Grp[c].s[yr]);
			if (BmassFlags.size) _copy_over(&grid_Stat[cell].gsize[yr][c], &_Gsize[c].s[yr]);
			if (BmassFlags.pr)	_copy_over(&grid_Stat[cell].gpr2[yr][c], &_Gpr[c].s[yr]);
		}
	}
  
  	if(BmassFlags.sppb) {
  		SppIndex s;
  		ForEachSpecies(s) {
  			_copy_over(&grid_Stat[cell].spp[yr][s], &_Spp[s].s[yr]);
  			if (BmassFlags.indv) _copy_over(&grid_Stat[cell].indv[yr][s], &_Indv[s].s[yr]);
  		}
  	}

	if(UseGrid && UseSeedDispersal) {
		SppIndex s;
		ForEachSpecies(s)
			_copy_over(&grid_Stat[cell].sreceived[yr][s], &_Sreceived[s].s[yr]);
	}

}

/***********************************************************/
void stat_Free_Accumulators( void ) {
	//frees all the memory allocated in stat_init_Accumulators()

  	int i, j;
  	for( i = 0; i < Globals.nCells; i++) {
  		for( j = 0; j < Globals.runModelYears; j++) {
  			if(BmassFlags.grpb) {
  				Mem_Free(grid_Stat[i].grp1[j]);
  				if (BmassFlags.size) Mem_Free(grid_Stat[i].gsize[j]);
  				if (BmassFlags.pr) Mem_Free(grid_Stat[i].gpr2[j]);
  			}
  			if(BmassFlags.sppb) {
  				Mem_Free(grid_Stat[i].spp[j]);
  				if(BmassFlags.indv) Mem_Free(grid_Stat[i].indv[j]);
  			}
			if(UseSeedDispersal && UseGrid)
				Mem_Free(grid_Stat[i].sreceived[j]);
  		}
  		
  		if (BmassFlags.dist) Mem_Free(grid_Stat[i].dist);
  		if (BmassFlags.ppt) Mem_Free(grid_Stat[i].ppt); 
  		if (BmassFlags.tmp) Mem_Free(grid_Stat[i].temp);
  		
  		if(BmassFlags.grpb) {
  			Mem_Free(grid_Stat[i].grp1);// gave grp and gpr numbers attached to them so I wouldn't mix them up lol... bad (confusing) variable names on part of the original creator.
  			if (BmassFlags.size) Mem_Free(grid_Stat[i].gsize);
  			if (BmassFlags.size) Mem_Free(grid_Stat[i].gpr2);
  		}
  		if (MortFlags.group) {
  			GrpIndex gp;
  			ForEachGroup(gp) {
  				Mem_Free(grid_Stat[i].gmort[gp]);
  				Mem_Free(grid_Stat[i].gestab[gp]);
  			}
  			Mem_Free(grid_Stat[i].gmort);
  			Mem_Free(grid_Stat[i].gestab);
  		}
  		if(BmassFlags.sppb) {
  			Mem_Free(grid_Stat[i].spp);
  			if(BmassFlags.indv) Mem_Free(grid_Stat[i].indv);
  		}
  		if (MortFlags.species) {
  			SppIndex sp;
  			ForEachSpecies(sp) {
  				Mem_Free(grid_Stat[i].smort[sp]);
  				Mem_Free(grid_Stat[i].sestab[sp]);
  			}
  			Mem_Free(grid_Stat[i].smort);
  			Mem_Free(grid_Stat[i].sestab);
  		}
		if(UseSeedDispersal && UseGrid)
			Mem_Free(grid_Stat[i].sreceived);
  	}
  	Mem_Free(grid_Stat);
  	stat_free_mem();
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
  		
  	if (BmassFlags.dist) Mem_Free(_Dist.s);
  	if (BmassFlags.ppt) Mem_Free(_Ppt.s); 
  	if (BmassFlags.tmp) Mem_Free(_Temp.s);
  		
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

/***********************************************************/
void stat_Collect_GMort ( void ) {
/* accumulated for the entire model run within
   Species_Update_Kills(), then collected
   here to compare among iterations.

   5/20/01
*/
    IntS rg, age;

    ForEachGroup(rg) {
      if (!RGroup[rg]->use_me) continue;
      _collect_add( _Gestab[rg].s, RGroup[rg]->estabs);
      for (age=0; age < GrpMaxAge(rg); age++)
        _collect_add( &_Gmort[rg].s[age],
              (double) RGroup[rg]->kills[age]);

    }

}

/***********************************************************/
void stat_Collect_SMort ( void ) {
/* accumulated for the entire model run within
   Species_Update_Kills(), then collected
   here to compare among iterations.

   5/20/01

*/
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

/***********************************************************/
void stat_Output_YrMorts( void ) {

  FILE *f = Globals.mort.fp_year;
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
  for(age=0; age < Globals.Max_Age; age++) {
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

/***********************************************************/
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
  for(age=0; age < Globals.Max_Age; age++) {
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


//This function will create individual grid cell output file and copy values to calculate for calculating grid cell avg values
/***********************************************************/
void stat_Output_AllBmassAvg() {
	  char buf[1024], tbuf[80], sep = BmassFlags.sep;
	  IntS yr;
	  GrpIndex rg;
	  SppIndex sp;
	  FILE *f;

	  if (!BmassFlags.summary) return;

	  f = OpenFile( Parm_name( F_BMassAvg), "w");

	  buf[0]='\0';

	  if (BmassFlags.header) {
	    _make_header(buf);
	    fprintf(f, "%s", buf);
	  }

	  for( yr=1; yr<= Globals.runModelYears; yr++) {
	    *buf = '\0';
	    if (BmassFlags.yr)
	      sprintf(buf, "%d%c", yr, sep);

	    if (BmassFlags.dist) {
	      sprintf(tbuf, "%ld%c", _Dist.s[yr-1].nobs, sep);
	      _Dist_grid_cell.s[yr-1].nobs = _Dist.s[yr-1].nobs;
	      strcat(buf, tbuf);
	    }

		if (BmassFlags.ppt)
		{
			RealF avg = _get_avg(&_Ppt.s[yr - 1]);
			RealF std = _get_std(&_Ppt.s[yr - 1]);
			sprintf(tbuf, "%f%c%f%c", avg, sep, std, sep);
			copyStruct(avg, std, &_Ppt_grid_cell.s[yr - 1]);
			strcat(buf, tbuf);
		}

		if (BmassFlags.pclass)
		{
			sprintf(tbuf, "\"NA\"%c", sep);
			strcat(buf, tbuf);
		}

		if (BmassFlags.tmp)
		{
			RealF avg = _get_avg(&_Temp.s[yr - 1]);
			RealF std = _get_std(&_Temp.s[yr - 1]);
			sprintf(tbuf, "%f%c%f%c", avg, sep, std, sep);
			copyStruct(avg,std, &_Temp_grid_cell.s[yr - 1]);
			strcat(buf, tbuf);
		}

		if (BmassFlags.grpb)
		{
			ForEachGroup(rg)
			{
				RealF avg = _get_avg(&_Grp[rg].s[yr - 1]);
				sprintf(tbuf, "%f%c", avg, sep);
				copyStruct(avg,0.0, &_Grp_grid_cell[rg].s[yr - 1]);
				strcat(buf, tbuf);

				if (BmassFlags.size)
				{
					RealF sizeAvg = _get_avg(&_Gsize[rg].s[yr - 1]);
					sprintf(tbuf, "%f%c",sizeAvg , sep);
					copyStruct(sizeAvg,0.0,&_Gsize_grid_cell[rg].s[yr - 1]);
					strcat(buf, tbuf);
				}

				if (BmassFlags.pr)
				{
					RealF prAvg = _get_avg(&_Gpr[rg].s[yr - 1]);
					RealF std = _get_std(&_Gpr[rg].s[yr - 1]);
					sprintf(tbuf, "%f%c%f%c", prAvg, sep,std , sep);
					copyStruct(prAvg,std,&_Gpr_grid_cell[rg].s[yr - 1]);
					strcat(buf, tbuf);
				}
			}
		}


		if (BmassFlags.sppb)
		{
			ForEachSpecies(sp)
			{
				RealF spAvg = _get_avg(&_Spp[sp].s[yr - 1]);
				sprintf(tbuf, "%f%c", spAvg, sep);
				copyStruct(spAvg,0.0,&_Spp_grid_cell[sp].s[yr - 1]);
				strcat(buf, tbuf);

				if (BmassFlags.indv)
				{
					RealF indvAvg = _get_avg(&_Indv[sp].s[yr - 1]);
					sprintf(tbuf, "%f%c", indvAvg, sep);
					copyStruct(indvAvg,0.0,&_Indv_grid_cell[sp].s[yr - 1]);
					strcat(buf, tbuf);
				}

			}
		}

	    fprintf( f, "%s\n", buf);
	  }  /* end of foreach year */
	  CloseFile(&f);

}

//This function will create grid cell avg values output file
void stat_Output_AllCellAvgBmass(const char * filename)
{

	char buf[1024], tbuf[80], sep = BmassFlags.sep;
	IntS yr;
	GrpIndex rg;
	SppIndex sp;
	FILE *f;

	if (!BmassFlags.summary)
		return;

	f = OpenFile(filename, "w");

	buf[0] = '\0';

	if (BmassFlags.header)
	{
		_make_header(buf);
		fprintf(f, "%s", buf);
	}

	for (yr = 1; yr <= Globals.runModelYears; yr++)
	{

		*buf = '\0';
		if (BmassFlags.yr)
			sprintf(buf, "%d%c", yr, sep);

		if (BmassFlags.dist)
		{
			sprintf(tbuf, "%ld%c", _Dist_grid_cell.s[yr - 1].nobs, sep);
			strcat(buf, tbuf);
		}

		if (BmassFlags.ppt)
		{
			sprintf(tbuf, "%f%c%f%c",
					_get_gridcell_avg(&_Ppt_grid_cell.s[yr - 1]), sep,
					_get_gridcell_std(&_Ppt_grid_cell.s[yr - 1]), sep);
			strcat(buf, tbuf);
		}

		if (BmassFlags.pclass)
		{
			sprintf(tbuf, "\"NA\"%c", sep);
			strcat(buf, tbuf);
		}

		if (BmassFlags.tmp)
		{
			sprintf(tbuf, "%f%c%f%c",
					_get_gridcell_avg(&_Temp_grid_cell.s[yr - 1]), sep,
					_get_gridcell_std(&_Temp_grid_cell.s[yr - 1]), sep);
			strcat(buf, tbuf);
		}

		if (BmassFlags.grpb)
		{
			ForEachGroup(rg)
			{
				sprintf(tbuf, "%f%c", _get_gridcell_avg(&_Grp_grid_cell[rg].s[yr - 1]), sep);
				strcat(buf, tbuf);

				if (BmassFlags.size)
				{
					sprintf(tbuf, "%f%c", _get_gridcell_avg(&_Gsize_grid_cell[rg].s[yr - 1]), sep);
					strcat(buf, tbuf);
				}

				if (BmassFlags.pr)
				{
					sprintf(tbuf, "%f%c%f%c",
							_get_gridcell_avg(&_Gpr_grid_cell[rg].s[yr - 1]), sep,
							_get_gridcell_std(&_Gpr_grid_cell[rg].s[yr - 1]), sep);
					strcat(buf, tbuf);
				}
			}
		}

		if (BmassFlags.sppb)
		{
			ForEachSpecies(sp)
			{
		        sprintf(tbuf, "%f%c",_get_gridcell_avg( &_Spp_grid_cell[sp].s[yr-1]), sep);
		        strcat( buf, tbuf);

				if (BmassFlags.indv)
				{
		          sprintf(tbuf, "%f%c", _get_gridcell_avg( &_Indv_grid_cell[sp].s[yr-1]), sep);
		          strcat( buf, tbuf);
				}

			}
		}

		fprintf(f, "%s\n", buf);
	} /* end of foreach year */
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
	_make_header_with_std(buf);
    fprintf(f, "%s", buf);
  }

  for( yr=1; yr<= Globals.runModelYears; yr++) {
    *buf = '\0';
    if (BmassFlags.yr)
      sprintf(buf, "%d%c", yr, sep);

    if (BmassFlags.dist) {
      sprintf(tbuf, "%ld%c", _Dist.s[yr-1].nobs,
              sep);
      strcat(buf, tbuf);
    }

    if (BmassFlags.ppt) {
      sprintf(tbuf, "%f%c%f%c",
              _get_avg(&_Ppt.s[yr-1]), sep,
              _get_std(&_Ppt.s[yr-1]), sep);
      strcat( buf, tbuf);
    }

    if (BmassFlags.pclass) {
      sprintf(tbuf, "\"NA\"%c", sep);
      strcat( buf, tbuf);
    }

    if (BmassFlags.tmp) {
      sprintf(tbuf, "%f%c%f%c",
              _get_avg(&_Temp.s[yr-1]), sep,
              _get_std(&_Temp.s[yr-1]), sep);
      strcat( buf, tbuf);
    }

    if (BmassFlags.grpb) {
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
      }
    }

		if (BmassFlags.sppb)
		{
			for ((sp) = 0; (sp) < Globals.sppCount - 1; (sp)++)
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


/***********************************************************/
void stat_Output_Seed_Dispersal(const char * filename, const char sep, Bool makeHeader) {
	//do stuff...
	char buf[1024], tbuf[80];
	IntS yr;
	SppIndex sp;
	FILE *f;

	f = OpenFile(filename, "w");

	if(makeHeader) {
		fprintf(f,"Year");
		ForEachSpecies(sp) {
			fprintf(f, "%c%s_prob", sep, Species[sp]->name);
			fprintf(f, "%c%s_std", sep, Species[sp]->name);
		}
		fprintf(f,"\n");
	}

	for( yr=1; yr<= Globals.runModelYears; yr++) {
		*buf = '\0';
		
		sprintf(buf, "%d%c", yr, sep);
		
		ForEachSpecies(sp) {
			sprintf(tbuf, "%f%c%f%c", _get_avg( &_Sreceived[sp].s[yr-1]), sep, _get_std( &_Sreceived[sp].s[yr-1]), sep);
			strcat(buf, tbuf);
		}

		fprintf(f, "%s\n", buf);
	}	
	CloseFile(&f);
}


/***********************************************************/
static RealF _get_avg( struct accumulators_st *p) {

	if (p->nobs == 0) return 0.0;

	return (RealF) (p->sum / (double) p->nobs);

}


/***********************************************************/
static RealF _get_std(struct accumulators_st *p)
{
	double s;

	/* Change: Jul 30 2016 (AKT)
	 * Fixing bug at std_dev calculation, as this function was giving wrong result.
	 * Whenever intermediate calcuation value of s coming negative result become NAN.
	 * So added check for negative number and floating point calculation.
	 * Also this std_dev calculation function is a derivation of actual mathematical std_dev calculation function,
	 * however while dividing intermediate s value with ( p->nobs * (p->nobs -1)) was wrong,
	 * it should be dividing by square of p->nobs (no of iteration).
	 */

	if (p->nobs <= 1)
		return 0.0;

	s = (p->nobs * p->sum_sq) - (p->sum * p->sum);

	if ((int) s <= 0)
		return 0.0;

	s = s / (p->nobs * p->nobs);

	return (RealF) sqrt(s);

}



static void copyStruct(RealF val,RealF std_val,struct accumulators_grid_cell_st *p)
{
	p->sum = p->sum + val;
	p->sum_std = p->sum_std + std_val;

}

/***********************************************************/
static RealF _get_gridcell_avg(struct accumulators_grid_cell_st *p)
{

	if (Globals.nCells == 0)
		return 0.0;
	RealF avg = (RealF) (p->sum / (double) Globals.nCells);
	return avg;
}


/***********************************************************/
static RealF _get_gridcell_std(struct accumulators_grid_cell_st *p)
{
	if (Globals.nCells == 0)
			return 0.0;
		RealF avg = (RealF) (p->sum_std / (double) Globals.nCells);
		return avg;
}



/***********************************************************/
static void _make_header_with_std( char *buf) {

  char fields[MAX_OUTFIELDS*2][MAX_FIELDLEN+1];
  char tbuf[80];
  GrpIndex rg;
  SppIndex sp;
  Int i, fc=0;


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


}

/***********************************************************/
static void _make_header( char *buf) {

  char fields[MAX_OUTFIELDS*2][MAX_FIELDLEN+1];
  char tbuf[80];
  GrpIndex rg;
  SppIndex sp;
  Int i, fc=0;


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
