/**
 * \file ST_params.c
 * \brief Reads and initializes the model parameters. 
 * 
 * Most of the parameters come from the input files and some are computed.
 * 
 * History
 * (6/15/2000) -- INITIAL CODING - cwb
 * (15-Apr-02)  -- added code to interface with SOILWAT (cwb)
 *                 only modified parm_Files_Init().
 * 
 * \author CWB (initial coding)
 * \author Chandler Haukap (author of this documentation)
 * \date 15 April 2002
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "ST_steppe.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sw_src/rands.h"
#include "sxw_funcs.h"
#include "ST_globals.h"
#include "sxw_vars.h"

/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
SppIndex species_New(void);

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
  void parm_Initialize(void);
  void parm_SetFirstName( char *s);
  void parm_SetName( char *s, int which);
  void parm_free_memory( void );
  void maxrgroupspecies_init(void);
  void files_init(void);

/*********** Locally Used Function Declarations ************/
/***********************************************************/
static void _env_init( void);
static void _plot_init( void);
static void _setNameLen(char *dest, char *src, Int len);
static void _rgroup_init( void);
static void _species_init( void);
static void _check_species( void);
static void _bmassflags_init( void);
static void _mortflags_init( void);
static void _model_init( void);
static void _rgroup_add1( char name[], RealF space, RealF density,
                      Int estab, RealF slow, Int stretch,
                      Int xres, Int estann, Int turnon,
                      Int styr,  RealF xgro, Int veg_prod_type, Int mort,
                      RealF biomass, RealF transpiration);
static void _rgroup_add2( char name[],
                      RealF nslope, RealF nint,
                      RealF wslope, RealF wint,
                      RealF dslope, RealF dint);
static void _rgroup_add_disturbance( char name[],  Int killyr, Int killfreq_startyr,RealF killfreq,
                      Int extirp, RealF prop_killed, RealF prop_recovered,RealF grazing_frq,RealF prop_grazing,Int grazingfreq_startyr);
static void _rgroup_addsucculent( char name[],
                               RealF wslope, RealF wint,
                               RealF dslope, RealF dint);

//static void _recover_names(void);

/************ Module Variable Declarations ******************/
/***********************************************************/
  #define NFILES 15

static char *_files[NFILES];
char *MyFileName;

/**************************************************************/
/* fdpierson: This function should only be called once, otherwise
 * memory leaks involving resource groups, species, and potentially
 * other variables will occur. */
/**************************************************************/
void parm_Initialize() {
/*======================================================*/
  _model_init();
  _env_init();
  _plot_init();
  _bmassflags_init();
  _mortflags_init();
  _rgroup_init();
  _species_init();
  _check_species();
}

/**************************************************************/
char *Parm_name( ST_FileIndex i) {
/*======================================================*/

  return _files[i];
}

/**************************************************************/
void parm_SetFirstName( char *s) {
/*======================================================*/

  if (_files[F_First]) Mem_Free( _files[F_First]);
  _files[F_First] = Str_Dup(s);

}

/**************************************************************/
void parm_SetName( char *s, int which) {
/*======================================================*/
// (DLM) - 5/28/2013 - added this function because I needed a way to change the name of some of the files (particularly the output files) easily from other modules

	if(which > (NFILES - 1) || which < 0)
		return;
	if (_files[which]) Mem_Free( _files[which]);
  	_files[which] = Str_Dup(s);

}

/**************************************************************/
void files_init( void ) {
/*======================================================*/
  /* 7-May-02 (cwb) added code to interface with SOILWAT */

  FILE *f;
  ST_FileIndex i;
  ST_FileIndex last = F_MaxRGroupSpecies;

  MyFileName = Parm_name(F_First);

  f = OpenFile(MyFileName, "r");

  for(i=F_Log; i <= last; i++) {
    if ( ! GetALine(f, inbuf)) break;
    _files[i] = Str_Dup(Str_TrimLeftQ(inbuf));
    //printf("FILES: %d : %s\n", i, _files[i]);
  }

  if ( i < last) {
    LogError(logfp, LOGFATAL, "%s: Too few input files specified",
                              MyFileName);
  }

  CloseFile(&f);

  if(!UseGrid) {	//if the gridded option has been specified, then the logfile has already been opened
  	if ( !strcmp("stdout", _files[F_Log]) )
    	logfp = stdout;
  	else
    	logfp = OpenFile(_files[F_Log], "w");
   }

}

/**************************************************************/
static void _model_init( void) {
/*======================================================*/
   FILE *f;
   int seed, years;
   char tmp[80];

   MyFileName = Parm_name(F_Model);
   f = OpenFile(MyFileName, "r");
   /* ----------------------------------------------------*/
   /* scan for the first line*/
   if (!GetALine(f, inbuf)) {
     LogError(logfp, LOGFATAL, "%s: No data found!\n", MyFileName);
   } else {
     sscanf( inbuf, "%s %d %d",
             tmp,
             &years,
             &seed);

     SuperGlobals.runModelYears = years;
     Globals->Max_Age = SuperGlobals.runModelYears;
     SuperGlobals.runModelIterations = atoi(tmp);
     if (SuperGlobals.runModelIterations < 1 ||
         SuperGlobals.runModelYears < 1 ) {
       LogError(logfp, LOGFATAL,"Invalid parameters for RunModelIterations "
               "or RunModelYears (%s)",
               MyFileName);
     }
     Globals->bmass.suffixwidth = Globals->mort.suffixwidth = strlen(tmp);
     SuperGlobals.randseed = (IntL) ((seed) ? -abs(seed) : 0);
   }
   
   CloseFile(&f);
}


/**************************************************************/
static void _env_init( void) {
/*======================================================*/
/* Read environmental parameters that drive
   the abiotic conditions.  These values go
   into the Globals structure
*/
   FILE *f;
   int x,
       index=0,
       nitems,
       use[3];

   MyFileName = Parm_name(F_Env);
   f = OpenFile(MyFileName, "r");

   while( GetALine(f, inbuf)) {

      switch(++index) {
        case 1:
            x=sscanf( inbuf, "%f %f %hu %hu %hu %hu %f %hu",
                      &Globals->ppt.avg, &Globals->ppt.std,
                      &Globals->ppt.min, &Globals->ppt.max,
                      &Globals->ppt.dry, &Globals->ppt.wet,
                      &Globals->gsppt_prop, &Globals->transp_window);
            nitems = 8;
            break;
        case 2:
            x=sscanf( inbuf, "%f %f %f %f %f",
                      &Globals->temp.avg, &Globals->temp.std,
                      &Globals->temp.min, &Globals->temp.max,
                      &Globals->temp.gstemp);
            nitems = 5;
            break;
        case 3:
            x=sscanf( inbuf, "%d %f %f %f %f ", &use[0],
                      &Globals->pat.occur, &Globals->pat.removal,
                      &Globals->pat.recol[Slope],
                      &Globals->pat.recol[Intcpt]);
            nitems = 4;
            break;
        case 4:
            x=sscanf( inbuf, "%d %f %hu %hu", &use[1],
                      &Globals->mound.occur,
                      &Globals->mound.minyr,
                      &Globals->mound.maxyr);
            nitems = 3;
            break;
        case 5:
            x=sscanf( inbuf, "%d %f %hu", &use[2],
                      &Globals->burrow.occur,
                      &Globals->burrow.minyr);
            nitems = 2;
            break;
        case 6:
            x=sscanf( inbuf, "%f %f %f %f %f %f",
                      &Globals->tempparm[CoolSeason][0],
                      &Globals->tempparm[CoolSeason][1],
                      &Globals->tempparm[CoolSeason][2],
                      &Globals->tempparm[WarmSeason][0],
                      &Globals->tempparm[WarmSeason][1],
                      &Globals->tempparm[WarmSeason][2]);
            nitems = 6;
            break;
      }
      if (x<nitems) {
         LogError(logfp, LOGFATAL, "%s: Invalid record %d",
                 MyFileName, index);
      }

   } /* end while*/

   Globals->pat.use    = itob(use[0]);
   Globals->mound.use  = itob(use[1]);
   Globals->burrow.use = itob(use[2]);

   CloseFile(&f);
}


/**************************************************************/
static void _plot_init( void) {
/*======================================================*/

   FILE *f;
   int x, nitems=1;

   MyFileName = Parm_name(F_Plot);
   f = OpenFile(MyFileName, "r");

   /* ----------------------------------------------------*/
   /* scan for the first line*/
   if (!GetALine(f, inbuf)) {
     LogError(logfp, LOGFATAL, "%s: No data found!\n", MyFileName);
   }

   x = sscanf( inbuf, " %f", &Globals->plotsize);
   if (x < nitems) {
     LogError(logfp, LOGFATAL, "%s: Incorrect number of fields",
                     MyFileName);
   }


   CloseFile(&f);
}


/**************************************************************/
static void _check_species( void) {
/*======================================================*/
  /* Make sure max_spp_estab <= max_spp */
  /* Also this block creates then links */
  /* species' numbers to their respective groups. */
  /* - and sets up mortality arrays for groups and spp (5/18/01)*/
  /* 10-Apr-03 - cwb - compute max_bmass for the group as the
   *      sum of the mature biomass of each species in the group.
   *      plot biomass = sum(g->relsize * g->max_density * g->max_bmass)
   * 7-Nov-03 (cwb) make sure no annuals are assigned to
   *      perennials groups.
   */
  IntS i, cnt,
       maxage, minage, /* placeholders */
       runyrs = SuperGlobals.runModelYears;  /* shorthand */
  SppIndex sp;
  GrpIndex rg;
  Bool tripped = FALSE;
  GroupType *g;
  SpeciesType *s;

  /* -------------------------------------------*/
   /* count and link species to their groups.
   * print a message if more specified than available */
  ForEachGroup(rg) { g = RGroup[rg];
    cnt = 0;
    ForEachSpecies(sp) {
      if (Species[sp]->res_grp == rg)
        g->species[cnt++] = sp;
    }
    g->max_spp = cnt;
    if (cnt < g->max_spp_estab) {
      tripped = TRUE;
      g->max_spp_estab = cnt;
      LogError(logfp, LOGNOTE, "Max_Spp_Estab > Number of Spp for %s",
              g->name);
    }
  }
  if (tripped) LogError(logfp, LOGNOTE,"Continuing.");

  /* -------------------------------------------*/
  /* determine max age for the species and
   * keep track for group's max age */
  /* and compute max g/m^2 for this group */
  /* 10-Apr-03 - also, compute max_bmass for the group. */
  ForEachGroup(rg) { g = RGroup[rg];
    maxage = 0; minage = 30000;
    g->max_bmass = 0.0;
    ForEachGroupSpp(sp, rg, i) {
      s = Species[sp];
      if (s->max_age == 0) s->max_age = runyrs+1;
/* maxage shouldn't be set to extirp due to age-independent mortality.
   extirp happens due to the flag, not age;
   but I hesitate to remove the code, so I'll comment it out for now.
      if ( g->extirp && g->extirp < s->max_age)
        s->max_age = g->extirp;
*/
      maxage = max(s->max_age, maxage);
      minage = min(s->max_age, minage);
      g->max_bmass += s->mature_biomass;
    }

    if (minage == 1 && maxage != 1) {
      LogError(logfp, LOGFATAL, "%s: Can't mix annuals and perennials within a group\n"
                      "Refer to the groups.in and species.in files\n",
                       RGroup[rg]->name);
    }
    g->max_age = maxage;
  }

  /* -------------------------------------------*/
  /* check out the definitions for SppMaxAge and GrpMaxAge */
  /* they're used here for some hoped-for sense of readability */

  if (MortFlags.species) {
    /* size each species kills list, age ==0 ok */
    ForEachSpecies(sp) {
      if ( Species[sp]->use_me) {
        Species[sp]->kills = (IntUS *) Mem_Calloc(SppMaxAge(sp),
                                               sizeof(IntUS),
                             "_check_species(Species.kills)");
      } else {
        Species[sp]->kills = NULL;
      }
    }
  } else {
    ForEachSpecies(sp)
      Species[sp]->kills = NULL;
  }

  if (MortFlags.group) {
    ForEachGroup(rg) {
      if (RGroup[rg]->use_me) {
        /* find max age of the group and size the kills array */
        ForEachGroupSpp( sp, rg, i) {
          RGroup[rg]->max_age = (Species[sp]->max_age)
                                ? max(Species[sp]->max_age,
                                      RGroup[rg]->max_age)
                                : Globals->Max_Age;
        }
        RGroup[rg]->kills = (IntUS *) Mem_Calloc(GrpMaxAge(rg),
                                                sizeof(IntUS),
                        "_check_species(RGroup.kills)");
      } else {
        RGroup[rg]->kills = NULL;
      }
    }
  } else {
    ForEachGroup(rg)
      RGroup[rg]->kills = NULL;
  }
}


/**************************************************************/
static void _bmassflags_init( void) {
/*======================================================*/
/* read in the flags to control the output quantities etc.
   It helps to have the file visible while reading this code.

 * 9-Dec-03 (cwb) added code to create output directory.

*/

   FILE *fin;
   Int x, i,
       nitems=16; /* number of items expected in first input line */

   /*   code      controls: */
   char u[5],  /* summary? if 'n' don't init and don't print */
        a[5],  /* print yearly files */
        h[5],  /* header line containing var names */
        f[5],  /* field separator */
        y[5],  /* year number */
        d[5],  /* disturbance type */
        p[5],  /* precip */
        c[5],  /* ppt class (wet/dry/norm) */
        t[5],  /* temperature */
        g[5],  /* biomass for the groups */
        q[5],  /* groups PR quotient (see C&L1990, p240) */
        r[5],  /* relative size for each group */
        w[5],  /* wildfire count */
        m[5],  /* prescribed fire count */
        s[5],  /* biomass for each species */
        n[5];  /* number of individuals for each species */
   char z;

   MyFileName = Parm_name(F_BMassFlag);
   fin = OpenFile(MyFileName, "r");

   if (!GetALine(fin, inbuf)) {
     LogError(logfp, LOGFATAL, "%s: No data found!\n", MyFileName);
   }

   x = sscanf( inbuf, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
                      u, a, h, f, y, d, p, c, t, g, q, r, w, m, s, n );

   /* don't bother initializing the rest if first flag is 'n' */
   BmassFlags.summary  = (Bool)(*u=='y'||*u=='Y');
   BmassFlags.yearly   = (Bool)(*a=='y'||*a=='Y');
   if (!(BmassFlags.summary || BmassFlags.yearly)) {
     return;
   }

   if (x < nitems) {
     LogError(logfp, LOGFATAL, "%s: Invalid number of parameters",
             MyFileName);
   }

   for( i=1; i<= nitems; i++) {
   /* to stick in a new item, just add it to the bottom
      and adjust nitems in definitions. */
     switch (i) {
       case 1:
            BmassFlags.header = (Bool)(*h=='y'||*h=='Y');
              break;
       case 2:
            z = (char) tolower((Int)*f);
            switch( z )  {
               case 't': BmassFlags.sep = '\t'; break;
               case 's': BmassFlags.sep = ' '; break;
               default : BmassFlags.sep = z; break;
            }
            break;
       case 3:
            BmassFlags.yr     = (Bool)(*y=='y'||*y=='Y');
            break;
       case 4:
            BmassFlags.dist   = (Bool)(*d=='y'||*d=='Y');
            break;
       case 5:
            BmassFlags.ppt    = (Bool)(*p=='y'||*p=='Y');
            break;
       case 6:
            BmassFlags.pclass = (Bool)(*c=='y'||*c=='Y');
           break;
       case 7:
            BmassFlags.tmp    = (Bool)(*t=='y'||*t=='Y');
            break;
       case 8:
            BmassFlags.grpb   = (Bool)(*g=='y'||*g=='Y');
            break;
       case 9:
            BmassFlags.pr     = (Bool)(*q=='y'||*q=='Y');
            break;
       case 10:
            BmassFlags.size   = (Bool)(*r=='y'||*r=='Y');
            break;
       case 11:
            BmassFlags.wildfire     = (Bool)(*w=='y'||*w=='Y');
            break;
       case 12:
            BmassFlags.prescribedfire     = (Bool)(*m=='y'||*m=='Y');
            break;
       case 13:
            BmassFlags.sppb   = (Bool)(*s=='y'||*s=='Y');
            break;
       case 14:
            BmassFlags.indv   = (Bool)(*n=='y'||*n=='Y');
            break;
       case 15:
            break;
     }

   }
   CloseFile(&fin);

   /* remove old output and/or create the output directories if needed */
   /* borrow inbuf for filenames */
   /* -- do avg file first, otherwise it may get deleted by match
    *    with BMassPre and then not be there for specific delete.
    */
   if (DirExists(DirName(Parm_name(F_BMassAvg)))) {
     strcpy(inbuf, Parm_name(F_BMassAvg));
     if (!RemoveFiles(inbuf) )
       LogError(logfp, LOGWARN, "Can't remove old average biomass output file %s\n%s",
                inbuf, strerror(errno) );

   } else if (!MkDir(DirName(Parm_name(F_BMassAvg))) ) {
     LogError(logfp, LOGFATAL,
              "Can't make output path for average biomass file: %s\n%s",
              DirName(Parm_name(F_BMassAvg)), strerror(errno));
   }

   if (DirExists(DirName(Parm_name(F_BMassPre)))) {
     strcpy(inbuf, Parm_name(F_BMassPre));
     strcat(inbuf, "*.csv");
     if (!RemoveFiles(inbuf) )
       LogError(logfp, LOGWARN, "Can't remove old biomass output files %s\n%s",
                inbuf, strerror(errno) );

   } else if (!MkDir(DirName(Parm_name(F_BMassPre))) ) {
       LogError(logfp, LOGFATAL,
                "Can't make output path for yearly biomass files: %s\n%s",
                DirName(Parm_name(F_BMassPre)), strerror(errno) );
   }

}

/**************************************************************/
static void _mortflags_init( void) {
/*======================================================*/
/* read in the flags to control the output quantities etc.
   It helps to have the file visible while reading this code.
*/

   FILE *fin;
   Int x, i,
       nitems=6; /* number of items expected in first input line */

   /*   code      controls: */
   char s[5],  /* summary stats? if 'n' don't init and don't print */
        y[5],  /* print yearly files */
        h[5],  /* header line containing var names */
        f[5],  /* field separator */
        g[5],  /* group data */
        k[5];  /* species data */
   char z;


   MyFileName = Parm_name(F_MortFlag);
   fin = OpenFile(MyFileName, "r");

   if (!GetALine(fin, inbuf)) {
     LogError(logfp, LOGFATAL, "%s No data found!\n", MyFileName);

   }

   x = sscanf( inbuf, "%s %s %s %s %s %s",
                      s, y, h, f, g, k );

   /* don't bother initializing the rest if first flag is 'n' */
   MortFlags.summary = (Bool)(*s=='y'||*s=='Y');
   MortFlags.yearly  = (Bool)(*y=='y'||*y=='Y');
   if (!(MortFlags.summary || MortFlags.yearly)) {
     MortFlags.header = MortFlags.group = MortFlags.species = FALSE;
     return;
   }

   if (x < nitems -2) {
     LogError(logfp, LOGFATAL,"%s: Invalid number of parameters",
             MyFileName);
   }

   for( i=3; i<= nitems ; i++) {  /* 3 accounts for first two flags */
   /* to stick in a new item, just add it to the bottom
      and adjust nitems in definitions. */
     switch (i) {
       case 3:
            MortFlags.header = (Bool)(*h=='y'||*h=='Y');
              break;
       case 4:
            MortFlags.group  = (Bool)(*g=='y'||*g=='Y');
              break;
       case 5:
            MortFlags.species = (Bool)(*k=='y'||*k=='Y');
              break;
       case 6:
            z = (char) tolower((Int)*f);
            switch( z )  {
               case 't': MortFlags.sep = '\t'; break;
               case 's': MortFlags.sep = ' '; break;
               default : MortFlags.sep = z; break;
            }
            break;
       default:
            break;
     }

    }
    CloseFile(&fin);

   /* remove old output and/or create the output directories if needed */
   /* borrow inbuf for filenames */
   /* -- do avg file first, otherwise it may get deleted by match
    *    with MortPre and then not be there for specific delete.
    */
   if (DirExists(DirName(Parm_name(F_MortAvg)))) {
     strcpy(inbuf, Parm_name(F_MortAvg));
     if (!RemoveFiles(inbuf) )
       LogError(logfp, LOGWARN, "Can't remove old average biomass output file %s\n%s",
                inbuf, strerror(errno) );

   } else if (!MkDir(DirName(Parm_name(F_MortAvg))) ) {
     LogError(logfp, LOGFATAL,
              "Can't make output path for average biomass file: %s\n%s",
              DirName(Parm_name(F_MortAvg)), strerror(errno));
   }


   if (DirExists(DirName(Parm_name(F_MortPre)))) {
     strcpy(inbuf, Parm_name(F_MortPre));
     strcat(inbuf, "*.csv");
     if (!RemoveFiles(inbuf) )
       LogError(logfp, LOGWARN, "Can't remove old biomass output files %s\n%s",
                inbuf, strerror(errno) );

   } else if (!MkDir(DirName(Parm_name(F_MortPre))) ) {
       LogError(logfp, LOGFATAL,
                "Can't make output path for yearly biomass files: %s\n%s",
                DirName(Parm_name(F_MortPre)), strerror(errno) );
   }

}

/**************************************************************/
static void _setNameLen(char *dest, char *src, Int len) {
/*======================================================*/
  strncpy(dest, src, len);
  dest[len] = '\0';
}

/**************************************************************/
void maxrgroupspecies_init( void) {
/*======================================================*/
    FILE *f;
    
    MyFileName = Parm_name(F_MaxRGroupSpecies);
    f = OpenFile(MyFileName, "r");
    
    /* These values determine the memory allocated for resource groups, species, and their names. */
    /* If these limits are exceeded, memory leaks will result. */

    /* Resource group limits */

    if (!GetALine(f, inbuf)) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum resource groups allowed.", MyFileName);
    }

    if (sscanf(inbuf, "%zu", &SuperGlobals.max_rgroups) != 1) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum resource groups allowed.", MyFileName);
    }

    if (!GetALine(f, inbuf)) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum resource group name length.", MyFileName);
    }

    if (sscanf(inbuf, "%zu", &SuperGlobals.max_groupnamelen) != 1) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum resource group name length.", MyFileName);
    }

    /* Species limits */

    if (!GetALine(f, inbuf)) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum species allowed per resource group.", MyFileName);
    }

    if (sscanf(inbuf, "%zu", &SuperGlobals.max_spp_per_grp) != 1) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum species allowed per resource group.", MyFileName);
    }

    if (!GetALine(f, inbuf)) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum individuals allowed per species.", MyFileName);
    }

    if (sscanf(inbuf, "%zu", &SuperGlobals.max_indivs_per_spp) != 1) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum individuals allowed per species.", MyFileName);
    }

    if (!GetALine(f, inbuf)) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum species name length.", MyFileName);
    }

    if (sscanf(inbuf, "%zu", &SuperGlobals.max_speciesnamelen) != 1) {
       LogError(logfp, LOGFATAL, "%s: Could not read maximum species name length.", MyFileName);
    }
    
    CloseFile(&f);
}


/**************************************************************/
/**************************************************************/
/*****************************************************
 * The *_RGroup_* functions read the input from the user
 * file for resource-group-level information.
 *
 *****************************************************/
static void _rgroup_init( void) {
/*======================================================*/
/* Read parameters for resource group
*/
   FILE *f;
   IntS x;
   Bool groupsok;

   /* temp vars to hold the group info*/
   char *name;

   /* input variables*/
   Int estab, stretch, xres, turnon, estann,
       styr, veg_prod_type, mort;
   RealF space, density, slow, 
         nslope, nint, wslope, wint, dslope, dint, xgro;
   /* input variables related to disturbances */
   Int extirp, killyr, killfreq_startyr, 
       grazingfreq_startyr;
   RealF  killfreq, prop_killed, prop_recovered,grazing_frq, prop_grazing, biomass,
        transpiration;

   MyFileName = Parm_name(F_RGroup);
   f = OpenFile(MyFileName, "r");
   
   name = (char *)Mem_Calloc(SuperGlobals.max_groupnamelen + 1, sizeof(char), "_rgroup_init");
   RGroup = (GroupType **)Mem_Calloc(SuperGlobals.max_rgroups, sizeof(GroupType *), "_rgroup_init");

   /* ------------------------------------------------------------*/
   /* Install all the defined groups, except for dry/wet/norm parms */
   groupsok = FALSE;
   while( GetALine(f,inbuf)) {
     if (!isnull(strstr(inbuf,"[end]"))) {
        groupsok = TRUE;
        break;
     }
     x=sscanf( inbuf, "%s %f %f %d %f %d %d %d %d %d %f %d %d %d %f %d %d %f %f"
               " %f %f %d %f %f", name, &space, &density, &estab, &slow, 
               &stretch, &xres, &estann, &turnon, &styr, &xgro, &veg_prod_type,
               &killyr, &killfreq_startyr, &killfreq, &extirp, &mort, 
               &prop_killed, &prop_recovered,&grazing_frq, &prop_grazing,
               &grazingfreq_startyr, &biomass, &transpiration);
     if (x < 22) {
       LogError(logfp, LOGFATAL, "%s: Too few columns in groups",
               MyFileName);
     }

    // Convert to SOILWAT2 vegetation index
    veg_prod_type = get_SW2_veg_index(veg_prod_type);

     _rgroup_add1( name, space, density, estab, slow, stretch, xres, estann, 
                   turnon, styr, xgro, veg_prod_type, mort, biomass, 
                   transpiration);

     _rgroup_add_disturbance(name, killyr, killfreq_startyr, killfreq,
                   extirp, prop_killed, prop_recovered,grazing_frq, prop_grazing, 
                   grazingfreq_startyr);
   }/* end while*/

   if (!groupsok) {
      LogError(logfp, LOGFATAL, "%s: Incomplete input in group definitions",
              MyFileName);
   }

   /* ------------------------------------------------------------*/
   /* Install  dry/wet/norm parms for defined groups */
   groupsok = FALSE;
   while( GetALine(f,inbuf)) {
     if (!isnull(strstr(inbuf,"[end]"))) {
        groupsok = TRUE;
        break;
     }
     x=sscanf( inbuf, "%s %f %f %f %f %f %f",
               name,
               &nslope, &nint, &wslope, &wint, &dslope, &dint);
     if (x != 7) {
       LogError(logfp, LOGFATAL, "%s: Wrong number of columns in groups' wet/dry parms",
               MyFileName);
     }
     _rgroup_add2( name, nslope, nint, wslope, wint, dslope, dint);
   }/* end while*/

   if (!groupsok) {
      LogError(logfp, LOGFATAL, "%s: Incomplete input in group definitions",
              MyFileName);
   }

   /* ------------------------------------------------------------*/
   /* Get succulent growth modifiers*/
   GetALine(f, inbuf);
   x=sscanf( inbuf, "%s %f %f %f %f",
             name, &wslope, &wint, &dslope, &dint);
   if (x < 5) {
     LogError(logfp, LOGFATAL,
          "%s: Too few values in succulent growth parameters",
           MyFileName);
   }
   _rgroup_addsucculent( name, wslope, wint, dslope, dint);

   /* ------------------------------------------------------------ */
   /* Get wildfire parameters */
   GetALine(f,inbuf);

   groupsok = FALSE;
   while(GetALine(f,inbuf)) {
     if (!isnull(strstr(inbuf,"[end]"))) {
        groupsok = TRUE;
        
        break;
     }

     x=sscanf( inbuf, "%u", &UseCheatgrassWildfire);
     if (x != 1) {
       LogError(logfp, LOGFATAL, "%s: Cheatgrass-Wildfire flag not read.",
               MyFileName);
     } 
   }/* end while*/

   Mem_Free(name);
   
   CloseFile(&f);
}

/**************************************************************/
static void _rgroup_add1( char name[], RealF space, RealF density,
                      Int estab, RealF slow, Int stretch,
                      Int xres, Int estann, Int turnon,
                      Int styr,  RealF xgro, Int veg_prod_type, Int mort,
                      RealF biomass, RealF transpiration) {
  GrpIndex rg;
  size_t len;
  
  rg = RGroup_New();
  
  len = strlen(name);
  
  _setNameLen(RGroup[rg]->name, name, len);
  RGroup[rg]->grp_num       = rg;
  RGroup[rg]->max_stretch   = (IntS) stretch;
  RGroup[rg]->max_spp_estab = (IntS) estab;
  // input of `density` is in units of [# / m2]; convert to units of [# / plot]
  RGroup[rg]->max_density   = density * Globals->plotsize; // density per plot
  RGroup[rg]->max_per_sqm   = density; // density per square-meter
  RGroup[rg]->use_mort      = itob(mort);
  RGroup[rg]->slowrate      = slow;
  RGroup[rg]->space = space;
  RGroup[rg]->min_res_req   = space;
  RGroup[rg]->est_annually  = itob(estann);
  RGroup[rg]->startyr       = styr;
  RGroup[rg]->xgrow         = xgro;
  RGroup[rg]->use_me        = itob(turnon);
  RGroup[rg]->veg_prod_type = veg_prod_type;
  RGroup[rg]->use_extra_res = itob(xres);
  RGroup[rg]->_bvt = biomass / transpiration;
}


/**************************************************************/
static void _rgroup_add2( char name[], RealF nslope, RealF nint, RealF wslope,
                          RealF wint, RealF dslope, RealF dint) {
  GrpIndex rg;

  rg = RGroup_Name2Index(name);
  if (rg < 0) {
    LogError(logfp, LOGFATAL, "%s: Mismatched name (%s) for succulents",
             MyFileName, name);
  }

  RGroup[rg]->ppt_slope[Ppt_Norm] = nslope;
  RGroup[rg]->ppt_intcpt[Ppt_Norm]= nint;
  RGroup[rg]->ppt_slope[Ppt_Wet]  = wslope;
  RGroup[rg]->ppt_intcpt[Ppt_Wet] = wint;
  RGroup[rg]->ppt_slope[Ppt_Dry]  = dslope;
  RGroup[rg]->ppt_intcpt[Ppt_Dry] = dint;
}


static void _rgroup_add_disturbance( char name[], Int killyr, Int killfreq_startyr, 
                      RealF killfreq, Int extirp, RealF prop_killed, 
                      RealF prop_recovered, RealF grazing_frq, RealF prop_grazing, 
                      Int grazingfreq_startyr) {
/*======================================================*/
   GrpIndex rg;

   rg = RGroup_Name2Index(name);
   if (rg < 0) {
     LogError(logfp, LOGFATAL, "%s: Mismatched name (%s) for disturbance",
             MyFileName, name);
   }

  RGroup[rg]->killyr        = killyr;
  RGroup[rg]->killfreq_startyr = killfreq_startyr;
  RGroup[rg]->killfreq      = killfreq;
  RGroup[rg]->extirp        = (IntS) extirp;
  RGroup[rg]->proportion_killed    = prop_killed;
  RGroup[rg]->proportion_recovered = prop_recovered;
  RGroup[rg]->grazingfrq           = grazing_frq;
  RGroup[rg]->proportion_grazing   = prop_grazing;
  RGroup[rg]->grazingfreq_startyr  = grazingfreq_startyr;

  RGroup[rg]->extirpated    = FALSE;
}

/**************************************************************/
static void _rgroup_addsucculent( char name[], RealF wslope, RealF wint,
                                  RealF dslope, RealF dint) {
/* This should only get called one time, ie, there should be*/
/* only one set of succulent parameters, even if there is more*/
/* than one succulent species to which they pertain.*/
   GrpIndex rg;

   rg = RGroup_Name2Index(name);
   if (rg < 0) {
     LogError(logfp, LOGFATAL, "%s: Mismatched name (%s) for succulents",
             MyFileName, name);
   }
   RGroup[rg]->succulent = TRUE;
   Succulent->growth[Slope]  = wslope;
   Succulent->growth[Intcpt] = wint;
   Succulent->mort[Slope]  = dslope;
   Succulent->mort[Intcpt] = dint;
}

/** 
 * \brief Read the species-specific inputs and initialize the \ref Species
 *        array.
 * 
 * This function will read the species.in files and allocate enough memory
 * for the number of species requested. It then adds the species to the
 * \ref Species array and initializes all variables in the corresponding 
 * \ref SpeciesType struct.
 * 
 * \sideeffect Memory will be allocated for multiple \ref SpeciesType structs
 *             and the memory will be populated according to the input
 *             parameters.
 * 
 * \ingroup SPECIES_PRIVATE
 */
static void _species_init( void) {
   FILE *f;
   Bool readspp = TRUE, sppok = TRUE;

   /* temp vars to hold the group info*/
   char *name;
   size_t len;
   Int x;      /* temp val */
   GrpIndex rg;
   SppIndex sp;
   IntS age,
       slow,
       dist,
       eind,
       vegi,
       temp,
       turnon,
       turnondispersal,
       viable,
       pseed;
   RealF irate, ratep, estab, minb, maxb, cohort, xdecay,
         p1, p2, p3, p4, HMAX, PMD, HSlope;
   float var;
   char clonal[5];

   MyFileName = Parm_name( F_Species);
   f = OpenFile(MyFileName, "r");

   name = (char *)Mem_Calloc(SuperGlobals.max_speciesnamelen + 1, sizeof(char), "_species_init");
   Species = (SpeciesType **)Mem_Calloc(MAX_SPECIES, sizeof(SpeciesType *), "_species_init");

   while( readspp) {
      if( ! GetALine(f, inbuf )) {sppok=FALSE;break;}
      if ( ! isnull( strstr(inbuf,"[end]")) ) {
         readspp=FALSE;
         continue;
      }

		x = sscanf(inbuf,
				"%s %hd %hd %hd %hd %hd %hd %hd %hd %s %f %f %f %f %f %f",
				name, &rg, &turnon, &age, &slow, &dist, &eind, &vegi, &temp, clonal, &irate, &ratep,  &estab,
				&minb, &maxb, &cohort);
      if (x != 16) {
        LogError(logfp, LOGFATAL, "%s: Wrong number of columns in species",
                MyFileName);
      }

      sp = species_New();
      
      len = strlen(name);
      
      _setNameLen(Species[sp]->name, name, len);
      Species[sp]->sp_num  = (SppIndex) sp;
      Species[sp]->res_grp = (GrpIndex) rg-1; /* file gives natural # */
      Species[sp]->max_age = (IntS) age;
      Species[sp]->intrin_rate = irate;
      Species[sp]->max_rate = irate * ratep;
      Species[sp]->max_slow = (IntS) slow;
      switch ( dist ) {
        case 1:
          Species[sp]->disturbclass = VerySensitive;   break;
        case 2:
          Species[sp]->disturbclass = Sensitive;       break;
        case 3:
          Species[sp]->disturbclass = Insensitive;     break;
        case 4:
          Species[sp]->disturbclass = VeryInsensitive; break;
        default:
          LogError(logfp, LOGFATAL, "%s: Incorrect disturbance class found",
                  MyFileName);
      }

      Species[sp]->seedling_estab_prob = estab;
      Species[sp]->seedling_estab_prob_old = estab;
      Species[sp]->max_seed_estab = (IntS) eind;
      Species[sp]->seedling_biomass = minb;
      Species[sp]->mature_biomass = maxb;
      Species[sp]->tempclass = (temp==1) ? WarmSeason
                             : (temp==2) ? CoolSeason
                             : NoSeason;
      Species[sp]->relseedlingsize = minb / maxb;
      Species[sp]->isclonal = (Bool)(clonal[0]=='y' || clonal[0]=='Y');
      Species[sp]->max_vegunits = (Species[sp]->isclonal)
                                ? vegi
                                : 0;
      Species[sp]->use_me = (RGroup[rg-1]->use_me) ? itob(turnon) : FALSE ;
      Species[sp]->received_prob = 0;
      Species[sp]->cohort_surv = cohort;
      // input of `pseed` is in units of [# / m2]; convert to units of [# / plot]
      Species[sp]->pseed = pseed * Globals->plotsize;
/*      Species[sp]->ann_mort_prob = (age > 0)
                                         ? -log(cohort)/age
                                         : 0.0;
         */
    }/* end while*/

   if (!sppok) {
      LogError(logfp, LOGFATAL, "%s: Incorrect/incomplete input",
              MyFileName);
   }



   /* ------------------------------------------------- */
   /* ------ read the additional annuals' establishment parameters */
   sppok = readspp = TRUE;
   while( readspp) {
     if( ! GetALine(f, inbuf )) {sppok=FALSE; break;}
     if ( ! isnull( strstr(inbuf,"[end]")) ) {
        readspp=FALSE;
        continue;
     }

     x=sscanf( inbuf, "%s %hd %f %hd %f",
               name, &viable, &xdecay, &pseed, &var);
     if (x < 5) {
       LogError(logfp, LOGFATAL, "%s: Too few columns in annual estab parms",
               MyFileName);
     }

     sp = Species_Name2Index(name);
     if (sp < 0) {
       LogError(logfp, LOGFATAL, "%s: Mismatched name (%s) for annual estab parms",
               MyFileName, name);
     }

     Species[sp]->viable_yrs = viable;
     Species[sp]->exp_decay  = xdecay;
     Species[sp]->seedprod = (IntUS *) Mem_Calloc( viable, sizeof(IntUS), "species_init()");
     Species[sp]->var = var;
     Species[sp]->pseed = pseed / Globals->plotsize;

     /*Calculate alpha and beta for each species based on mean (pestab) and variance (var)*/
     Species[sp]->alpha = ((pow(Species[sp]->seedling_estab_prob, 2)
                            - pow(Species[sp]->seedling_estab_prob, 3))
                           / Species[sp]->var)
                          - Species[sp]->seedling_estab_prob;
     Species[sp]->beta = (Species[sp]->alpha / Species[sp]->seedling_estab_prob)
                          - Species[sp]->alpha;

     /*If the following two conditions are met, the beta distribution is bimodal or nearly bimodal,
      * which is not the desired outcome. The variance (s->var) and mean (pestab) should be adjusted
      * to obtain an unimodal beta-distribution with density > 0 */
     if (Species[sp]->alpha < 1) {
        LogError(logfp, LOGWARN, "Species %s, alpha less than 1: %f \n", 
                 Species[sp]->name, Species[sp]->alpha);
     }
     if (Species[sp]->beta < 1) {
        LogError(logfp, LOGWARN, "Species %s, beta less than 1: %f \n", 
                 Species[sp]->name, Species[sp]->beta);
     }
   } /* end while readspp*/

   if (!sppok) {
      LogError(logfp, LOGFATAL, "%s: Incorrect/incomplete input in annual estab parms",
              MyFileName);
   }


   /* ------------------------------------------------- */
   /* ------ read the vegetative propagation parameters */
   sppok = readspp = TRUE;
   while( readspp) {
      if( ! GetALine(f, inbuf )) {sppok=FALSE; break;}
      if ( ! isnull( strstr(inbuf,"[end]")) ) {
         readspp=FALSE;
         continue;
      }

      x=sscanf( inbuf, "%s %f %f %f %f",
                name, &p1, &p2, &p3, &p4);
      if (x < 5) {
        LogError(logfp, LOGFATAL, "%s: Too few columns in species probs",
                MyFileName);
      }

      sp = Species_Name2Index(name);
      if (sp < 0) {
        LogError(logfp, LOGFATAL, "%s: Mismatched name (%s) for species probs",
                MyFileName, name);
      }

      Species[sp]->prob_veggrow[NoResources] = p1;
      Species[sp]->prob_veggrow[Slow]        = p2;
      Species[sp]->prob_veggrow[Intrinsic]   = p3;
      Species[sp]->prob_veggrow[Disturbance] = p4;
    } /* end while readspp*/
   if (!sppok) {
      LogError(logfp, LOGFATAL, "%s: Incorrect/incomplete input in probs",
              MyFileName);
   }

  /* ------------------------------------------------- */
  /* ------- read the seed dispersal parameters ------ */
  sppok = readspp = TRUE;
  while( readspp ) {
    if( !GetALine(f, inbuf) ) {
      sppok=FALSE; break;
    }
    if( !isnull(strstr(inbuf,"[end]")) ) {
      readspp=FALSE;
      continue;
    }

    x = sscanf( inbuf, "%s %hd %f %f %f %f",
                name, &turnondispersal, &p1, &HMAX, &PMD, &HSlope); 
    if(x < 6) {
      LogError(logfp, LOGFATAL, "%s: Too few columns in species seed dispersal inputs", MyFileName);
    }

    sp = Species_Name2Index(name);
    if(sp < 0){
      LogError(logfp, LOGFATAL, "%s: Mismatched name (%s) for species seed dispersal inputs", MyFileName, name);
    }

    Species[sp]->use_dispersal = itob(turnondispersal);
    Species[sp]->minReproductiveSize = p1;
    Species[sp]->maxHeight = HMAX;
    Species[sp]->maxDispersalProbability = PMD;
    Species[sp]->heightSlope = HSlope;
  }
  if(!sppok) {
	  LogError(logfp, LOGFATAL, "%s: Incorrect/incomplete input in species seed dispersal input", MyFileName);
  }
   
   Mem_Free(name);

   CloseFile(&f);
}


/*static void _recover_names(void) {
	int i, last = NFILES - 1;

	last--; // have to save sxw.in for later //

	for (i = 0; i <= last; i++) {
		Mem_Free(_files[i]);
	}

}*/

/**************************************************************/
void parm_free_memory( void ) {
	//function to free memory allocated in this module
	int i;
	for(i = 0; i < NFILES; i++)
		Mem_Free(_files[i]);
}



#ifdef DEBUG_MEM
#include "sw_src/myMemory.h"
/**************************************************************/
void Parm_SetMemoryRefs( void) {
/*======================================================*/
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
  ST_FileIndex i;
  GrpIndex rg;
  SppIndex sp;

  for(i=F_First; i<= F_MortAvg; i++)
    NoteMemoryRef(_files[i]);

  ForEachGroup(rg)
    NoteMemoryRef( RGroup[rg]->kills);
  ForEachSpecies(sp)
    NoteMemoryRef( Species[sp]->kills);

  NoteMemoryRef(_files[F_SXW]);
}

#endif
