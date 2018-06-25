/********************************************************/
/********************************************************/
/*  Source file: output.c
 *  Type: module
 *  Application: STEPPE - plant community dynamics simulator
 *  Purpose: All functions pertaining to producing output
 *           should go here. */
/*  History */
/*     (6/15/2000) -- INITIAL CODING - cwb */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "filefuncs.h"


/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
#ifdef STEPWAT
  #include "SW_Model.h"
  extern SW_MODEL SW_Model;
  extern Bool UseSoilwat;
#endif

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
  void output_Bmass_Yearly( Int year );
  void output_Mort_Yearly( void );



/**************************************************************/
/**************************************************************/
void output_Bmass_Yearly( Int year ) {
/*======================================================*/
/* note that year is only printed, not used as index, so
 * we don't need to decrement to make base0
 */
  char fields[MAX_OUTFIELDS][MAX_FIELDLEN+1];
  GrpIndex rg;
  SppIndex sp;
  int i, fc=0;
  char s[MAX_FIELDLEN];
  char filename[FILENAME_MAX];


  if (!BmassFlags.yearly) return;

  if (BmassFlags.yr) {
      if (UseSoilwat)
        sprintf(fields[fc++], "%d", SW_Model.year);
      else
        sprintf(fields[fc++], "%d", year);
  }

  if (BmassFlags.dist) {
    switch (Plot.disturbance) {
      case NoDisturb: strcpy(s, "None"); break;
      case FecalPat: strcpy(s, "Pat"); break;
      case AntMound: strcpy(s, "Mound"); break;
      case Burrow:   strcpy(s, "Burrow"); break;
    default: strcpy(s, "Unknown!"); break;
    }
    sprintf(fields[fc++], "%s", s);
  }

  if (BmassFlags.ppt)
    sprintf(fields[fc++], "%d", Env.ppt);

  if (BmassFlags.pclass) {
    switch (Env.wet_dry) {
      case Ppt_Norm: strcpy(s, "Normal"); break;
      case Ppt_Wet: strcpy(s, "Wet"); break;
      case Ppt_Dry: strcpy(s, "Dry"); break;
    default: strcpy(s, "Unknown!"); break;
    }
    sprintf(fields[fc++], "%s", s);
  }

  if (BmassFlags.tmp)
    sprintf(fields[fc++], "%0.1f", Env.temp);

  if (BmassFlags.grpb) {
    ForEachGroup(rg) {
      sprintf(fields[fc++], "%f", RGroup_GetBiomass(rg));
      if (BmassFlags.size)
        sprintf(fields[fc++],"%f", RGroup[rg]->relsize);
      if (BmassFlags.pr)
        sprintf(fields[fc++],"%f", RGroup[rg]->pr);
    }
  }

  if (BmassFlags.sppb) {
    ForEachSpecies(sp) {
      sprintf(fields[fc++], "%f", Species_GetBiomass(sp));
      if (BmassFlags.indv)
        sprintf(fields[fc++],"%d", Species[sp]->est_count);
    }
  }
  sprintf(filename, "%s%0*d.csv", Parm_name(F_BMassPre),
                                 Globals.bmass.suffixwidth,
                                 Globals.currIter);
  Globals.bmass.fp_year = OpenFile(filename, "a");

  /* Write data line to already opened file */
  for (i=0; i< fc-1; i++) {
	  fprintf(Globals.bmass.fp_year,"%s%c", fields[i], BmassFlags.sep);
  }

  if (i) fprintf(Globals.bmass.fp_year,"%s\n", fields[i]);
  fflush(Globals.bmass.fp_year);
  CloseFile(&Globals.bmass.fp_year);
}


/***********************************************************/
void output_Mort_Yearly( void ) {
/*======================================================*/

	IntS age, rg, sp;
	char filename[FILENAME_MAX];

	sprintf(filename, "%s%0*d.csv", Parm_name(F_MortPre), Globals.mort.suffixwidth, Globals.currIter);
	Globals.mort.fp_year = OpenFile(filename, "a");
	FILE *f = Globals.mort.fp_year;

	if (!MortFlags.yearly)
		return;

	/* Note: Header line already printed */

	/* Print a line of establishments */
	fprintf(f, "(Estabs)");
	if (MortFlags.group) {
		ForEachGroup(rg)
			fprintf(Globals.mort.fp_year, "%c%d", MortFlags.sep, RGroup[rg]->estabs);
	}
	if (MortFlags.species) {
		ForEachSpecies(sp)
			fprintf(f, "%c%d", MortFlags.sep, Species[sp]->estabs);
	}
	fprintf(f, "\n");


  /* now print the kill data */
	for (age = 0; age < Globals.Max_Age; age++) {
		fprintf(f, "%d", age + 1);
		if (MortFlags.group) {
			ForEachGroup(rg)
			{
				if (age < GrpMaxAge(rg) && RGroup[rg]->use_me)
					fprintf(f, "%c%d", MortFlags.sep, RGroup[rg]->kills[age]);
				else
					fprintf(f, "%c", MortFlags.sep);
			}
		}
		if (MortFlags.species) {
			ForEachSpecies(sp)
			{
				/*adding fix for bus error while on/off some species:
				 *  Reason proper species use-me boolean values were not use in check so added the same
				 *  Modify By: Ashish
				 */
				if (age < SppMaxAge(sp) && Species[sp]->use_me && RGroup[Species[sp]->res_grp]->use_me)
					fprintf(f, "%c%d", MortFlags.sep, Species[sp]->kills[age]);
				else
					fprintf(f, "%c", MortFlags.sep);
			}
		}
		fprintf(f, "\n");
	}

	CloseFile(&Globals.mort.fp_year);
}
