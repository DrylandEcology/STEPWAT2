/**
 * \file ST_output.c
 * \brief Outputs mortality or biomass on a yearly time step. 
 * 
 * This file differs from \ref ST_stats.c because it does not accumulate values.
 * It simply prints the values to the given year.
 * 
 * \ingroup OUTPUT
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <string.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"


/******** Modular External Function Declarations ***********/
/* -- truly global functions are declared in functions.h --*/
/***********************************************************/
#include "sw_src/SW_Model.h"
extern SW_MODEL SW_Model;

/*------------------------------------------------------*/
/* Modular functions only used on one or two specific   */
/* places; that is, they are not generally useful       */
/* (like C++ friend functions) but have to be declared. */
  void output_Bmass_Yearly( Int year );
  void output_Mort_Yearly( void );



/**
 * \brief Outputs the current year's values to the file denoted in 
 *        [Globals.bmass.fp_year](\ref Globals)
 * 
 * \param year is the year that these values are being printed. This is 1 indexed.
 * 
 * \ingroup OUTPUT
 */
void output_Bmass_Yearly( Int year ) {
  char **fields;
  GrpIndex rg;
  SppIndex sp;
  int i, fc = 0;
  char *s;
  char filename[FILENAME_MAX];
  
  if (!BmassFlags.yearly) return;

  fields = (char **)Mem_Calloc(MAX_OUTFIELDS, sizeof(char *), "output_Bmass_Yearly");
  s = (char *)Mem_Calloc(MAX_FIELDLEN + 1, sizeof(char), "output_Bmass_Yearly");
  
  for (i = 0; i < MAX_OUTFIELDS; i++) {
      fields[i] = (char *)Mem_Calloc(MAX_FIELDLEN + 1, sizeof(char), "output_Bmass_Yearly");
  }
  
  if(Globals->currYear == 1) // At year one we need a header.
  {
  /* --------- Begin setting up header ------- */

    if (BmassFlags.yr)
      strcpy(fields[fc++], "Year");
    if (BmassFlags.dist)
      strcpy(fields[fc++], "Disturbs");
    if (BmassFlags.ppt) {
      strcpy(fields[fc++], "PPT");
    }
    if (BmassFlags.pclass)
      strcpy(fields[fc++], "PPTClass");
    if (BmassFlags.tmp) {
      strcpy(fields[fc++], "Temp");
    }
    if (BmassFlags.grpb) {
      if(BmassFlags.wildfire){
        strcpy(fields[fc++], "Wildfire");
      }
      ForEachGroup(rg) 
      {
        strcpy(fields[fc++], RGroup[rg]->name);
        if (BmassFlags.size) {
          strcpy(fields[fc], RGroup[rg]->name);
          strcat(fields[fc++], "_RSize");
        }
        if (BmassFlags.pr) {
          strcpy(fields[fc], RGroup[rg]->name);
          strcat(fields[fc++],"_PR");
        }
        if(BmassFlags.prescribedfire){
          strcpy(fields[fc], RGroup[rg]->name);
          strcat(fields[fc++],"_PFire");
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
    sprintf(filename, "%s%0*d.csv", Parm_name(F_BMassPre),
                                 Globals->bmass.suffixwidth,
                                 Globals->currIter);
    Globals->bmass.fp_year = OpenFile(filename, "a");

    /* Write data line to already opened file */
    for (i=0; i< fc-1; i++) {
      fprintf(Globals->bmass.fp_year,"%s%c", fields[i], BmassFlags.sep);
    }

    if (i) fprintf(Globals->bmass.fp_year,"%s\n", fields[i]);
    fflush(Globals->bmass.fp_year);
    CloseFile(&Globals->bmass.fp_year);

    fc = 0; //reset fc for first line of data.
  }
  /* ------------- end setting up header -------------- */

  if (BmassFlags.yr) {
    sprintf(fields[fc++], "%d", SW_Model.year);
  }

  if (BmassFlags.dist) {
    switch (Plot->disturbance) {
      case NoDisturb: strcpy(s, "None"); break;
      case FecalPat: strcpy(s, "Pat"); break;
      case AntMound: strcpy(s, "Mound"); break;
      case Burrow:   strcpy(s, "Burrow"); break;
    default: strcpy(s, "Unknown!"); break;
    }
    sprintf(fields[fc++], "%s", s);
  }

  if (BmassFlags.ppt)
    sprintf(fields[fc++], "%d", Env->ppt);

  if (BmassFlags.pclass) {
    switch (Env->wet_dry) {
      case Ppt_Norm: strcpy(s, "Normal"); break;
      case Ppt_Wet: strcpy(s, "Wet"); break;
      case Ppt_Dry: strcpy(s, "Dry"); break;
    default: strcpy(s, "Unknown!"); break;
    }
    sprintf(fields[fc++], "%s", s);
  }

  if (BmassFlags.tmp)
    sprintf(fields[fc++], "%0.1f", Env->temp);

  if (BmassFlags.grpb) {
    if (BmassFlags.wildfire)
        sprintf(fields[fc++],"%d", RGroup[0]->wildfire);
    ForEachGroup(rg) {
      sprintf(fields[fc++], "%f", RGroup_GetBiomass(rg));
      if (BmassFlags.size)
        sprintf(fields[fc++],"%f", getRGroupRelsize(rg));
      if (BmassFlags.pr)
        sprintf(fields[fc++],"%f", RGroup[rg]->pr);
      if (BmassFlags.prescribedfire)
        sprintf(fields[fc++],"%d", RGroup[rg]->prescribedfire);
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
                                 Globals->bmass.suffixwidth,
                                 Globals->currIter);
  Globals->bmass.fp_year = OpenFile(filename, "a");

  /* Write data line to already opened file */
  for (i=0; i< fc-1; i++) {
	  fprintf(Globals->bmass.fp_year,"%s%c", fields[i], BmassFlags.sep);
  }

  if (i) fprintf(Globals->bmass.fp_year,"%s\n", fields[i]);
  fflush(Globals->bmass.fp_year);
  CloseFile(&Globals->bmass.fp_year);
  
  for (i = 0; i < MAX_OUTFIELDS; i++) {
      Mem_Free(fields[i]);
  }
  
  Mem_Free(s);
  Mem_Free(fields);
}


/**
 * \brief Outputs the current year's values to the file specified in
 *        [Globals->mort.fp_year](\ref Globals)
 * 
 * Prints mortality values. These values are indexed by age at death.
 * 
 * \ingroup OUTPUT
 */
void output_Mort_Yearly( void ) {
	IntS age, rg, sp;
	char filename[FILENAME_MAX];

	sprintf(filename, "%s%0*d.csv", Parm_name(F_MortPre), Globals->mort.suffixwidth, Globals->currIter);
	Globals->mort.fp_year = OpenFile(filename, "a");
	FILE *f = Globals->mort.fp_year;

	if (!MortFlags.yearly)
		return;

    /* Print header line */
	fprintf(f, "Age");

	if (MortFlags.group) {
        ForEachGroup(rg) fprintf(f, "%c%s", MortFlags.sep, RGroup[rg]->name);
	}

	if (MortFlags.species) {
        ForEachSpecies(sp) fprintf(f, "%c%s", MortFlags.sep, Species[sp]->name);
	}

	/* Header line is now complete */
	fprintf(f, "\n");

	/* Print a line of establishments */
	fprintf(f, "(Estabs)");
	if (MortFlags.group) {
		ForEachGroup(rg)
			fprintf(Globals->mort.fp_year, "%c%d", MortFlags.sep, RGroup[rg]->estabs);
	}
	if (MortFlags.species) {
		ForEachSpecies(sp)
			fprintf(f, "%c%d", MortFlags.sep, Species[sp]->estabs);
	}
	fprintf(f, "\n");


  /* now print the kill data */
	for (age = 0; age < Globals->Max_Age; age++) {
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

	CloseFile(&Globals->mort.fp_year);
}
