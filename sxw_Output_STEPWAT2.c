/********************************************************/
/********************************************************/
/*  Source file: sxw_Output_STEPWAT2.c
  Type: module
  Application: SOILWAT - soilwater dynamics simulator
  Purpose: define `get_XXX` functions for STEPWAT2 interface
    see SW_Output_core.c and SW_Output.h

  History:
  2018 June 04 (drs) moved output formatter `get_XXX` functions from
    previous `SW_Output.c` to dedicated `sxw_Output_STEPWAT2.c`
 */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sw_src/generic.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sw_src/Times.h"

#include "sw_src/SW_Carbon.h"
#include "sw_src/SW_Defines.h"
#include "sw_src/SW_Files.h"
#include "sw_src/SW_Model.h"
#include "sw_src/SW_Site.h"
#include "sw_src/SW_SoilWater.h"
#include "sw_src/SW_Times.h"
#include "sw_src/SW_Weather.h"
#include "sw_src/SW_VegEstab.h"
#include "sw_src/SW_VegProd.h"

#include "sw_src/SW_Output.h"

#include "sxw.h"
#include "ST_globals.h"


/* =================================================== */
/*                  Global Variables                   */
/* --------------------------------------------------- */
extern SW_SITE SW_Site;
extern SW_SOILWAT SW_Soilwat;
extern SW_MODEL SW_Model;
extern SW_WEATHER SW_Weather;
extern SW_VEGPROD SW_VegProd;
extern SW_VEGESTAB SW_VegEstab;
extern SW_CARBON SW_Carbon;
extern SW_OUTPUT SW_Output[SW_OUTNKEYS]; // defined in `SW_Output_core.c`

extern SXW_t SXW; // structure to store values in and pass back to STEPPE
extern SXW_avg SXW_AVG;
//extern ModelType Globals; // defined in `ST_Main.c`, but externed by `ST_globals.h`
extern Bool isPartialSoilwatOutput; // defined in `SW_Output_core.c`
extern Bool storeAllIterations; // defined in `SW_Output_core.c`
extern char sw_outstr_agg[OUTSTRLEN]; // defined in `SW_Output_core.c`

extern char _Sep; // defined in `SW_Output_core.c`: output delimiter
extern TimeInt tOffset; // defined in `SW_Output_core.c`: 1 or 0 means we're writing previous or current period
extern Bool bFlush_output; // defined in `SW_Output_core.c`: process partial period ?
extern char sw_outstr[OUTSTRLEN]; // defined in `SW_Output_core.c`: string with formatted output which is to be written to output files



/* =================================================== */
/* =================================================== */
/*             Private Functions                       */
/* --------------------------------------------------- */

static void _create_filename_ST(char *str, char *flag, int iteration, char *filename);
static void _create_csv_file_ST(int iteration, OutPeriod pd);

/** Splits a filename such as `name.ext` into its two parts `name` and `ext`;
		appends `flag` and, if positive, `iteration` to `name` with `_` as
		separator; and returns the full name concatenated

		\return `name_flagiteration.ext`
*/
static void _create_filename_ST(char *str, char *flag, int iteration, char *filename) {
	char *basename;
	char *ext;
	char *fileDup = (char *)malloc(strlen(str) + 1);

	// Determine basename and file extension
	strcpy(fileDup, str); // copy file name to new variable
	basename = strtok(fileDup, ".");
	ext = strtok(NULL, ".");
	free(fileDup);

	// Put new file together
	if (iteration > 0) {
		sprintf(filename, "%s_%s%d.%s", basename, flag, iteration, ext);
	} else {
		sprintf(filename, "%s_%s.%s", basename, flag, ext);
	}
}


/**
  \fn void _create_csv_file_ST(int iteration, OutPeriod pd)

  Creates `csv` output files for specified time step depending on
  `-o` and `-i` flags, for `STEPWAT2`.

  If `-i` flag is used, then this function creates a file for each `iteration`
  with the file name containing the value of `iteration`.

  If `-o` flag is used, then this function creates only one set of output files.

  \param iteration. Current iteration value that is used for the file name
    if -i flag used in STEPWAT2. Set to a negative value otherwise.
  \param pd. The output time step.
*/
/***********************************************************/
static void _create_csv_file_ST(int iteration, OutPeriod pd)
{
	char filename[FILENAME_MAX];

	if (iteration <= 0)
	{ // STEPWAT2: aggregated values over all iterations or SOILWAT2-standalone
		if (SW_OutFiles.make_regular) {
			// PROGRAMMER Note: `eOutputDaily + pd` is not very elegant and assumes
			// a specific order of `SW_FileIndex` --> fix and create something that
			// allows subsetting such as `eOutputFile[pd]` or append time period to
			// a basename, etc.
			_create_filename_ST(SW_F_name(eOutputDaily + pd), "agg", 0, filename);
			SW_OutFiles.fp_reg_agg[pd] = OpenFile(filename, "w");
		}

		if (SW_OutFiles.make_soil) {
			_create_filename_ST(SW_F_name(eOutputDaily_soil + pd), "agg", 0, filename);
			SW_OutFiles.fp_soil_agg[pd] = OpenFile(filename, "w");
		}

	} else
	{ // STEPWAT2: storing values for every iteration
		if (iteration > 1) {
			// close files from previous iteration
			if (SW_OutFiles.make_regular) {
				CloseFile(&SW_OutFiles.fp_reg[pd]);
			}
			if (SW_OutFiles.make_soil) {
				CloseFile(&SW_OutFiles.fp_soil[pd]);
			}
		}

		if (SW_OutFiles.make_regular) {
			_create_filename_ST(SW_F_name(eOutputDaily + pd), "rep", iteration, filename);
			SW_OutFiles.fp_reg[pd] = OpenFile(filename, "w");
		}

		if (SW_OutFiles.make_soil) {
			_create_filename_ST(SW_F_name(eOutputDaily_soil + pd), "rep", iteration, filename);
			SW_OutFiles.fp_soil[pd] = OpenFile(filename, "w");
		}
	}
	#endif
}



/* =================================================== */
/* =================================================== */
/*             Function Definitions                    */
/*             (declared in SW_Output.h)               */
/* --------------------------------------------------- */


void SW_OUT_create_summary_files(void) {
	OutPeriod p;

		if (!isPartialSoilwatOutput)
		{
			ForEachOutPeriod(p) {
				if (use_OutPeriod[p]) {
					_create_csv_file_ST(-1, p);

					write_headers_to_csv(p, SW_OutFiles.fp_reg_agg[p],
						SW_OutFiles.fp_soil_agg[p], swTRUE);
				}
			}
		}
}

void SW_OUT_create_iteration_files(void) {
	OutPeriod p;

	if (storeAllIterations) {
		ForEachOutPeriod(p) {
			if (use_OutPeriod[p]) {
				_create_csv_file_ST(Globals.currIter + 1, p);

				write_headers_to_csv(p, SW_OutFiles.fp_reg[p],
					SW_OutFiles.fp_soil[p], swFALSE);
			}
		}
	}
}



void get_co2effects(OutPeriod pd) {
	SW_VEGPROD *v = &SW_VegProd;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;

	RealD biomass_total = SW_MISSING, biolive_total = SW_MISSING;
	RealD biomass_grass = SW_MISSING, biomass_shrub = SW_MISSING,
		biomass_tree = SW_MISSING, biomass_forb = SW_MISSING;
	RealD biolive_grass = SW_MISSING, biolive_shrub = SW_MISSING,
		biolive_tree = SW_MISSING, biolive_forb = SW_MISSING;

	// Grab the multipliers that were just used
	// No averaging or summing required
	RealD bio_mult_grass = v->veg[SW_GRASS].co2_multipliers[BIO_INDEX][SW_Model.simyear];
	RealD bio_mult_shrub = v->veg[SW_SHRUB].co2_multipliers[BIO_INDEX][SW_Model.simyear];
	RealD bio_mult_tree = v->veg[SW_TREES].co2_multipliers[BIO_INDEX][SW_Model.simyear];
	RealD bio_mult_forb = v->veg[SW_FORBS].co2_multipliers[BIO_INDEX][SW_Model.simyear];
	RealD wue_mult_grass = v->veg[SW_GRASS].co2_multipliers[WUE_INDEX][SW_Model.simyear];
	RealD wue_mult_shrub = v->veg[SW_SHRUB].co2_multipliers[WUE_INDEX][SW_Model.simyear];
	RealD wue_mult_tree = v->veg[SW_TREES].co2_multipliers[WUE_INDEX][SW_Model.simyear];
	RealD wue_mult_forb = v->veg[SW_FORBS].co2_multipliers[WUE_INDEX][SW_Model.simyear];

	if ((isPartialSoilwatOutput == FALSE &&
		Globals.currIter == Globals.runModelIterations) || storeAllIterations) {
		get_outstrleader2(pd);
	}

	switch(pd) {
		case eSW_Day:
			biomass_grass = v->dysum.veg[SW_GRASS].biomass;
			biomass_shrub = v->dysum.veg[SW_SHRUB].biomass;
			biomass_tree = v->dysum.veg[SW_TREES].biomass;
			biomass_forb = v->dysum.veg[SW_FORBS].biomass;
			biolive_grass = v->dysum.veg[SW_GRASS].biolive;
			biolive_shrub = v->dysum.veg[SW_SHRUB].biolive;
			biolive_tree = v->dysum.veg[SW_TREES].biolive;
			biolive_forb = v->dysum.veg[SW_FORBS].biolive;
			biomass_total = biomass_grass + biomass_shrub + biomass_tree + biomass_forb;
			biolive_total = biolive_grass + biolive_shrub + biolive_tree + biolive_forb;
			break;

		case eSW_Week:
			biomass_grass = v->wkavg.veg[SW_GRASS].biomass;
			biomass_shrub = v->wkavg.veg[SW_SHRUB].biomass;
			biomass_tree = v->wkavg.veg[SW_TREES].biomass;
			biomass_forb = v->wkavg.veg[SW_FORBS].biomass;
			biolive_grass = v->wkavg.veg[SW_GRASS].biolive;
			biolive_shrub = v->wkavg.veg[SW_SHRUB].biolive;
			biolive_tree = v->wkavg.veg[SW_TREES].biolive;
			biolive_forb = v->wkavg.veg[SW_FORBS].biolive;
			biomass_total = biomass_grass + biomass_shrub + biomass_tree + biomass_forb;
			biolive_total = biolive_grass + biolive_shrub + biolive_tree + biolive_forb;
			break;

		case eSW_Month:
			biomass_grass = v->moavg.veg[SW_GRASS].biomass;
			biomass_shrub = v->moavg.veg[SW_SHRUB].biomass;
			biomass_tree = v->moavg.veg[SW_TREES].biomass;
			biomass_forb = v->moavg.veg[SW_FORBS].biomass;
			biolive_grass = v->moavg.veg[SW_GRASS].biolive;
			biolive_shrub = v->moavg.veg[SW_SHRUB].biolive;
			biolive_tree = v->moavg.veg[SW_TREES].biolive;
			biolive_forb = v->moavg.veg[SW_FORBS].biolive;
			biomass_total = biomass_grass + biomass_shrub + biomass_tree + biomass_forb;
			biolive_total = biolive_grass + biolive_shrub + biolive_tree + biolive_forb;
			break;

		case eSW_Year:
			biomass_grass = v->yravg.veg[SW_GRASS].biomass;
			biomass_shrub = v->yravg.veg[SW_SHRUB].biomass;
			biomass_tree = v->yravg.veg[SW_TREES].biomass;
			biomass_forb = v->yravg.veg[SW_FORBS].biomass;
			biolive_grass = v->yravg.veg[SW_GRASS].biolive;
			biolive_shrub = v->yravg.veg[SW_SHRUB].biolive;
			biolive_tree = v->yravg.veg[SW_TREES].biolive;
			biolive_forb = v->yravg.veg[SW_FORBS].biolive;
			biomass_total = biomass_grass + biomass_shrub + biomass_tree + biomass_forb;
			biolive_total = biolive_grass + biolive_shrub + biolive_tree + biolive_forb;
			break;
	}

		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				break;
			case eSW_Week:
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				break;
			case eSW_Year:
				p = Globals.currYear-1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd

			float old_biomass_grass = SXW_AVG.biomass_grass_avg[ind0];
			float old_biomass_shrub = SXW_AVG.biomass_shrub_avg[ind0];
			float old_biomass_tree = SXW_AVG.biomass_tree_avg[ind0];
			float old_biomass_forb = SXW_AVG.biomass_forb_avg[ind0];
			float old_biomass_total = SXW_AVG.biomass_total_avg[ind0];

			float old_biolive_grass = SXW_AVG.biolive_grass_avg[ind0];
			float old_biolive_shrub = SXW_AVG.biolive_shrub_avg[ind0];
			float old_biolive_tree = SXW_AVG.biolive_tree_avg[ind0];
			float old_biolive_forb = SXW_AVG.biolive_forb_avg[ind0];
			float old_biolive_total = SXW_AVG.biolive_total_avg[ind0];

			float old_bio_mult_grass = SXW_AVG.bio_mult_grass_avg[ind0];
			float old_bio_mult_shrub = SXW_AVG.bio_mult_shrub_avg[ind0];
			float old_bio_mult_tree = SXW_AVG.bio_mult_tree_avg[ind0];
			float old_bio_mult_forb = SXW_AVG.bio_mult_forb_avg[ind0];

			float old_wue_mult_grass = SXW_AVG.wue_mult_grass_avg[ind0];
			float old_wue_mult_shrub = SXW_AVG.wue_mult_shrub_avg[ind0];
			float old_wue_mult_tree = SXW_AVG.wue_mult_tree_avg[ind0];
			float old_wue_mult_forb = SXW_AVG.wue_mult_forb_avg[ind0];


			SXW_AVG.biomass_grass_avg[ind0] = get_running_avg(SXW_AVG.biomass_grass_avg[ind0], biomass_grass);
			SXW_AVG.biomass_shrub_avg[ind0] = get_running_avg(SXW_AVG.biomass_shrub_avg[ind0], biomass_shrub);
			SXW_AVG.biomass_tree_avg[ind0] = get_running_avg(SXW_AVG.biomass_tree_avg[ind0], biomass_tree);
			SXW_AVG.biomass_forb_avg[ind0] = get_running_avg(SXW_AVG.biomass_forb_avg[ind0], biomass_forb);
			SXW_AVG.biomass_total_avg[ind0] = get_running_avg(SXW_AVG.biomass_total_avg[ind0], biomass_total);

			SXW_AVG.biolive_grass_avg[ind0] = get_running_avg(SXW_AVG.biolive_grass_avg[ind0], biolive_grass);
			SXW_AVG.biolive_shrub_avg[ind0] = get_running_avg(SXW_AVG.biolive_shrub_avg[ind0], biolive_shrub);
			SXW_AVG.biolive_tree_avg[ind0] = get_running_avg(SXW_AVG.biolive_tree_avg[ind0], biolive_tree);
			SXW_AVG.biolive_forb_avg[ind0] = get_running_avg(SXW_AVG.biolive_forb_avg[ind0], biolive_forb);
			SXW_AVG.biolive_total_avg[ind0] = get_running_avg(SXW_AVG.biolive_total_avg[ind0], biolive_total);

			SXW_AVG.bio_mult_grass_avg[ind0] = get_running_avg(SXW_AVG.bio_mult_grass_avg[ind0], bio_mult_grass);
			SXW_AVG.bio_mult_shrub_avg[ind0] = get_running_avg(SXW_AVG.bio_mult_shrub_avg[ind0], bio_mult_shrub);
			SXW_AVG.bio_mult_tree_avg[ind0] = get_running_avg(SXW_AVG.bio_mult_tree_avg[ind0], bio_mult_tree);
			SXW_AVG.bio_mult_forb_avg[ind0] = get_running_avg(SXW_AVG.bio_mult_forb_avg[ind0], bio_mult_forb);

			SXW_AVG.wue_mult_grass_avg[ind0] = get_running_avg(SXW_AVG.wue_mult_grass_avg[ind0], wue_mult_grass);
			SXW_AVG.wue_mult_shrub_avg[ind0] = get_running_avg(SXW_AVG.wue_mult_shrub_avg[ind0], wue_mult_shrub);
			SXW_AVG.wue_mult_tree_avg[ind0] = get_running_avg(SXW_AVG.wue_mult_tree_avg[ind0], wue_mult_tree);
			SXW_AVG.wue_mult_forb_avg[ind0] = get_running_avg(SXW_AVG.wue_mult_forb_avg[ind0], wue_mult_forb);

			// ---------------------------

			SXW_AVG.biomass_grass_avg[ind1] += get_running_sqr(old_biomass_grass, biomass_grass, SXW_AVG.biomass_grass_avg[ind0]);
			SXW_AVG.biomass_shrub_avg[ind1] += get_running_sqr(old_biomass_shrub, biomass_shrub, SXW_AVG.biomass_shrub_avg[ind0]);
			SXW_AVG.biomass_tree_avg[ind1] += get_running_sqr(old_biomass_tree, biomass_tree, SXW_AVG.biomass_tree_avg[ind0]);
			SXW_AVG.biomass_forb_avg[ind1] += get_running_sqr(old_biomass_forb, biomass_forb, SXW_AVG.biomass_forb_avg[ind0]);
			SXW_AVG.biomass_total_avg[ind1] += get_running_sqr(old_biomass_total, biomass_total, SXW_AVG.biomass_total_avg[ind0]);

			SXW_AVG.biolive_grass_avg[ind1] += get_running_sqr(old_biolive_grass, biolive_grass, SXW_AVG.biolive_grass_avg[ind0]);
			SXW_AVG.biolive_shrub_avg[ind1] += get_running_sqr(old_biolive_shrub, biolive_shrub, SXW_AVG.biolive_shrub_avg[ind0]);
			SXW_AVG.biolive_tree_avg[ind1] += get_running_sqr(old_biolive_tree, biolive_tree, SXW_AVG.biolive_tree_avg[ind0]);
			SXW_AVG.biolive_forb_avg[ind1] += get_running_sqr(old_biolive_forb, biolive_forb, SXW_AVG.biolive_forb_avg[ind0]);
			SXW_AVG.biolive_total_avg[ind1] += get_running_sqr(old_biolive_total, biolive_total, SXW_AVG.biolive_total_avg[ind0]);

			SXW_AVG.bio_mult_grass_avg[ind1] += get_running_sqr(old_bio_mult_grass, bio_mult_grass, SXW_AVG.bio_mult_grass_avg[ind0]);
			SXW_AVG.bio_mult_shrub_avg[ind1] += get_running_sqr(old_bio_mult_shrub, bio_mult_shrub, SXW_AVG.bio_mult_shrub_avg[ind0]);
			SXW_AVG.bio_mult_tree_avg[ind1] += get_running_sqr(old_bio_mult_tree, bio_mult_tree, SXW_AVG.bio_mult_tree_avg[ind0]);
			SXW_AVG.bio_mult_forb_avg[ind1] += get_running_sqr(old_bio_mult_forb, bio_mult_forb, SXW_AVG.bio_mult_forb_avg[ind0]);

			SXW_AVG.wue_mult_grass_avg[ind1] += get_running_sqr(old_wue_mult_grass, wue_mult_grass, SXW_AVG.wue_mult_grass_avg[ind0]);
			SXW_AVG.wue_mult_shrub_avg[ind1] += get_running_sqr(old_wue_mult_shrub, wue_mult_shrub, SXW_AVG.wue_mult_shrub_avg[ind0]);
			SXW_AVG.wue_mult_tree_avg[ind1] += get_running_sqr(old_wue_mult_tree, wue_mult_tree, SXW_AVG.wue_mult_tree_avg[ind0]);
			SXW_AVG.wue_mult_forb_avg[ind1] += get_running_sqr(old_wue_mult_forb, wue_mult_forb, SXW_AVG.wue_mult_forb_avg[ind0]);


			if(Globals.currIter == Globals.runModelIterations){
				float std_biomass_grass = sqrt(SXW_AVG.biomass_grass_avg[ind1] / Globals.currIter);
				float std_biomass_shrub = sqrt(SXW_AVG.biomass_shrub_avg[ind1] / Globals.currIter);
				float std_biomass_tree = sqrt(SXW_AVG.biomass_tree_avg[ind1] / Globals.currIter);
				float std_biomass_forb = sqrt(SXW_AVG.biomass_forb_avg[ind1] / Globals.currIter);
				float std_biomass_total = sqrt(SXW_AVG.biomass_total_avg[ind1] / Globals.currIter);

				float std_biolive_grass = sqrt(SXW_AVG.biolive_grass_avg[ind1] / Globals.currIter);
				float std_biolive_shrub = sqrt(SXW_AVG.biolive_shrub_avg[ind1] / Globals.currIter);
				float std_biolive_tree = sqrt(SXW_AVG.biolive_tree_avg[ind1] / Globals.currIter);
				float std_biolive_forb = sqrt(SXW_AVG.biolive_forb_avg[ind1] / Globals.currIter);
				float std_biolive_total = sqrt(SXW_AVG.biolive_total_avg[ind1] / Globals.currIter);

				float std_bio_mult_grass = sqrt(SXW_AVG.bio_mult_grass_avg[ind1] / Globals.currIter);
				float std_bio_mult_shrub = sqrt(SXW_AVG.bio_mult_shrub_avg[ind1] / Globals.currIter);
				float std_bio_mult_tree = sqrt(SXW_AVG.bio_mult_tree_avg[ind1] / Globals.currIter);
				float std_bio_mult_forb = sqrt(SXW_AVG.bio_mult_forb_avg[ind1] / Globals.currIter);

				float std_wue_mult_grass = sqrt(SXW_AVG.wue_mult_grass_avg[ind1] / Globals.currIter);
				float std_wue_mult_shrub = sqrt(SXW_AVG.wue_mult_shrub_avg[ind1] / Globals.currIter);
				float std_wue_mult_tree = sqrt(SXW_AVG.wue_mult_tree_avg[ind1] / Globals.currIter);
				float std_wue_mult_forb = sqrt(SXW_AVG.wue_mult_forb_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f",
				_Sep, SXW_AVG.biomass_grass_avg[ind0], _Sep, std_biomass_grass,
				_Sep, SXW_AVG.biomass_shrub_avg[ind0], _Sep, std_biomass_shrub,
				_Sep, SXW_AVG.biomass_tree_avg[ind0], _Sep, std_biomass_tree,
				_Sep, SXW_AVG.biomass_forb_avg[ind0], _Sep, std_biomass_forb,
				_Sep, SXW_AVG.biomass_total_avg[ind0], _Sep, std_biomass_total,
				_Sep, SXW_AVG.biolive_grass_avg[ind0], _Sep, std_biolive_grass,
				_Sep, SXW_AVG.biolive_shrub_avg[ind0], _Sep, std_biolive_shrub,
				_Sep, SXW_AVG.biolive_tree_avg[ind0], _Sep, std_biolive_tree,
				_Sep, SXW_AVG.biolive_forb_avg[ind0], _Sep, std_biolive_forb,
				_Sep, SXW_AVG.biolive_total_avg[ind0], _Sep, std_biolive_total,
				_Sep, SXW_AVG.bio_mult_grass_avg[ind0], _Sep, std_bio_mult_grass,
				_Sep, SXW_AVG.bio_mult_shrub_avg[ind0], _Sep, std_bio_mult_shrub,
				_Sep, SXW_AVG.bio_mult_tree_avg[ind0], _Sep, std_bio_mult_tree,
				_Sep, SXW_AVG.bio_mult_forb_avg[ind0], _Sep, std_bio_mult_forb,
				_Sep, SXW_AVG.wue_mult_grass_avg[ind0], _Sep, std_wue_mult_grass,
				_Sep, SXW_AVG.wue_mult_shrub_avg[ind0], _Sep, std_wue_mult_shrub,
				_Sep, SXW_AVG.wue_mult_tree_avg[ind0], _Sep, std_wue_mult_tree,
				_Sep, SXW_AVG.wue_mult_forb_avg[ind0], _Sep, std_wue_mult_forb);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f",
			_Sep, biomass_grass,
			_Sep, biomass_shrub,
			_Sep, biomass_tree,
			_Sep, biomass_forb,
			_Sep, biomass_total,
			_Sep, biolive_grass,
			_Sep, biolive_shrub,
			_Sep, biolive_tree,
			_Sep, biolive_forb,
			_Sep, biolive_total,
			_Sep, bio_mult_grass,
			_Sep, bio_mult_shrub,
			_Sep, bio_mult_tree,
			_Sep, bio_mult_forb,
			_Sep, wue_mult_grass,
			_Sep, wue_mult_shrub,
			_Sep, wue_mult_tree,
			_Sep, wue_mult_forb);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_estab(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* the establishment check produces, for each species in
	 * the given set, a day of year >=0 that the species
	 * established itself in the current year.  The output
	 * will be a single row of numbers for each year.  Each
	 * column represents a species in the order it was entered
	 * in the estabs.in file.  The value will be the day that
	 * the species established, or 0 if it didn't establish
	 * this year.
	 */
	SW_VEGESTAB *v = &SW_VegEstab;
	IntU i;

	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations) {
		get_outstrleader2(pd);
	}

	for (i = 0; i < v->count; i++)
	{
			switch (pd)
			{
				case eSW_Day:
					p = SW_Model.doy-1;
					break;
				case eSW_Week:
					p = SW_Model.week-tOffset;
					break;
				case eSW_Month:
					p = SW_Model.month-tOffset;
					break;
				case eSW_Year:
					p = 0; // Iypc requires 0 for yearly timeperiod
					break;
			}

			if (isPartialSoilwatOutput == FALSE)
			{
				int
					ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
					ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
				float old_val = SXW_AVG.estab_avg[ind0];

				SXW_AVG.estab_avg[ind0] = get_running_avg(SXW_AVG.estab_avg[ind0], v->parms[i]->estab_doy);

				SXW_AVG.estab_avg[ind1] += get_running_sqr(old_val, v->parms[i]->estab_doy, SXW_AVG.estab_avg[ind0]);
				if(Globals.currIter == Globals.runModelIterations){
					float std_estab = sqrt(SXW_AVG.estab_avg[ind1] / Globals.currIter);

					sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.estab_avg[ind0], _Sep, std_estab);
					strcat(sw_outstr, str);
				}
			}
			if(storeAllIterations){
				sprintf(str_iters, "%c%d", _Sep, v->parms[i]->estab_doy);
				strcat(sw_outstr_agg, str_iters);
			}
	}
}

void get_temp(OutPeriod pd)
{
	SW_WEATHER *v = &SW_Weather;
	RealD v_avg = SW_MISSING;
	RealD v_min = SW_MISSING, v_max = SW_MISSING;
	RealD surfaceTempVal = SW_MISSING;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	switch (pd)
	{
	case eSW_Day:
		v_max = v->dysum.temp_max;
		v_min = v->dysum.temp_min;
		v_avg = v->dysum.temp_avg;
		surfaceTempVal = v->dysum.surfaceTemp;
		break;

	case eSW_Week:
		v_max = v->wkavg.temp_max;
		v_min = v->wkavg.temp_min;
		v_avg = v->wkavg.temp_avg;
		surfaceTempVal = v->wkavg.surfaceTemp;
		break;

	case eSW_Month:
		v_max = v->moavg.temp_max;
		v_min = v->moavg.temp_min;
		v_avg = v->moavg.temp_avg;
		surfaceTempVal = v->moavg.surfaceTemp;
		break;

	case eSW_Year:
		v_max = v->yravg.temp_max;
		v_min = v->yravg.temp_min;
		v_avg = v->yravg.temp_avg;
		surfaceTempVal = v->yravg.surfaceTemp;
		break;
	}

	switch (pd)
	{
		case eSW_Day:
			p = SW_Model.doy-1;
			break;
		case eSW_Week:
			p = SW_Model.week-tOffset;
			break;
		case eSW_Month:
			p = SW_Model.month-tOffset;
			break;
		case eSW_Year:
			p = Globals.currYear-1;
			break;
	}

	if (isPartialSoilwatOutput == FALSE)
	{
		int
			ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
			ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
		float old_val_temp_max = SXW_AVG.max_temp_avg[ind0];
		float old_val_temp_min = SXW_AVG.min_temp_avg[ind0];
		float old_val_temp_avg = SXW_AVG.avg_temp_avg[ind0];
		int old_val_surface = SXW_AVG.surfaceTemp_avg[ind0];

		SXW_AVG.max_temp_avg[ind0] = get_running_avg(SXW_AVG.max_temp_avg[ind0], v_max);
		SXW_AVG.min_temp_avg[ind0] = get_running_avg(SXW_AVG.min_temp_avg[ind0], v_min);
		SXW_AVG.avg_temp_avg[ind0] = get_running_avg(SXW_AVG.avg_temp_avg[ind0], v_avg);
		SXW_AVG.surfaceTemp_avg[ind0] = get_running_avg(SXW_AVG.surfaceTemp_avg[ind0], surfaceTempVal);

		SXW_AVG.max_temp_avg[ind1] += get_running_sqr(old_val_temp_max, v_max, SXW_AVG.max_temp_avg[ind0]);
		SXW_AVG.min_temp_avg[ind1] += get_running_sqr(old_val_temp_min, v_min, SXW_AVG.min_temp_avg[ind0]);
		SXW_AVG.avg_temp_avg[ind1] += get_running_sqr(old_val_temp_avg, v_avg, SXW_AVG.avg_temp_avg[ind0]);
		SXW_AVG.surfaceTemp_avg[ind1] += get_running_sqr(old_val_surface, surfaceTempVal, SXW_AVG.surfaceTemp_avg[ind0]);


   if (Globals.currIter == Globals.runModelIterations){
		 float std_temp_max = sqrt(SXW_AVG.max_temp_avg[ind1] / Globals.currIter);
		 float std_temp_min = sqrt(SXW_AVG.min_temp_avg[ind1] / Globals.currIter);
		 float std_temp_avg = sqrt(SXW_AVG.avg_temp_avg[ind1] / Globals.currIter);
		 float std_surface = sqrt(SXW_AVG.surfaceTemp_avg[ind1] / Globals.currIter);

	   sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, SXW_AVG.max_temp_avg[ind0], _Sep, std_temp_max,
		 	_Sep, SXW_AVG.min_temp_avg[ind0], _Sep, std_temp_min,
		  _Sep, SXW_AVG.avg_temp_avg[ind0], _Sep, std_temp_avg,
			_Sep, SXW_AVG.surfaceTemp_avg[ind0], _Sep, std_surface);
	   strcat(sw_outstr, str);
    }
	}

	if(storeAllIterations){
		sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, v_max, _Sep, v_min, _Sep,
			v_avg, _Sep, surfaceTempVal);
		strcat(sw_outstr_agg, str_iters);
	}

	SXW.temp = v_avg;
	SXW.surfaceTemp = surfaceTempVal;
}

void get_precip(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* 	20091015 (drs) ppt is divided into rain and snow and all three values are output into precip */
	SW_WEATHER *v = &SW_Weather;
	RealD val_ppt = SW_MISSING, val_rain = SW_MISSING, val_snow = SW_MISSING,
			val_snowmelt = SW_MISSING, val_snowloss = SW_MISSING;
	TimeInt p = 0;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	switch(pd)
	{
	case eSW_Day:
		val_ppt = v->dysum.ppt;
		val_rain = v->dysum.rain;
		val_snow = v->dysum.snow;
		val_snowmelt = v->dysum.snowmelt;
		val_snowloss = v->dysum.snowloss;
		break;

	case eSW_Week:
		val_ppt = v->wkavg.ppt;
		val_rain = v->wkavg.rain;
		val_snow = v->wkavg.snow;
		val_snowmelt = v->wkavg.snowmelt;
		val_snowloss = v->wkavg.snowloss;
		break;

	case eSW_Month:
		val_ppt = v->moavg.ppt;
		val_rain = v->moavg.rain;
		val_snow = v->moavg.snow;
		val_snowmelt = v->moavg.snowmelt;
		val_snowloss = v->moavg.snowloss;
		break;

	case eSW_Year:
		val_ppt = v->yravg.ppt;
		val_rain = v->yravg.rain;
		val_snow = v->yravg.snow;
		val_snowmelt = v->yravg.snowmelt;
		val_snowloss = v->yravg.snowloss;
		break;
	}

	switch (pd)
	{
		case eSW_Day:
			p = SW_Model.doy-1;
			break;
		case eSW_Week:
			p = SW_Model.week-tOffset;
			break;
		case eSW_Month:
			p = SW_Model.month-tOffset;
			break;
		case eSW_Year:
			p = 0; // Iypc requires 0 for yearly timeperiod
			break;
	}
	if(isPartialSoilwatOutput == FALSE)
	{
		int
			ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
			ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
		float old_ppt = SXW_AVG.ppt_avg[ind0];
		float old_rain = SXW_AVG.val_rain_avg[ind0];
		float old_snow = SXW_AVG.val_snow_avg[ind0];
		float old_snowmelt = SXW_AVG.val_snowmelt_avg[ind0];
		float old_snowloss = SXW_AVG.val_snowloss_avg[ind0];

		SXW_AVG.ppt_avg[ind0] = get_running_avg(SXW_AVG.ppt_avg[ind0], val_ppt);
		SXW_AVG.val_rain_avg[ind0] = get_running_avg(SXW_AVG.val_rain_avg[ind0], val_rain);
		SXW_AVG.val_snow_avg[ind0] = get_running_avg(SXW_AVG.val_snow_avg[ind0], val_snow);
		SXW_AVG.val_snowmelt_avg[ind0] = get_running_avg(SXW_AVG.val_snowmelt_avg[ind0], val_snowmelt);
		SXW_AVG.val_snowloss_avg[ind0] = get_running_avg(SXW_AVG.val_snowloss_avg[ind0], val_snowloss);

		SXW_AVG.ppt_avg[ind1] += get_running_sqr(old_ppt, val_ppt, SXW_AVG.ppt_avg[ind0]);
		SXW_AVG.val_rain_avg[ind1] += get_running_sqr(old_rain, val_rain, SXW_AVG.val_rain_avg[ind0]);
		SXW_AVG.val_snow_avg[ind1] += get_running_sqr(old_snow, val_snow, SXW_AVG.val_snow_avg[ind0]);
		SXW_AVG.val_snowmelt_avg[ind1] += get_running_sqr(old_snowmelt, val_snowmelt, SXW_AVG.val_snowmelt_avg[ind0]);
		SXW_AVG.val_snowloss_avg[ind1] += get_running_sqr(old_snowloss, val_snowloss, SXW_AVG.val_snowloss_avg[ind0]);

		if(Globals.currIter == Globals.runModelIterations){
			float std_ppt = sqrt(SXW_AVG.ppt_avg[ind1] / Globals.currIter);
			float std_rain = sqrt(SXW_AVG.val_rain_avg[ind1] / Globals.currIter);
			float std_snow = sqrt(SXW_AVG.val_snow_avg[ind1] / Globals.currIter);
			float std_snowmelt = sqrt(SXW_AVG.val_snowmelt_avg[ind1] / Globals.currIter);
			float std_snowloss = sqrt(SXW_AVG.val_snowloss_avg[ind1] / Globals.currIter);

			sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
			  _Sep, SXW_AVG.ppt_avg[ind0], _Sep, std_ppt, _Sep,
				SXW_AVG.val_rain_avg[ind0], _Sep, std_rain, _Sep, SXW_AVG.val_snow_avg[ind0], _Sep, std_snow,
				_Sep, SXW_AVG.val_snowmelt_avg[ind0],
				_Sep, std_snowmelt, _Sep, SXW_AVG.val_snowloss_avg[ind0], _Sep, std_snowloss);
			strcat(sw_outstr, str);
		}
	}
	if(storeAllIterations){
		sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, val_ppt, _Sep,
			val_rain, _Sep, val_snow, _Sep, val_snowmelt, _Sep, val_snowloss);
		strcat(sw_outstr_agg, str_iters);
	}
	SXW.ppt = val_ppt;
}

void get_vwcBulk(OutPeriod pd)
{
	/* --------------------------------------------------- */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	RealD *val = (RealD *) malloc(sizeof(RealD) * SW_Site.n_layers);
	ForEachSoilLayer(i)
		val[i] = SW_MISSING;

	  char str[OUTSTRLEN];
		char str_iters[OUTSTRLEN];
		TimeInt p = 0;
	  if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	  {
	  	get_outstrleader2(pd);
	  }

	switch (pd)
	{ /* vwcBulk at this point is identical to swcBulk */
	case eSW_Day:
		ForEachSoilLayer(i)
			val[i] = v->dysum.vwcBulk[i] / SW_Site.lyr[i]->width;
		break;

	case eSW_Week:
		ForEachSoilLayer(i)
			val[i] = v->wkavg.vwcBulk[i] / SW_Site.lyr[i]->width;
		break;

	case eSW_Month:
		ForEachSoilLayer(i)
			val[i] = v->moavg.vwcBulk[i] / SW_Site.lyr[i]->width;
		break;

	case eSW_Year:
		ForEachSoilLayer(i)
			val[i] = v->yravg.vwcBulk[i] / SW_Site.lyr[i]->width;
		break;
	}

ForEachSoilLayer(i){
	switch (pd)
	{
		case eSW_Day:
			p = SW_Model.doy-1;
			break;
		case eSW_Week:
			p = SW_Model.week-tOffset;
			break;
		case eSW_Month:
			p = SW_Model.month-tOffset;
			break;
		case eSW_Year:
			p = 0; // Iypc/Iylp require 0 for yearly timeperiod
			break;
	}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.vwcbulk_avg[indl0];

			SXW_AVG.vwcbulk_avg[indl0] = get_running_avg(SXW_AVG.vwcbulk_avg[indl0], val[i]);
			SXW_AVG.vwcbulk_avg[indl1] += get_running_sqr(old_val, val[i], SXW_AVG.vwcbulk_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std_vwcbulk = sqrt(SXW_AVG.vwcbulk_avg[indl1] / Globals.currIter);
				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.vwcbulk_avg[indl0], _Sep, std_vwcbulk);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val[i]);
			strcat(sw_outstr_agg, str_iters);
		}
	}

	free(val);
}

void get_vwcMatric(OutPeriod pd)
{
	/* --------------------------------------------------- */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	RealD convert;
	RealD *val = (RealD *) malloc(sizeof(RealD) * SW_Site.n_layers);
	ForEachSoilLayer(i)
		val[i] = SW_MISSING;

	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	/* vwcMatric at this point is identical to swcBulk */
	switch (pd)
	{
	case eSW_Day:
		ForEachSoilLayer(i)
		{
			convert = 1. / (1. - SW_Site.lyr[i]->fractionVolBulk_gravel)
					/ SW_Site.lyr[i]->width;
			val[i] = v->dysum.vwcMatric[i] * convert;
		}
		break;

	case eSW_Week:
		ForEachSoilLayer(i)
		{
			convert = 1. / (1. - SW_Site.lyr[i]->fractionVolBulk_gravel)
					/ SW_Site.lyr[i]->width;
			val[i] = v->wkavg.vwcMatric[i] * convert;
		}
		break;

	case eSW_Month:
		ForEachSoilLayer(i)
		{
			convert = 1. / (1. - SW_Site.lyr[i]->fractionVolBulk_gravel)
					/ SW_Site.lyr[i]->width;
			val[i] = v->moavg.vwcMatric[i] * convert;
		}
		break;

	case eSW_Year:
		ForEachSoilLayer(i)
		{
			convert = 1. / (1. - SW_Site.lyr[i]->fractionVolBulk_gravel)
					/ SW_Site.lyr[i]->width;
			val[i] = v->yravg.vwcMatric[i] * convert;
		}
		break;
	}

switch (pd)
{
	case eSW_Day:
		p = SW_Model.doy-1;
		break;
	case eSW_Week:
		p = SW_Model.week-tOffset;
		break;
	case eSW_Month:
		p = SW_Model.month-tOffset;
		break;
	case eSW_Year:
		p = 0; // Iypc requires 0 for yearly timeperiod
		break;
}
ForEachSoilLayer(i){
	if (isPartialSoilwatOutput == FALSE)
	{
		int
			indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
			indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
		float old_val = SXW_AVG.vwcmatric_avg[indl0];

		SXW_AVG.vwcmatric_avg[indl0] = get_running_avg(SXW_AVG.vwcmatric_avg[indl0], val[i]);
		SXW_AVG.vwcmatric_avg[indl1] += get_running_sqr(old_val, val[i], SXW_AVG.vwcmatric_avg[indl0]);

		if(Globals.currIter == Globals.runModelIterations){
			float std_vwcmatric = sqrt(SXW_AVG.vwcmatric_avg[indl1] / Globals.currIter);
			sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.vwcmatric_avg[indl0], _Sep, std_vwcmatric);
			strcat(sw_outstr, str);
		}
	}
	if(storeAllIterations){
		sprintf(str_iters, "%c%7.6f", _Sep, val[i]);
		strcat(sw_outstr_agg, str_iters);
	}
}
	free(val);
}

/**
**/
void get_swa(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* added 21-Oct-03, cwb */
	TimeInt p = 0;
	int j = 0;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealF val[NVEGTYPES][MAX_LAYERS]; // need 2D array for values

	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;

		if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
			get_outstrleader(pd);

		ForEachSoilLayer(i)
		{
			ForEachVegType(j){
				switch (pd)
				{
					case eSW_Day:
						p = SW_Model.doy-1;
						val[j][i] = v->dysum.SWA_VegType[j][i];
						break;
					case eSW_Week:
						p = SW_Model.week-tOffset;
						val[j][i] = v->wkavg.SWA_VegType[j][i];
						break;
					case eSW_Month:
						p = SW_Model.month-tOffset;
						val[j][i] = v->moavg.SWA_VegType[j][i];
						break;
					case eSW_Year:
						p = Globals.currYear - 1;
						val[j][i] = v->yravg.SWA_VegType[j][i];
						break;
				}
				SXW.sum_dSWA_repartitioned[Ivlp(j,i,p)] = val[j][i];
			}

			if(storeAllIterations){
				sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f",_Sep, val[0][i], _Sep, val[1][i], _Sep, val[2][i], _Sep, val[3][i]);
				strcat(sw_outstr_agg, str_iters);
			}

			if (isPartialSoilwatOutput == FALSE)
			{
				int
					indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
					indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
				// get old average for use in running square
				float old_tree = SXW.SWA_tree_avg[indl0];
				float old_shrub = SXW.SWA_shrub_avg[indl0];
				float old_forb = SXW.SWA_forb_avg[indl0];
				float old_grass = SXW.SWA_grass_avg[indl0];

				// get running average over all iterations
				SXW.SWA_tree_avg[indl0] = get_running_avg(SXW.SWA_tree_avg[indl0], val[0][i]);
				SXW.SWA_shrub_avg[indl0] = get_running_avg(SXW.SWA_shrub_avg[indl0], val[1][i]);
				SXW.SWA_forb_avg[indl0] = get_running_avg(SXW.SWA_forb_avg[indl0], val[2][i]);
				SXW.SWA_grass_avg[indl0] = get_running_avg(SXW.SWA_grass_avg[indl0], val[3][i]);
				//if(i == 7 && p == 0 && Globals.currYear-1 == 0)
					//printf("SXW.sum_dSWA_repartitioned[Ivlp(3,i,p)]: %f\n", val[3][i]);
				// get running square over all iterations. Going to be used to get standard deviation
				SXW.SWA_tree_avg[indl1] = get_running_sqr(old_tree, val[0][i], SXW.SWA_tree_avg[indl0]);
				SXW.SWA_shrub_avg[indl1] = get_running_sqr(old_shrub, val[1][i], SXW.SWA_shrub_avg[indl0]);
				SXW.SWA_forb_avg[indl1] = get_running_sqr(old_forb, val[2][i], SXW.SWA_forb_avg[indl0]);
				SXW.SWA_grass_avg[indl1] = get_running_sqr(old_grass, val[3][i], SXW.SWA_grass_avg[indl0]);

				// divide by number of iterations at end to store average
				if(Globals.currIter == Globals.runModelIterations){
					//if(i == 7 && p == 0 && Globals.currYear-1 == 0)
						//printf("SXW.SWA_grass_avg[Iylp(%d,%d,%d,0)]: %f\n", Globals.currYear-1,i,p, SXW.SWA_grass_avg[indl0]);
					// get standard deviation
					float std_forb = sqrt(SXW.SWA_forb_avg[indl1] / Globals.currIter);
					float std_tree = sqrt(SXW.SWA_tree_avg[indl1] / Globals.currIter);
					float std_shrub = sqrt(SXW.SWA_shrub_avg[indl1] / Globals.currIter);
					float std_grass = sqrt(SXW.SWA_grass_avg[indl1] / Globals.currIter);

					sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",_Sep, SXW.SWA_tree_avg[indl0],
					 	_Sep, std_tree, _Sep, SXW.SWA_shrub_avg[indl0], _Sep, std_shrub, _Sep, SXW.SWA_forb_avg[indl0],
						_Sep, std_forb, _Sep, SXW.SWA_grass_avg[indl0], _Sep, std_grass);
					strcat(sw_outstr, str);
				}
				if (bFlush_output) p++;
			}
		}

}


void get_swcBulk(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* added 21-Oct-03, cwb */
	TimeInt p = 0;
	RealD val = SW_MISSING;
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;

	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
		get_outstrleader2(pd);


	ForEachSoilLayer(i)
	{
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy - 1;
				val = v->dysum.swcBulk[i];
				break; // print current but as index
			case eSW_Week:
				p = SW_Model.week - 1;
				val = v->wkavg.swcBulk[i];
				break;// print previous to current
			case eSW_Month:
				p = SW_Model.month-1;
				val = v->moavg.swcBulk[i];
				break;// print previous to current
			case eSW_Year:
				p = Globals.currYear - 1;
				val = v->yravg.swcBulk[i];
				break;
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = 0.;
			old_val = SXW_AVG.swc_avg[indl0];

			SXW_AVG.swc_avg[indl0] = get_running_avg(SXW_AVG.swc_avg[indl0], val);
			SXW_AVG.swc_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.swc_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float swc_std = sqrt(SXW_AVG.swc_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.swc_avg[indl0], _Sep, swc_std);
				strcat(sw_outstr, str);
			}

		}
		if (bFlush_output) p++;
		SXW.swc[Ilp(i,p)] = val; // SXW.swc[Ilp(layer,timeperiod)]
	}

}

void get_swpMatric(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* can't take arithmetic average of swp because it's
	 * exponential.  At this time (until I remember to look
	 * up whether harmonic or some other average is better
	 * and fix this) we're not averaging swp but converting
	 * the averaged swc.  This also avoids converting for
	 * each day.
	 *
	 * added 12-Oct-03, cwb */

	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;

	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}


	ForEachSoilLayer(i)
	{
	switch (pd)
	{ /* swpMatric at this point is identical to swcBulk */
	case eSW_Day:
		val = SW_SWCbulk2SWPmatric(SW_Site.lyr[i]->fractionVolBulk_gravel,
				v->dysum.swpMatric[i], i);
		p = SW_Model.doy-1;
		break;
	case eSW_Week:
		p = SW_Model.week-tOffset;
		val = SW_SWCbulk2SWPmatric(SW_Site.lyr[i]->fractionVolBulk_gravel,
				v->wkavg.swpMatric[i], i);
		break;
	case eSW_Month:
		p = SW_Model.month-tOffset;
		val = SW_SWCbulk2SWPmatric(SW_Site.lyr[i]->fractionVolBulk_gravel,
				v->moavg.swpMatric[i], i);
		break;
	case eSW_Year:
		p = Globals.currYear - 1;
		val = SW_SWCbulk2SWPmatric(SW_Site.lyr[i]->fractionVolBulk_gravel,
				v->yravg.swpMatric[i], i);
		break;
	}
	if (isPartialSoilwatOutput == FALSE)
	{
		int
			indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
			indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
		float old_val = SXW_AVG.swpmatric_avg[indl0];

		SXW_AVG.swpmatric_avg[indl0] = get_running_avg(SXW_AVG.swpmatric_avg[indl0], val);
		SXW_AVG.swpmatric_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.swpmatric_avg[indl0]);

		if(Globals.currIter == Globals.runModelIterations){
			float std = sqrt(SXW_AVG.swpmatric_avg[indl1] / Globals.currIter);

			sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.swpmatric_avg[indl0], _Sep, std);
			strcat(sw_outstr, str);
		}
	}
	if(storeAllIterations){
		sprintf(str_iters, "%c%7.6f", _Sep, val);
		strcat(sw_outstr_agg, str_iters);
	}
}
}

void get_swaBulk(OutPeriod pd)
{
	/* --------------------------------------------------- */

	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	RealD val = SW_MISSING;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}


	ForEachSoilLayer(i)
	{
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.swaBulk[i];
				break;
			case eSW_Week:
				val = v->wkavg.swaBulk[i];
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val = v->moavg.swaBulk[i];
				break;
			case eSW_Year:
				val = v->yravg.swaBulk[i];
				p = Globals.currYear - 1;
				break;
		}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.swabulk_avg[indl0];

			SXW_AVG.swabulk_avg[indl0] = get_running_avg(SXW_AVG.swabulk_avg[indl0], val);
			SXW_AVG.swabulk_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.swabulk_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.swabulk_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.swabulk_avg[indl0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}

void get_swaMatric(OutPeriod pd)
{
	/* --------------------------------------------------- */

	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	RealD convert;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	ForEachSoilLayer(i)
	{
		convert = 1. / (1. - SW_Site.lyr[i]->fractionVolBulk_gravel);
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.swaMatric[i] * convert;
				break;
			case eSW_Week:
				val = v->wkavg.swaMatric[i] * convert;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				val = v->moavg.swaMatric[i] * convert;
				p = SW_Model.month-tOffset;
				break;
			case eSW_Year:
				val = v->yravg.swaMatric[i] * convert;
				p = Globals.currYear - 1;
				break;
		}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.swamatric_avg[indl0];

			SXW_AVG.swamatric_avg[indl0] = get_running_avg(SXW_AVG.swamatric_avg[indl0], val);
			SXW_AVG.swamatric_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.swamatric_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.swamatric_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.swamatric_avg[indl0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}

void get_surfaceWater(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val_surfacewater = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}


		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val_surfacewater = v->dysum.surfaceWater;
				break;
			case eSW_Week:
				val_surfacewater = v->wkavg.surfaceWater;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val_surfacewater = v->moavg.surfaceWater;
				break;
			case eSW_Year:
				val_surfacewater = v->yravg.surfaceWater;
				p = Globals.currYear - 1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val = SXW_AVG.surfacewater_avg[ind0];

			SXW_AVG.surfacewater_avg[ind0] = get_running_avg(SXW_AVG.surfacewater_avg[ind0], val_surfacewater);
			SXW_AVG.surfacewater_avg[ind1] += get_running_sqr(old_val, val_surfacewater, SXW_AVG.surfacewater_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.surfacewater_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.surfacewater_avg[ind0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val_surfacewater);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_runoffrunon(OutPeriod pd) {
  /* --------------------------------------------------- */
  /* (12/13/2012) (clk) Added function to output runoff variables */

	SW_WEATHER *w = &SW_Weather;
	RealD val_netRunoff = SW_MISSING, val_surfaceRunoff = SW_MISSING,
		val_surfaceRunon = SW_MISSING, val_snowRunoff = SW_MISSING;

	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}


		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val_surfaceRunoff = w->dysum.surfaceRunoff;
				val_surfaceRunon = w->dysum.surfaceRunon;
				val_snowRunoff = w->dysum.snowRunoff;
				break;
			case eSW_Week:
				val_surfaceRunoff = w->wkavg.surfaceRunoff;
				val_surfaceRunon = w->wkavg.surfaceRunon;
				val_snowRunoff = w->wkavg.snowRunoff;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val_surfaceRunoff = w->moavg.surfaceRunoff;
				val_surfaceRunon = w->moavg.surfaceRunon;
				val_snowRunoff = w->moavg.snowRunoff;
				break;
			case eSW_Year:
				p = Globals.currYear-1;
				val_surfaceRunoff = w->yravg.surfaceRunoff;
				val_surfaceRunon = w->yravg.surfaceRunon;
				val_snowRunoff = w->yravg.snowRunoff;
				break;
		}
		val_netRunoff = val_surfaceRunoff + val_snowRunoff - val_surfaceRunon;

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val_total = SXW_AVG.runoff_total_avg[ind0];
			float old_val_surface_runoff = SXW_AVG.surface_runoff_avg[ind0];
			float old_val_surface_runon = SXW_AVG.surface_runon_avg[ind0];
			float old_val_snow = SXW_AVG.runoff_snow_avg[ind0];

			SXW_AVG.runoff_total_avg[ind0] = get_running_avg(SXW_AVG.runoff_total_avg[ind0], val_netRunoff);
			SXW_AVG.runoff_total_avg[ind1] += get_running_sqr(old_val_total, val_netRunoff, SXW_AVG.runoff_total_avg[ind0]);

			SXW_AVG.surface_runoff_avg[ind0] = get_running_avg(SXW_AVG.surface_runoff_avg[ind0], val_surfaceRunoff);
			SXW_AVG.surface_runoff_avg[ind1] += get_running_sqr(old_val_surface_runoff, val_surfaceRunoff, SXW_AVG.surface_runoff_avg[ind0]);

			SXW_AVG.surface_runon_avg[ind0] = get_running_avg(SXW_AVG.surface_runon_avg[ind0], val_surfaceRunon);
			SXW_AVG.surface_runon_avg[ind1] += get_running_sqr(old_val_surface_runon, val_surfaceRunon, SXW_AVG.surface_runon_avg[ind0]);


			SXW_AVG.runoff_snow_avg[ind0] = get_running_avg(SXW_AVG.runoff_snow_avg[ind0], val_snowRunoff);
			SXW_AVG.runoff_snow_avg[ind1] += get_running_sqr(old_val_snow, val_snowRunoff, SXW_AVG.runoff_snow_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std_total = sqrt(SXW_AVG.runoff_total_avg[ind1] / Globals.currIter);
				float std_surface_runoff = sqrt(SXW_AVG.surface_runoff_avg[ind1] / Globals.currIter);
				float std_surface_runon = sqrt(SXW_AVG.surface_runon_avg[ind1] / Globals.currIter);
				float std_snow = sqrt(SXW_AVG.runoff_snow_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, SXW_AVG.runoff_total_avg[ind0], _Sep, std_total,
								_Sep, SXW_AVG.surface_runoff_avg[ind0], _Sep, std_surface_runoff,
								_Sep, SXW_AVG.runoff_snow_avg[ind0], _Sep, std_snow,
								_Sep, SXW_AVG.surface_runon_avg[ind0], _Sep, std_surface_runon);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, val_netRunoff,
		      _Sep, val_surfaceRunoff, _Sep, val_snowRunoff, _Sep, val_surfaceRunon);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_transp(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* 10-May-02 (cwb) Added conditional code to interface
	 *           with STEPPE.
	 */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	RealF *val = (RealF *) malloc(sizeof(RealF) * SW_Site.n_layers); // changed val_total to val
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
		get_outstrleader2(pd);

	RealF *val_total = (RealF *) malloc(sizeof(RealF) * SW_Site.n_layers);
	RealF *val_tree = (RealF *) malloc(sizeof(RealF) * SW_Site.n_layers);
	RealF *val_forb = (RealF *) malloc(sizeof(RealF) * SW_Site.n_layers);
	RealF *val_grass = (RealF *) malloc(sizeof(RealF) * SW_Site.n_layers);
	RealF *val_shrub = (RealF *) malloc(sizeof(RealF) * SW_Site.n_layers);

	ForEachSoilLayer(i){
		val[i] = 0;
		val_total[i] = 0;
		val_tree[i] = 0;
		val_forb[i] = 0;
		val_grass[i] = 0;
		val_shrub[i] = 0;
	}

	get_outstrleader2(pd);
	/* total transpiration */
	ForEachSoilLayer(i)
	{
		switch (pd)
		{
		case eSW_Day:
			val[i] = v->dysum.transp_total[i];
			break;
		case eSW_Week:
			val[i] = v->wkavg.transp_total[i];
			break;
		case eSW_Month:
			val[i] = v->moavg.transp_total[i];
			break;
		case eSW_Year:
			val[i] = v->yravg.transp_total[i];
			break;
		}
	}

	ForEachSoilLayer(i)
	{
		val_total[i] = val[i];
	}

	/* tree-component transpiration */
	ForEachSoilLayer(i)
	{
		switch (pd)
		{
		case eSW_Day:
			val[i] = v->dysum.transp[SW_TREES][i];
			break;
		case eSW_Week:
			val[i] = v->wkavg.transp[SW_TREES][i];
			break;
		case eSW_Month:
			val[i] = v->moavg.transp[SW_TREES][i];
			break;
		case eSW_Year:
			val[i] = v->yravg.transp[SW_TREES][i];
			break;
		}
	}

ForEachSoilLayer(i)
{
	val_tree[i] = val[i];
}

	/* shrub-component transpiration */
	ForEachSoilLayer(i)
	{
		switch (pd)
		{
		case eSW_Day:
			val[i] = v->dysum.transp[SW_SHRUB][i];
			break;
		case eSW_Week:
			val[i] = v->wkavg.transp[SW_SHRUB][i];
			break;
		case eSW_Month:
			val[i] = v->moavg.transp[SW_SHRUB][i];
			break;
		case eSW_Year:
			val[i] = v->yravg.transp[SW_SHRUB][i];
			break;
		}
	}

ForEachSoilLayer(i)
{
	val_shrub[i] = val[i];
}

	/* forb-component transpiration */
	ForEachSoilLayer(i)
	{
		switch (pd)
		{
		case eSW_Day:
			val[i] = v->dysum.transp[SW_FORBS][i];
			break;
		case eSW_Week:
			val[i] = v->wkavg.transp[SW_FORBS][i];
			break;
		case eSW_Month:
			val[i] = v->moavg.transp[SW_FORBS][i];
			break;
		case eSW_Year:
			val[i] = v->yravg.transp[SW_FORBS][i];
			break;
		}
	}

ForEachSoilLayer(i)
{
	val_forb[i] = val[i];
}

	/* grass-component transpiration */
	ForEachSoilLayer(i)
	{
		switch (pd)
		{
		case eSW_Day:
			val[i] = v->dysum.transp[SW_GRASS][i];
			break;
		case eSW_Week:
			val[i] = v->wkavg.transp[SW_GRASS][i];
			break;
		case eSW_Month:
			val[i] = v->moavg.transp[SW_GRASS][i];
			break;
		case eSW_Year:
			val[i] = v->yravg.transp[SW_GRASS][i];
			break;
		}
	}

ForEachSoilLayer(i)
{
	val_grass[i] = val[i];
}

	switch (pd)
	{
		case eSW_Day:
			p = SW_Model.doy - 1;
			break; /* print current but as index */
		case eSW_Week:
			p = SW_Model.week - tOffset;
			break; /* print previous to current */
		case eSW_Month:
			p = SW_Model.month - tOffset;
			break; /* print previous to current */
		case eSW_Year:
			p = 0; // Iypc requires 0 for yearly timeperiod
			break;
	}
	if (bFlush_output) p++;

  ForEachSoilLayer(i)
  {
  /* Pass monthly transpiration values to STEPWAT2 as resources: the
     function `_transp_contribution_by_group` deals with these monthly x layer
     values */
  if (pd == eSW_Month) {
    SXW.transpTotal[Ilp(i,p)] = val_total[i];
    SXW.transpTrees[Ilp(i,p)] = val_tree[i];
    SXW.transpShrubs[Ilp(i,p)] = val_shrub[i];
    SXW.transpForbs[Ilp(i,p)] = val_forb[i];
    SXW.transpGrasses[Ilp(i,p)] = val_grass[i];

    //printf("Tshrubs: bFlush_output=%d pd=%d t=%d lyr=%d ilp=%d T=%.3f\n", bFlush_output, pd, p, i, Ilp(i,p), SXW.transpShrubs[Ilp(i,p)]);
  }

	if (isPartialSoilwatOutput == FALSE)
	{
		int
			indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
			indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
		float old_total = SXW.transpTotal_avg[indl0];
		float old_tree = SXW.transpTrees_avg[indl0];
		float old_shrub = SXW.transpShrubs_avg[indl0];
		float old_forb = SXW.transpForbs_avg[indl0];
		float old_grass = SXW.transpGrasses_avg[indl0];

		// for the average over iteration we need to include year too so values are not overlapping.
		// [Iylp(Globals.currYear-1,i,p)] is a new macro defined in sxw.h that represents year, layer, timeperiod
		SXW.transpTotal_avg[indl0] = get_running_avg(SXW.transpTotal_avg[indl0], SXW.transpTotal[Ilp(i,p)]);
		SXW.transpTrees_avg[indl0] = get_running_avg(SXW.transpTrees_avg[indl0], SXW.transpTrees[Ilp(i,p)]);
		SXW.transpShrubs_avg[indl0] = get_running_avg(SXW.transpShrubs_avg[indl0], SXW.transpShrubs[Ilp(i,p)]);
		SXW.transpForbs_avg[indl0] = get_running_avg(SXW.transpForbs_avg[indl0], SXW.transpForbs[Ilp(i,p)]);
		SXW.transpGrasses_avg[indl0] = get_running_avg(SXW.transpGrasses_avg[indl0], SXW.transpGrasses[Ilp(i,p)]);

		SXW.transpTotal_avg[indl1] = get_running_sqr(old_total, SXW.transpTotal[Ilp(i,p)], SXW.transpTotal_avg[indl0]);
		SXW.transpTrees_avg[indl1] = get_running_sqr(old_tree, SXW.transpTrees[Ilp(i,p)], SXW.transpTrees_avg[indl0]);
		SXW.transpShrubs_avg[indl1] = get_running_sqr(old_shrub, SXW.transpShrubs[Ilp(i,p)], SXW.transpShrubs_avg[indl0]);
		SXW.transpForbs_avg[indl1] = get_running_sqr(old_forb, SXW.transpForbs[Ilp(i,p)], SXW.transpForbs_avg[indl0]);
		SXW.transpGrasses_avg[indl1] = get_running_sqr(old_grass, SXW.transpGrasses[Ilp(i,p)], SXW.transpGrasses_avg[indl0]);

		// if last iteration need to divide by number of iterations to get average over all iterations
		if(Globals.currIter == Globals.runModelIterations){
			float std_total = sqrt(SXW.transpTotal_avg[indl1] / Globals.currIter);
			float std_trees = sqrt(SXW.transpTrees_avg[indl1] / Globals.currIter);
			float std_shrubs = sqrt(SXW.transpShrubs_avg[indl1] / Globals.currIter);
			float std_forbs = sqrt(SXW.transpForbs_avg[indl1] / Globals.currIter);
			float std_grasses = sqrt(SXW.transpGrasses_avg[indl1] / Globals.currIter);

			sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
				_Sep, SXW.transpTotal_avg[indl0], _Sep, std_total, _Sep, SXW.transpTrees_avg[indl0], _Sep,
				std_trees, _Sep, SXW.transpShrubs_avg[indl0], _Sep, std_shrubs, _Sep, SXW.transpForbs_avg[indl0],
				_Sep, std_forbs, _Sep, SXW.transpGrasses_avg[indl0], _Sep, std_grasses);
			strcat(sw_outstr, str);
		}
	}

	if(storeAllIterations){
		sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
			_Sep, val_total[i], _Sep, val_tree[i], _Sep, val_shrub[i], _Sep, val_forb[i], _Sep, val_grass[i]);
		strcat(sw_outstr_agg, str_iters);
	}
}
free(val_total);
free(val_tree);
free(val_forb);
free(val_grass);
free(val_shrub);
	free(val);
}


void get_evapSoil(OutPeriod pd)
{
	/* --------------------------------------------------- */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	ForEachEvapLayer(i)
	{
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.evap[i];
				break;
			case eSW_Week:
				val = v->wkavg.evap[i];
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				val = v->moavg.evap[i];
				p = SW_Model.month-tOffset;
				break;
			case eSW_Year:
				val = v->yravg.evap[i];
				p = Globals.currYear - 1;
				break;
		}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.evapsoil_avg[indl0];

			SXW_AVG.evapsoil_avg[indl0] = get_running_avg(SXW_AVG.evapsoil_avg[indl0], val);
			SXW_AVG.evapsoil_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.evapsoil_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.evapsoil_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.evapsoil_avg[indl0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}

void get_evapSurface(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val_tot = SW_MISSING, val_tree = SW_MISSING, val_forb = SW_MISSING,
			val_shrub = SW_MISSING, val_grass = SW_MISSING, val_litter =
					SW_MISSING, val_water = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val_tot = v->dysum.total_evap;
				val_tree = v->dysum.evap_veg[SW_TREES];
				val_forb = v->dysum.evap_veg[SW_FORBS];
				val_shrub = v->dysum.evap_veg[SW_SHRUB];
				val_grass = v->dysum.evap_veg[SW_GRASS];
				val_litter = v->dysum.litter_evap;
				val_water = v->dysum.surfaceWater_evap;
				break;
			case eSW_Week:
				p = SW_Model.week-tOffset;
				val_tot = v->wkavg.total_evap;
				val_tree = v->wkavg.evap_veg[SW_TREES];
				val_forb = v->wkavg.evap_veg[SW_FORBS];
				val_shrub = v->wkavg.evap_veg[SW_SHRUB];
				val_grass = v->wkavg.evap_veg[SW_GRASS];
				val_litter = v->wkavg.litter_evap;
				val_water = v->wkavg.surfaceWater_evap;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val_tot = v->moavg.total_evap;
				val_tree = v->moavg.evap_veg[SW_TREES];
				val_forb = v->moavg.evap_veg[SW_FORBS];
				val_shrub = v->moavg.evap_veg[SW_SHRUB];
				val_grass = v->moavg.evap_veg[SW_GRASS];
				val_litter = v->moavg.litter_evap;
				val_water = v->moavg.surfaceWater_evap;
				break;
			case eSW_Year:
				p = Globals.currYear - 1;
				val_tot = v->yravg.total_evap;
				val_tree = v->yravg.evap_veg[SW_TREES];
				val_forb = v->yravg.evap_veg[SW_FORBS];
				val_shrub = v->yravg.evap_veg[SW_SHRUB];
				val_grass = v->yravg.evap_veg[SW_GRASS];
				val_litter = v->yravg.litter_evap;
				val_water = v->yravg.surfaceWater_evap;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val_total = SXW_AVG.evapsurface_total_avg[ind0];
			float old_val_tree = SXW_AVG.evapsurface_tree_avg[ind0];
			float old_val_forb = SXW_AVG.evapsurface_forb_avg[ind0];
			float old_val_shrub = SXW_AVG.evapsurface_shrub_avg[ind0];
			float old_val_grass = SXW_AVG.evapsurface_grass_avg[ind0];
			float old_val_litter = SXW_AVG.evapsurface_litter_avg[ind0];
			float old_val_water = SXW_AVG.evapsurface_water_avg[ind0];


			SXW_AVG.evapsurface_total_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_total_avg[ind0], val_tot);
			SXW_AVG.evapsurface_total_avg[ind1] += get_running_sqr(old_val_total, val_tot, SXW_AVG.evapsurface_total_avg[ind0]);

			SXW_AVG.evapsurface_tree_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_tree_avg[ind0], val_tree);
			SXW_AVG.evapsurface_tree_avg[ind1] += get_running_sqr(old_val_tree, val_tree, SXW_AVG.evapsurface_tree_avg[ind0]);

			SXW_AVG.evapsurface_forb_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_forb_avg[ind0], val_forb);
			SXW_AVG.evapsurface_forb_avg[ind1] += get_running_sqr(old_val_forb, val_forb, SXW_AVG.evapsurface_forb_avg[ind0]);

			SXW_AVG.evapsurface_shrub_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_shrub_avg[ind0], val_shrub);
			SXW_AVG.evapsurface_shrub_avg[ind1] += get_running_sqr(old_val_shrub, val_shrub, SXW_AVG.evapsurface_shrub_avg[ind0]);

			SXW_AVG.evapsurface_grass_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_grass_avg[ind0], val_grass);
			SXW_AVG.evapsurface_grass_avg[ind1] += get_running_sqr(old_val_grass, val_grass, SXW_AVG.evapsurface_grass_avg[ind0]);

			SXW_AVG.evapsurface_litter_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_litter_avg[ind0], val_litter);
			SXW_AVG.evapsurface_litter_avg[ind1] += get_running_sqr(old_val_litter, val_litter, SXW_AVG.evapsurface_litter_avg[ind0]);

			SXW_AVG.evapsurface_water_avg[ind0] = get_running_avg(SXW_AVG.evapsurface_water_avg[ind0], val_water);
			SXW_AVG.evapsurface_water_avg[ind1] += get_running_sqr(old_val_water, val_water, SXW_AVG.evapsurface_water_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std_total = sqrt(SXW_AVG.evapsurface_total_avg[ind1] / Globals.currIter);
				float std_tree = sqrt(SXW_AVG.evapsurface_tree_avg[ind1] / Globals.currIter);
				float std_forb = sqrt(SXW_AVG.evapsurface_forb_avg[ind1] / Globals.currIter);
				float std_shrub = sqrt(SXW_AVG.evapsurface_shrub_avg[ind1] / Globals.currIter);
				float std_grass = sqrt(SXW_AVG.evapsurface_grass_avg[ind1] / Globals.currIter);
				float std_litter = sqrt(SXW_AVG.evapsurface_litter_avg[ind1] / Globals.currIter);
				float std_water = sqrt(SXW_AVG.evapsurface_water_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
				 				_Sep, SXW_AVG.evapsurface_total_avg[ind0], _Sep, std_total,
								_Sep, SXW_AVG.evapsurface_tree_avg[ind0], _Sep, std_tree,
								_Sep, SXW_AVG.evapsurface_shrub_avg[ind0], _Sep, std_shrub,
								_Sep, SXW_AVG.evapsurface_forb_avg[ind0], _Sep, std_forb,
								_Sep, SXW_AVG.evapsurface_grass_avg[ind0], _Sep, std_grass,
								_Sep, SXW_AVG.evapsurface_litter_avg[ind0], _Sep, std_litter,
								_Sep, SXW_AVG.evapsurface_water_avg[ind0], _Sep, std_water);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
				_Sep, val_tot, _Sep, val_tree, _Sep, val_shrub, _Sep, val_forb, _Sep, val_grass, _Sep, val_litter, _Sep, val_water);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_interception(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val_tot = SW_MISSING, val_tree = SW_MISSING, val_forb = SW_MISSING,
			val_shrub = SW_MISSING, val_grass = SW_MISSING, val_litter =
					SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

			switch (pd)
			{
				case eSW_Day:
					p = SW_Model.doy-1;
					val_tot = v->dysum.total_int;
					val_tree = v->dysum.int_veg[SW_TREES];
					val_forb = v->dysum.int_veg[SW_FORBS];
					val_shrub = v->dysum.int_veg[SW_SHRUB];
					val_grass = v->dysum.int_veg[SW_GRASS];
					val_litter = v->dysum.litter_int;
					break;
				case eSW_Week:
					p = SW_Model.week-tOffset;
					val_tot = v->wkavg.total_int;
					val_tree = v->wkavg.int_veg[SW_TREES];
					val_forb = v->wkavg.int_veg[SW_FORBS];
					val_shrub = v->wkavg.int_veg[SW_SHRUB];
					val_grass = v->wkavg.int_veg[SW_GRASS];
					val_litter = v->wkavg.litter_int;
					break;
				case eSW_Month:
					p = SW_Model.month-tOffset;
					val_tot = v->moavg.total_int;
					val_tree = v->moavg.int_veg[SW_TREES];
					val_forb = v->moavg.int_veg[SW_FORBS];
					val_shrub = v->moavg.int_veg[SW_SHRUB];
					val_grass = v->moavg.int_veg[SW_GRASS];
					val_litter = v->moavg.litter_int;
					break;
				case eSW_Year:
					val_tot = v->yravg.total_int;
					val_tree = v->yravg.int_veg[SW_TREES];
					val_forb = v->yravg.int_veg[SW_FORBS];
					val_shrub = v->yravg.int_veg[SW_SHRUB];
					val_grass = v->yravg.int_veg[SW_GRASS];
					val_litter = v->yravg.litter_int;
					break;
			}

			if (isPartialSoilwatOutput == FALSE)
			{
				int
					ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
					ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
				float old_val_total, old_val_tree, old_val_forb, old_val_shrub, old_val_grass, old_val_litter = 0.;

				old_val_total = SXW_AVG.interception_total_avg[ind0];
				old_val_tree = SXW_AVG.interception_tree_avg[ind0];
				old_val_shrub = SXW_AVG.interception_shrub_avg[ind0];
				old_val_forb = SXW_AVG.interception_forb_avg[ind0];
				old_val_grass = SXW_AVG.interception_grass_avg[ind0];
				old_val_litter = SXW_AVG.interception_litter_avg[ind0];

				SXW_AVG.interception_total_avg[ind0] = get_running_avg(SXW_AVG.interception_total_avg[ind0], val_tot);
				SXW_AVG.interception_total_avg[ind1] += get_running_sqr(old_val_total, val_tot, SXW_AVG.interception_total_avg[ind0]);

				SXW_AVG.interception_tree_avg[ind0] = get_running_avg(SXW_AVG.interception_tree_avg[ind0], val_tree);
				SXW_AVG.interception_tree_avg[ind1] += get_running_sqr(old_val_tree, val_tree, SXW_AVG.interception_tree_avg[ind0]);

				SXW_AVG.interception_forb_avg[ind0] = get_running_avg(SXW_AVG.interception_forb_avg[ind0], val_forb);
				SXW_AVG.interception_forb_avg[ind1] += get_running_sqr(old_val_forb, val_forb, SXW_AVG.interception_forb_avg[ind0]);

				SXW_AVG.interception_shrub_avg[ind0] = get_running_avg(SXW_AVG.interception_shrub_avg[ind0], val_shrub);
				SXW_AVG.interception_shrub_avg[ind1] += get_running_sqr(old_val_shrub, val_shrub, SXW_AVG.interception_shrub_avg[ind0]);

				SXW_AVG.interception_grass_avg[ind0] = get_running_avg(SXW_AVG.interception_grass_avg[ind0], val_grass);
				SXW_AVG.interception_grass_avg[ind1] += get_running_sqr(old_val_grass, val_grass, SXW_AVG.interception_grass_avg[ind0]);

				SXW_AVG.interception_litter_avg[ind0] = get_running_avg(SXW_AVG.interception_litter_avg[ind0], val_litter);
				SXW_AVG.interception_litter_avg[ind1] += get_running_sqr(old_val_litter, val_litter, SXW_AVG.interception_litter_avg[ind0]);

				if(Globals.currIter == Globals.runModelIterations){
					float std_total = 0., std_tree = 0., std_forb = 0., std_shrub = 0., std_grass = 0., std_litter = 0.;

					std_total = sqrt(SXW_AVG.interception_total_avg[ind1] / Globals.currIter);
					std_tree = sqrt(SXW_AVG.interception_tree_avg[ind1] / Globals.currIter);
					std_forb = sqrt(SXW_AVG.interception_forb_avg[ind1] / Globals.currIter);
					std_shrub = sqrt(SXW_AVG.interception_shrub_avg[ind1] / Globals.currIter);
					std_grass = sqrt(SXW_AVG.interception_grass_avg[ind1] / Globals.currIter);
					std_litter = sqrt(SXW_AVG.interception_litter_avg[ind1] / Globals.currIter);

					sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
					 				_Sep, SXW_AVG.interception_total_avg[ind0], _Sep, std_total,
									_Sep, SXW_AVG.interception_tree_avg[ind0], _Sep, std_tree,
									_Sep, SXW_AVG.interception_shrub_avg[ind0], _Sep, std_shrub,
									_Sep, SXW_AVG.interception_forb_avg[ind0], _Sep, std_forb,
									_Sep, SXW_AVG.interception_grass_avg[ind0], _Sep, std_grass,
									_Sep, SXW_AVG.interception_litter_avg[ind0], _Sep, std_litter);
					strcat(sw_outstr, str);
				}
			}
			if(storeAllIterations){
				sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, val_tot, _Sep, val_tree, _Sep, val_shrub, _Sep, val_forb, _Sep, val_grass, _Sep, val_litter);
				strcat(sw_outstr_agg, str_iters);
			}
}

void get_soilinf(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* 20100202 (drs) added */
	/* 20110219 (drs) added runoff */
	/* 12/13/2012	(clk)	moved runoff, now named snowRunoff, to get_runoffrunon(); */
	SW_WEATHER *v = &SW_Weather;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val_inf = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val_inf = v->dysum.soil_inf;
				break;
			case eSW_Week:
				val_inf = v->wkavg.soil_inf;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val_inf = v->moavg.soil_inf;
				break;
			case eSW_Year:
				p = Globals.currYear - 1;
				val_inf = v->yravg.soil_inf;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val = SXW_AVG.soilinfilt_avg[ind0];

			SXW_AVG.soilinfilt_avg[ind0] = get_running_avg(SXW_AVG.soilinfilt_avg[ind0], val_inf);
			SXW_AVG.soilinfilt_avg[ind1] += get_running_sqr(old_val, val_inf, SXW_AVG.soilinfilt_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.soilinfilt_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.soilinfilt_avg[ind0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val_inf);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_lyrdrain(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* 20100202 (drs) added */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	for (i = 0; i < SW_Site.n_layers - 1; i++)
	{
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.lyrdrain[i];
				break;
			case eSW_Week:
				val = v->wkavg.lyrdrain[i];
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				val = v->moavg.lyrdrain[i];
				p = SW_Model.month-tOffset;
				break;
			case eSW_Year:
				p = Globals.currYear - 1;
				val = v->yravg.lyrdrain[i];
				break;
		}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.lyrdrain_avg[indl0];

			SXW_AVG.lyrdrain_avg[indl0] = get_running_avg(SXW_AVG.lyrdrain_avg[indl0], val);
			SXW_AVG.lyrdrain_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.lyrdrain_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.lyrdrain_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.lyrdrain_avg[indl0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}

void get_hydred(OutPeriod pd)
{
	/* --------------------------------------------------- */
	/* 20101020 (drs) added */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val_total = SW_MISSING;
	RealD val_tree = SW_MISSING;
	RealD val_shrub = SW_MISSING;
	RealD val_forb = SW_MISSING;
	RealD val_grass = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	ForEachSoilLayer(i)
	{
			switch (pd)
			{
			case eSW_Day:
				val_total = v->dysum.hydred_total[i];
				val_tree = v->dysum.hydred[SW_TREES][i];
				val_shrub = v->dysum.hydred[SW_SHRUB][i];
				val_grass = v->dysum.hydred[SW_GRASS][i];
				val_forb = v->dysum.hydred[SW_FORBS][i];
				p = SW_Model.doy-1;
				break;
			case eSW_Week:
				val_total = v->wkavg.hydred_total[i];
				val_tree = v->wkavg.hydred[SW_TREES][i];
				val_shrub = v->wkavg.hydred[SW_SHRUB][i];
				val_grass = v->wkavg.hydred[SW_GRASS][i];
				val_forb = v->wkavg.hydred[SW_FORBS][i];
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				val_total = v->moavg.hydred_total[i];
				val_tree = v->moavg.hydred[SW_TREES][i];
				val_shrub = v->moavg.hydred[SW_SHRUB][i];
				val_grass = v->moavg.hydred[SW_GRASS][i];
				val_forb = v->moavg.hydred[SW_FORBS][i];
				p = SW_Model.month-tOffset;
				break;
			case eSW_Year:
				val_total = v->yravg.hydred_total[i];
				val_tree = v->yravg.hydred[SW_TREES][i];
				val_shrub = v->yravg.hydred[SW_SHRUB][i];
				val_grass = v->yravg.hydred[SW_GRASS][i];
				val_forb = v->yravg.hydred[SW_FORBS][i];
				p = Globals.currYear - 1;
				break;
			}
		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val_total = SXW_AVG.hydred_total_avg[indl0];
			float old_val_tree = SXW_AVG.hydred_tree_avg[indl0];
			float old_val_forb = SXW_AVG.hydred_forb_avg[indl0];
			float old_val_shrub = SXW_AVG.hydred_shrub_avg[indl0];
			float old_val_grass = SXW_AVG.hydred_grass_avg[indl0];

			SXW_AVG.hydred_total_avg[indl0] = get_running_avg(SXW_AVG.hydred_total_avg[indl0], val_total);
			SXW_AVG.hydred_tree_avg[indl0] = get_running_avg(SXW_AVG.hydred_tree_avg[indl0], val_tree);
			SXW_AVG.hydred_shrub_avg[indl0] = get_running_avg(SXW_AVG.hydred_shrub_avg[indl0], val_shrub);
			SXW_AVG.hydred_forb_avg[indl0] = get_running_avg(SXW_AVG.hydred_forb_avg[indl0], val_forb);
			SXW_AVG.hydred_grass_avg[indl0] = get_running_avg(SXW_AVG.hydred_grass_avg[indl0], val_grass);

			SXW_AVG.hydred_total_avg[indl1] += get_running_sqr(old_val_total, val_total, SXW_AVG.hydred_total_avg[indl0]);
			SXW_AVG.hydred_tree_avg[indl1] += get_running_sqr(old_val_tree, val_tree, SXW_AVG.hydred_tree_avg[indl0]);
			SXW_AVG.hydred_shrub_avg[indl1] += get_running_sqr(old_val_shrub, val_shrub, SXW_AVG.hydred_shrub_avg[indl0]);
			SXW_AVG.hydred_forb_avg[indl1] += get_running_sqr(old_val_forb, val_forb, SXW_AVG.hydred_forb_avg[indl0]);
			SXW_AVG.hydred_grass_avg[indl1] += get_running_sqr(old_val_grass, val_grass, SXW_AVG.hydred_grass_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std_total = sqrt(SXW_AVG.hydred_total_avg[indl1] / Globals.currIter);
				float std_tree = sqrt(SXW_AVG.hydred_tree_avg[indl1] / Globals.currIter);
				float std_forb = sqrt(SXW_AVG.hydred_forb_avg[indl1] / Globals.currIter);
				float std_shrub = sqrt(SXW_AVG.hydred_shrub_avg[indl1] / Globals.currIter);
				float std_grass = sqrt(SXW_AVG.hydred_grass_avg[indl1] / Globals.currIter);


				sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f",
								_Sep, SXW_AVG.hydred_total_avg[indl0], _Sep, std_total,
								_Sep, SXW_AVG.hydred_tree_avg[indl0], _Sep, std_tree,
								_Sep, SXW_AVG.hydred_shrub_avg[indl0], _Sep, std_shrub,
								_Sep, SXW_AVG.hydred_forb_avg[indl0], _Sep, std_forb,
								_Sep, SXW_AVG.hydred_grass_avg[indl0], _Sep, std_grass
							);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, val_total, _Sep, val_tree, _Sep, val_shrub, _Sep, val_forb, _Sep, val_grass);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}

void get_aet(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
		get_outstrleader2(pd);

	switch (pd)
	{
	case eSW_Day:
		val = v->dysum.aet;
		break;
	case eSW_Week:
		val = v->wkavg.aet;
		break;
	case eSW_Month:
		val = v->moavg.aet;
		break;
	case eSW_Year:
		val = v->yravg.aet;
		break;
	}

	switch (pd)
	{
		case eSW_Day:
			p = SW_Model.doy-1;
			break;
		case eSW_Week:
			p = SW_Model.week-tOffset;
			break;
		case eSW_Month:
			p = SW_Model.month-tOffset;
			break;
		case eSW_Year:
			p = Globals.currYear - 1;
			break;
	}

	if (isPartialSoilwatOutput == FALSE)
	{
		int
			ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
			ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
		float old_val = SXW_AVG.aet_avg[ind0];

		//running aet_avg
		SXW_AVG.aet_avg[ind0] = get_running_avg(SXW_AVG.aet_avg[ind0], val);
		SXW_AVG.aet_avg[ind1] += get_running_sqr(old_val, val, SXW_AVG.aet_avg[ind0]);

		if(Globals.currIter == Globals.runModelIterations){
			float std_aet = sqrt(SXW_AVG.aet_avg[ind1] / Globals.currIter);

			sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.aet_avg[ind0], _Sep, std_aet);
			strcat(sw_outstr, str);
		}
	}

	if(storeAllIterations){
		sprintf(str_iters, "%c%7.6f", _Sep, val);
		strcat(sw_outstr_agg, str_iters);
	}
		SXW.aet += val;
}

void get_pet(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;

	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;

	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.pet;
				break;
			case eSW_Week:
				val = v->wkavg.pet;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val = v->moavg.pet;
				break;
			case eSW_Year:
				val = v->yravg.pet;
				p = Globals.currYear-1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val = SXW_AVG.pet_avg[ind0];

			SXW_AVG.pet_avg[ind0] = get_running_avg(SXW_AVG.pet_avg[ind0], val);
			SXW_AVG.pet_avg[ind1] += get_running_sqr(old_val, val, SXW_AVG.pet_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.pet_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.pet_avg[ind0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_wetdays(OutPeriod pd)
{
	/* --------------------------------------------------- */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

	int val = 99;
	ForEachSoilLayer(i){
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = (v->is_wet[i]) ? 1 : 0;
				break;
			case eSW_Week:
				val = (int) v->wkavg.wetdays[i];
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val = (int) v->moavg.wetdays[i];
				break;
			case eSW_Year:
				val = (int) v->yravg.wetdays[i];
				p = Globals.currYear - 1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.wetday_avg[indl0];

			SXW_AVG.wetday_avg[indl0] = get_running_avg(SXW_AVG.wetday_avg[indl0], val);
			SXW_AVG.wetday_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.wetday_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.wetday_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%i%c%i", _Sep, (int)SXW_AVG.wetday_avg[indl0], _Sep, (int)std); // cast to int for proper output format
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%i", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}

void get_snowpack(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val_swe = SW_MISSING, val_depth = SW_MISSING;
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val_swe = v->dysum.snowpack;
				val_depth = v->dysum.snowdepth;
				break;
			case eSW_Week:
				val_swe = v->wkavg.snowpack;
				val_depth = v->wkavg.snowdepth;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val_swe = v->moavg.snowpack;
				val_depth = v->moavg.snowdepth;
				break;
			case eSW_Year:
				val_swe = v->yravg.snowpack;
				val_depth = v->yravg.snowdepth;
				p = Globals.currYear - 1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val_swe = SXW_AVG.snowpack_water_eqv_avg[ind0];
			float old_val_depth = SXW_AVG.snowpack_depth_avg[ind0];

			SXW_AVG.snowpack_water_eqv_avg[ind0] = get_running_avg(SXW_AVG.snowpack_water_eqv_avg[ind0], val_swe);
			SXW_AVG.snowpack_water_eqv_avg[ind1] += get_running_sqr(old_val_swe, val_swe, SXW_AVG.snowpack_water_eqv_avg[ind0]);

			SXW_AVG.snowpack_depth_avg[ind0] = get_running_avg(SXW_AVG.snowpack_depth_avg[ind0], val_depth);
			SXW_AVG.snowpack_depth_avg[ind1] += get_running_sqr(old_val_depth, val_depth, SXW_AVG.snowpack_depth_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std_swe = sqrt(SXW_AVG.snowpack_water_eqv_avg[ind1] / Globals.currIter);
				float std_depth = sqrt(SXW_AVG.snowpack_depth_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f%c%7.6f%c%7.6f", _Sep, SXW_AVG.snowpack_water_eqv_avg[ind0], _Sep, std_swe,
								_Sep, SXW_AVG.snowpack_depth_avg[ind0], _Sep, std_depth);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f%c%7.6f", _Sep, val_swe, _Sep, val_depth);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_deepswc(OutPeriod pd)
{
	/* --------------------------------------------------- */
	SW_SOILWAT *v = &SW_Soilwat;
	char str[OUTSTRLEN];
	char str_iters[OUTSTRLEN];
	RealD val = SW_MISSING;
	TimeInt p = 0;
	if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
	{
		get_outstrleader2(pd);
	}

		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.deep;
				break;
			case eSW_Week:
				val = v->wkavg.deep;
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val = v->moavg.pet;
				break;
			case eSW_Year:
				val = v->yravg.deep;
				p = Globals.currYear - 1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				ind0 = Iypc(Globals.currYear - 1, p, 0, pd), // index for mean
				ind1 = Iypc(Globals.currYear - 1, p, 1, pd); // index for sd
			float old_val = SXW_AVG.deepswc_avg[ind0];

			SXW_AVG.deepswc_avg[ind0] = get_running_avg(SXW_AVG.deepswc_avg[ind0], val);
			SXW_AVG.deepswc_avg[ind1] += get_running_sqr(old_val, val, SXW_AVG.deepswc_avg[ind0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.deepswc_avg[ind1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.deepswc_avg[ind0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
}

void get_soiltemp(OutPeriod pd)
{
	/* --------------------------------------------------- */
	LyrIndex i;
	SW_SOILWAT *v = &SW_Soilwat;
		char str[OUTSTRLEN];
		char str_iters[OUTSTRLEN];
		RealD val = SW_MISSING;
		TimeInt p = 0;
		if ((isPartialSoilwatOutput == FALSE && Globals.currIter == Globals.runModelIterations) || storeAllIterations)
		{
			get_outstrleader2(pd);
		}

	ForEachSoilLayer(i){
		switch (pd)
		{
			case eSW_Day:
				p = SW_Model.doy-1;
				val = v->dysum.sTemp[i];
				break;
			case eSW_Week:
				val = v->wkavg.sTemp[i];
				p = SW_Model.week-tOffset;
				break;
			case eSW_Month:
				p = SW_Model.month-tOffset;
				val = v->moavg.sTemp[i];
				break;
			case eSW_Year:
				val = v->yravg.sTemp[i];
				p = Globals.currYear - 1;
				break;
		}

		if (isPartialSoilwatOutput == FALSE)
		{
			int
				indl0 = Iylp(Globals.currYear - 1, i, p, pd, 0), // index for mean
				indl1 = Iylp(Globals.currYear - 1, i, p, pd, 1); // index for sd
			float old_val = SXW_AVG.soiltemp_avg[indl0];

			SXW_AVG.soiltemp_avg[indl0] = get_running_avg(SXW_AVG.soiltemp_avg[indl0], val);
			SXW_AVG.soiltemp_avg[indl1] += get_running_sqr(old_val, val, SXW_AVG.soiltemp_avg[indl0]);

			if(Globals.currIter == Globals.runModelIterations){
				float std = sqrt(SXW_AVG.soiltemp_avg[indl1] / Globals.currIter);

				sprintf(str, "%c%7.6f%c%7.6f", _Sep, SXW_AVG.soiltemp_avg[indl0], _Sep, std);
				strcat(sw_outstr, str);
			}
		}
		if(storeAllIterations){
			sprintf(str_iters, "%c%7.6f", _Sep, val);
			strcat(sw_outstr_agg, str_iters);
		}
	}
}
