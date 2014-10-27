/********************************************************/
/********************************************************/
/*  Source file: sxw_environs.c
 *
/*  Type: module
 *
/*  Purpose: Precip and temperature values come from SOILWAT
 *           to be added to STEPPE.
 *
/*  Dependency:  sxw.c
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *
/*  History:
/*     (21-May-2002) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include "generic.h"
#include "ST_steppe.h"
/*#include "ST_globals.h"*/
#include "SW_Defines.h"
#include "sxw.h"
#include "sxw_module.h"
#include "SW_Model.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_Weather.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
extern SXW_t SXW;

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_SOILWAT SW_Soilwat;
extern SW_WEATHER SW_Weather;

extern EnvType Env;

/*************** Local Variable Declarations ***************/
/***********************************************************/

/*************** Local Function Declarations ***************/
/***********************************************************/

/***********************************************************/
/****************** Begin Function Code ********************/
/***********************************************************/

void _sxw_set_environs(void) {
/*======================================================*/
  /* run this after SOILWAT has completed.
   * need to convert ppt from cm to mm
  */
	//Integer conversion always does floor so .5 is added
	Env.ppt   = (IntS) (SXW.ppt * 10 + .5);
	Env.temp  = SXW.temp;
}


