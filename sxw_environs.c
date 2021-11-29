/********************************************************/
/********************************************************/
/*  Source file: sxw_environs.c
 *  Type: module
 *  Purpose: Precip and temperature values come from SOILWAT
 *           to be added to STEPPE.
 *  Dependency:  sxw.c
 *  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model. */
/*  History:
 *     (21-May-2002) -- INITIAL CODING - cwb */
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <stdio.h>
#include "sw_src/generic.h"
#include "ST_steppe.h"
#include "ST_globals.h" // externs `*Env`
#include "sw_src/SW_Defines.h"
#include "sxw.h" // externs `*SXW`
#include "sxw_module.h"
#include "sw_src/SW_Model.h" // externs SW_Model
#include "sw_src/SW_Site.h" // externs SW_Site
#include "sw_src/SW_SoilWater.h" // externs SW_Soilwat
#include "sw_src/SW_Weather.h" // externs SW_Weather



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
	Env->ppt   = (IntS) (SXW->ppt * 10 + .5);
	Env->temp  = SXW->temp;
}


