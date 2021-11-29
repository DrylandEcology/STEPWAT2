/**
 * \file: sxw_environs.c
 * \brief Generates [SOILWAT2](\ref sw_src) precipitation information and 
 *        passes it to [Steppe](\ref STEPPE).
 * 
 * \author CWB (inital programming)
 * \date 21 May 2002
 * \ingroup ENVIRONMENT
 * \ingroup SXW_PRIVATE
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

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


