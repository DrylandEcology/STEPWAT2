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
#include "sxw.h"
#include "sxw_module.h"
#include "sw_src/SW_Weather.h"

/*************** Global Variable Declarations ***************/
/***********************************************************/
extern SXW_t* SXW;

extern SW_SITE SW_Site;
extern SW_MODEL SW_Model;
extern SW_SOILWAT SW_Soilwat;
extern SW_WEATHER SW_Weather;

extern EnvType *Env;

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


