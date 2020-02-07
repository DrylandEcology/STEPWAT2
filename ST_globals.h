/********************************************************/
/********************************************************/
/*  Source file: ST_globals.h
 *  Type: header
 *  Application: STEPPE - plant community dynamics simulator
 *  Purpose: Defines the globally available variables to
 *           access/manipulate the various "objects".
 *           Except for main.c, all modules requiring
 *           access to a module will reference this file.
 *           The main.c module actually declares these
 *           variables.
 *  History:
 *     (6/15/2000) -- INITIAL CODING - cwb
 */
/********************************************************/
/********************************************************/


#include "ST_defines.h"
#include "ST_functions.h"

extern SpeciesType  **Species;
extern GroupType    **RGroup;
extern SucculentType  *Succulent;
extern EnvType        *Env;
extern PlotType       *Plot;
extern ModelType      *Globals;
extern BmassFlagsType BmassFlags;
extern MortFlagsType  MortFlags;
extern GlobalType     SuperGlobals;

extern Bool UseGrid;
extern Bool UseCheatgrassWildfire;