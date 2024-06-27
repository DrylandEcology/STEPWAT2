
/**
 * \file ST_globals.h
 * \brief Defines the globally available variables.
 * 
 * Defines the globally available variables to access/manipulate the various
 * "objects". Except for \ref ST_main.c, all modules requiring access to a
 * global variable will reference this file. \ref ST_main.c module actually 
 * declares these variables.
 * 
 * \author CWB (initial programming)
 * \date 15 June 2000
 * \author Chandler Haukap (changed all variables to pointers)
 * \date August 2019
 * \ingroup STEPPE
 */

#ifndef GLOBALS_H
#define GLOBALS_H

#include "ST_defines.h"
#include "ST_functions.h"
#include "sw_src/include/SW_datastructs.h"

extern SpeciesType  **Species;
extern GroupType    **RGroup;
extern SucculentType  *Succulent;
extern EnvType        *Env;
extern PlotType       *Plot;
extern ModelType      *Globals;
extern BmassFlagsType BmassFlags;
extern MortFlagsType  MortFlags;
extern GlobalType     SuperGlobals;
extern SW_DOMAIN      SoilWatDomain;
extern SW_ALL SoilWatAll;
extern SW_OUTPUT_POINTERS SoilWatOutputPtrs[SW_OUTNKEYS];
extern LOG_INFO LogInfo;

extern Bool UseGrid;
extern Bool UseProgressBar;


#endif
