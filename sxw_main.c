/**
 * \file sxw_main.c
 * \brief A deprecated testing file.
 *
 * This file exists only for reference. It is completely deprecated and will
 * not compile. At some point it should be removed.
 *
 * \author CWB (initial coding)
 * \date 25 October 2002
 */

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ST_steppe.h"
#include "sw_src/filefuncs.h"
#include "sw_src/myMemory.h"
#include "sxw_funcs.h"

/* for the chdir function */
#ifdef __BCC__
  #include <dir.h>
#else
  #include <unistd.h>
#endif


/************ External Variable Definitions  ***************/
/*              see ST_globals.h                       */
/***********************************************************/
char errstr[1024];
char inbuf[1024];
FILE *logfp;
int logged;  /* indicator that err file was written to */
SpeciesType  **Species;
GroupType    **RGroup;
SucculentType  Succulent;
EnvType        Env;
PlotType       Plot;
ModelType      Globals;
//FilesType      Files;
BmassFlagsType BmassFlags;
MortFlagsType  MortFlags;


/******************** Begin Code ***************************/
/***********************************************************/
void main( int argc, char **argv) {

  Parms_Initialize( 0);

  SXW_Init();

}
