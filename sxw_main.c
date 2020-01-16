/********************************************************/
/********************************************************/
/*  Source file: sxw_main.c
 *
/*  Type: testing module
 *
/*  Purpose: Contains main() for suite of code to test
 *           the sxw module.
 *
/*  Calls:
 *
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *
/*  History:
/*     (25-Oct-2002) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

/* =================================================== */
/*                INCLUDES / DEFINES                   */
/* --------------------------------------------------- */


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ST_steppe.h"
#include "filefuncs.h"
#include "myMemory.h"
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
