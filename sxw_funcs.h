/********************************************************/
/********************************************************/
/*  Source file: sxw_funcs.h
 *  Type: header
 *  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *  Purpose: Separates the declarations of the STEPWAT
 *           function declarations. */
/*  History:
 *     (14-Apr-2002) -- INITIAL CODING - cwb */
/********************************************************/
/********************************************************/

#ifndef SXW_FUNCS_DEF
#define SXW_FUNCS_DEF

#include "sxw.h"

RealF SXW_GetTranspiration( GrpIndex rg);
int get_SW2_veg_index(int veg_prod_type);

void SXW_Init( Bool init_SW, char *f_roots );
void SXW_Reset(char* SOILWAT_file);
void SXW_Run_SOILWAT (void);
void SXW_InitPlot (void);
void SXW_PrintDebug(Bool cleanup) ;

#ifdef DEBUG_MEM
 void SXW_SetMemoryRefs(void);
#endif


#endif
