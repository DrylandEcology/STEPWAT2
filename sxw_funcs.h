/**
 * \file sxw_funcs.h
 * \brief Declares some of the \ref SXW functions. 
 * 
 * Ideally this file would be merged into \ref sxw.h to put all exported 
 * functions in one place.
 * 
 * \author CWB (initial programming)
 * \date 14 April 2002
 * \author Chandler Haukap (author of this documentation)
 * \date 16 January 2020
 * \ingroup SXW
 */

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
