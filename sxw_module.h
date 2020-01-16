/**
 * \file sxw_module.h
 * \brief Contains declarations for relevant \ref SXW functions from several
 *        different source files.
 * 
 * Ideally this file would be incorporated into \ref sxw.h to make one header
 * file for the entire SXW module. 
 * 
 * \author CWB (initial coding)
 * \date 22 May 2002
 * \ingroup SXW
 */

#ifndef SXW_MODULE_DEF
#define SXW_MODULE_DEF

#include "sw_src/SW_Control.h"
#include "sw_src/SW_Model.h"
#include "sw_src/SW_VegProd.h"
#include "sw_src/SW_SoilWater.h"
#include "sw_src/SW_Files.h"

/* some macros for the production conversion array */
#define PC_Bmass 0
#define PC_Litter 1
#define PC_Live 2


/* These functions are found in sxw_resource.c */
void _sxw_root_phen(void);
void _sxw_update_resource(void);
void _sxw_update_root_tables( RealF sizes[] );


/* These functions are found in sxw_soilwat.c */
void  _sxw_sw_setup(RealF sizes[]);
void  _sxw_sw_run(void);
void  _sxw_sw_clear_transp(void);

/* These functions are found in sxw_environs.c */
void _sxw_set_environs(void);

/* testing code-- see sxw_tester.c */
void _sxw_test(void);

//sql
void connect(char *debugout);
void createTables(void);
void disconnect(void);
void insertInfo(void);
void insertRootsXphen(double * _rootsXphen);
void insertSXWPhen(void);
void insertSXWProd(void);
void insertInputVars(void);
void insertInputProd(void);
void insertInputSoils(void);
void insertOutputVars(RealF * _resource_cur, RealF added_transp);
void insertRgroupInfo(RealF * _resource_cur);
void insertOutputProd(SW_VEGPROD *v);
void insertRootsSum(RealD * _roots_active_sum);
void insertRootsRelative(RealD * _roots_active_rel);
void insertTranspiration(void);
void insertSWCBulk(void);


#endif
