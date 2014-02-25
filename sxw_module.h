/********************************************************/
/********************************************************/
/*  Source file: sxw_module.h
/*  Type: header
/*  Application: STEPWAT - plant community dynamics simulator
 *               coupled with the  SOILWAT model.
 *
/*  Purpose: Contains declarations relevant for the SXW_
 *           "module" made up of several source files.
 *
/*  Applies to: sxw.c sxw_steppe.c sxw_soilwat.c
 *
/*  History:
/*     (22-May-2002) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/

#ifndef SXW_MODULE_DEF
#define SXW_MODULE_DEF


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




#endif
