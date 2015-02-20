/********************************************************/
/********************************************************/
/*  Source file: ST_steppe.h
/*  Type: header
/*  Application: STEPPE - plant community dynamics simulator
/*  Purpose: This header file puts the commonly needed
 *           declarations in one place.
/*  History:
/*     (6/15/2000) -- INITIAL CODING - cwb
/*
/********************************************************/
/********************************************************/


#include "ST_defines.h"
#include "ST_functions.h"

#include "generic.h"

#define DFLT_FIRSTFILE "files.in"

void ST_connect(char *stdbName);
void ST_disconnect(void);
void insertIndivKill(int IndivID, int KillTypeID);
void insertIndivYearInfo(IndivType *ind);
void insertIndiv(IndivType *ind);
void insertSpecieYearInfo(SppIndex s);
void insertRGroupYearInfo(GrpIndex g);
