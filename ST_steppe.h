/**
 * \file ST_steppe.h
 * \brief Contains some function definitions from the [steppe](\ref STEPPE)
 *        module.
 * 
 * These functions are declared here and defined in \ref ST_sql.c. This file 
 * should be renamed ST_sql.c, but file renaming is outside the scope of my
 * surrent issue.
 * 
 * \author CWB (author of the code)
 * \author Chandler Haukap (author of this file description)
 * \date 15 June 2000
 * \ingroup SQL
 */

#ifndef STEPPE_H
#define STEPPE_H

#include "ST_functions.h"

#define DFLT_FIRSTFILE "files.in"

void ST_connect(char *stdbName);
void ST_disconnect(void);
void insertIndivKill(int IndivID, int KillTypeID);
void insertIndivYearInfo(IndivType *ind);
void insertIndiv(IndivType *ind);
void insertSpecieYearInfo(SppIndex s);
void insertRGroupYearInfo(GrpIndex g);

#endif