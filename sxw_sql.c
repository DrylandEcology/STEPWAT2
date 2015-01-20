/*
 * sxw_sql.c
 *
 *  Created on: Jan 15, 2015
 *      Author: Ryan J. Murphy
 */

#include <stdio.h>
#include <string.h>
#include <sqlite3.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "SW_Defines.h"
#include "sxw_module.h"
#include "sxw.h"

extern SW_MODEL SW_Model;
extern SXW_t SXW;
extern SW_SITE SW_Site;
extern SW_VEGPROD SW_VegProd;
extern RealD *_phen;
extern RealF _prod_conv[MAX_MONTHS][3];

static sqlite3 *db;
static char sql[1024];

static sqlite3_stmt * stmt_InputVars;
static sqlite3_stmt * stmt_InputProd;
static sqlite3_stmt * stmt_InputSoils;
static sqlite3_stmt * stmt_OutVars;
static sqlite3_stmt * stmt_OutRgroup;
static sqlite3_stmt * stmt_OutProd;
static sqlite3_stmt * stmt_OutRootsSum;
static sqlite3_stmt * stmt_OutRootsRel;
static sqlite3_stmt * stmt_OutTransp;
static sqlite3_stmt * stmt_OutSWC;

static void beginTransaction(void);
static void endTransaction(void);
static void prepareStatements(void);
static void finalizeStatements(void);

static void beginTransaction() {
	char * sErrMsg = 0;
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);
}

static void endTransaction() {
	char * sErrMsg = 0;
	sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &sErrMsg);
}

static int callback(void *NotUsed, int argc, char **argv, char **azColName) {
	int i;
	for (i = 0; i < argc; i++) {
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	return 0;
}

static void sqlcheck(int rc, char *zErrMsg) {
	if (rc != SQLITE_OK) {
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	} else {
		//fprintf(stdout, "Table created successfully\n");
	}
}

void connect(char *debugout) {
	int rc;
	char name[256] = { 0 };
	char *zErrMsg = 0;
	strcat(name, debugout);
	strcat(name, ".sqlite3");

	remove(name);
	rc = sqlite3_open(name, &db);
	if (rc) {
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return;
	}

	sqlite3_exec(db, "PRAGMA synchronous = OFF", NULL, NULL, &zErrMsg);
	sqlite3_exec(db, "PRAGMA journal_mode = MEMORY", NULL, NULL, &zErrMsg);
}

void disconnect() {
	finalizeStatements();
	sqlite3_close(db);
}

static void insertRgroups(void) {
	int rc;
	GrpIndex g;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	ForEachGroup(g)
	{
		sql[0] = 0;
		sprintf(sql, "INSERT INTO RGroups (ID, NAME, VegProdType) VALUES (%d, '%s', %d);", g+1, RGroup[g]->name, RGroup[g]->veg_prod_type);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
		sqlcheck(rc, zErrMsg);
	}
	endTransaction();
}

void insertSXWPhen(void) {
	int rc,m;
	GrpIndex g;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	ForEachGroup(g)
	{
		for(m=0;m<12;m++) {
			sql[0] = 0;
			sprintf(sql, "INSERT INTO sxwphen (RGroupID, Month, GrowthPCT) VALUES (%d, %d, %f);", g+1, m+1,_phen[Igp(g,m)]);
			rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
			sqlcheck(rc, zErrMsg);
		}
	}
	endTransaction();
}

void insertSXWProd(void) {
	int rc,m;
	GrpIndex g;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	for(m=0;m<12;m++) {
		sql[0] = 0;
		sprintf(sql, "INSERT INTO sxwprod (Month, BMASS, LITTER) VALUES (%d, %f, %f);", m+1, _prod_conv[m][PC_Bmass], _prod_conv[m][PC_Litter]);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
		sqlcheck(rc, zErrMsg);
	}
	endTransaction();
}

void insertInfo() {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql, "INSERT INTO info (StartYear, Years, Iterations, RGroups, TranspirationLayers, SoilLayers, PlotSize) VALUES (%d, %d, %d, %d, %d, %d, %f);", SW_Model.startyr, Globals.runModelYears, Globals.runModelIterations, Globals.grpCount, SXW.NTrLyrs, SXW.NSoLyrs, Globals.plotsize);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void insertRootsXphenRow(GrpIndex g, int Layer, double Jan, double Feb, double March, double April, double May, double June, double July, double Aug, double Sept, double Oct, double November, double December) {
	int rc;
	char *zErrMsg = 0;

	sql[0] = 0;

	sprintf(sql, "INSERT INTO rootsXphen (RGroupID,Layer,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", g, Layer, Jan, Feb, March, April, May, June, July, Aug, Sept, Oct, November, December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
}

void insertRootsXphen(double * _rootsXphen) {
	int r, l, p;
	int nLyrs;
	double m[12];

	beginTransaction();
	ForEachGroup(r)
	{
		nLyrs = getNTranspLayers(RGroup[r]->veg_prod_type);
		for (l = 0; l < nLyrs; l++) {
			for(p=0;p<12;p++)
			{
				m[p] = _rootsXphen[Iglp(r, l, p)];
			}
			insertRootsXphenRow(r+1, l+1, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
		}
	}
	endTransaction();
}

static void insertSXWinputVarsRow(int year, int iter, double fracGrass, double fracShrub, double fracTree, double fracForb, double fracBareGround) {
	//int rc;
	//char *zErrMsg = 0;
	//sql[0] = 0;

	sqlite3_bind_int(stmt_InputVars, 1, year);
	sqlite3_bind_int(stmt_InputVars, 2, iter);
	sqlite3_bind_double(stmt_InputVars, 3, fracGrass);
	sqlite3_bind_double(stmt_InputVars, 4, fracShrub);
	sqlite3_bind_double(stmt_InputVars, 5, fracTree);
	sqlite3_bind_double(stmt_InputVars, 6, fracForb);
	sqlite3_bind_double(stmt_InputVars, 7, fracBareGround);

	sqlite3_step(stmt_InputVars);
	sqlite3_clear_bindings(stmt_InputVars);
	sqlite3_reset(stmt_InputVars);

	//sprintf(sql, "INSERT INTO sxwInputVars (Year,Iteration,FracGrass,FracShrub,FracTree,FracForb,FracBareGround) VALUES (%d, %d, %.4f, %.4f, %.4f, %.4f, %.4f);", year, iter, fracGrass, fracShrub, fracTree, fracForb, fracBareGround);
	//rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	//sqlcheck(rc, zErrMsg);
}

void insertInputVars() {
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;
	SW_VEGPROD *v = &SW_VegProd;

	beginTransaction();
	insertSXWinputVarsRow(Year, Iteration, v->fractionGrass, v->fractionShrub, v->fractionTree, v->fractionForb, v->fractionBareGround);
	endTransaction();
}

static void insertSXWinputProdRow(int year, int iter, int VegProdType, int Month, double Litter, double Biomass, double PLive, double LAI_conv) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwInputProd (Year,Iteration,VegProdType,Month,Litter,Biomass,PLive,LAI_conv) VALUES (%d, %d, %d, %d, %.3f, %.3f, %.3f, %.3f);", year, iter, VegProdType, Month, Litter, Biomass, PLive, LAI_conv);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_InputProd, 1, year);
	sqlite3_bind_int(stmt_InputProd, 2, iter);
	sqlite3_bind_int(stmt_InputProd, 3, VegProdType);
	sqlite3_bind_int(stmt_InputProd, 4, Month);
	sqlite3_bind_double(stmt_InputProd, 5, Litter);
	sqlite3_bind_double(stmt_InputProd, 6, Biomass);
	sqlite3_bind_double(stmt_InputProd, 7, PLive);
	sqlite3_bind_double(stmt_InputProd, 8, LAI_conv);

	sqlite3_step(stmt_InputProd);
	sqlite3_clear_bindings(stmt_InputProd);
	sqlite3_reset(stmt_InputProd);
}

void insertInputProd() {
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;
	int p;
	SW_VEGPROD *v = &SW_VegProd;

	beginTransaction();
	ForEachTrPeriod(p) {
		insertSXWinputProdRow(Year, Iteration, 1, p+1, v->tree.litter[p], v->tree.biomass[p], v->tree.pct_live[p], v->tree.lai_conv[p]);
		insertSXWinputProdRow(Year, Iteration, 2, p+1, v->shrub.litter[p], v->shrub.biomass[p], v->shrub.pct_live[p], v->shrub.lai_conv[p]);
		insertSXWinputProdRow(Year, Iteration, 3, p+1, v->grass.litter[p], v->grass.biomass[p], v->grass.pct_live[p], v->grass.lai_conv[p]);
		insertSXWinputProdRow(Year, Iteration, 4, p+1, v->forb.litter[p], v->forb.biomass[p], v->forb.pct_live[p], v->forb.lai_conv[p]);
	}
	endTransaction();
}

static void insertSXWinputSoilsRow(int year, int iter, int Layer, double Tree_trco, double Shrub_trco, double Grass_trco, double Forb_trco) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwInputSoils (Year,Iteration,Layer,Tree_trco,Shrub_trco,Grass_trco,Forb_trco) VALUES (%d, %d, %d, %.3f, %.3f, %.3f, %.3f);", year, iter, Layer, Tree_trco,Shrub_trco,Grass_trco,Forb_trco);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_InputSoils, 1, year);
	sqlite3_bind_int(stmt_InputSoils, 2, iter);
	sqlite3_bind_int(stmt_InputSoils, 3, Layer);
	sqlite3_bind_double(stmt_InputSoils, 4, Tree_trco);
	sqlite3_bind_double(stmt_InputSoils, 5, Shrub_trco);
	sqlite3_bind_double(stmt_InputSoils, 6, Grass_trco);
	sqlite3_bind_double(stmt_InputSoils, 7, Forb_trco);

	sqlite3_step(stmt_InputSoils);
	sqlite3_clear_bindings(stmt_InputSoils);
	sqlite3_reset(stmt_InputSoils);
}

void insertInputSoils() {
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;
	int l;
	SW_SITE *s = &SW_Site;

	beginTransaction();
	ForEachSoilLayer(l)
	{
		insertSXWinputSoilsRow(Year, Iteration, l+1, s->lyr[l]->transp_coeff_tree, s->lyr[l]->transp_coeff_shrub, s->lyr[l]->transp_coeff_grass, s->lyr[l]->transp_coeff_forb);
	}
	endTransaction();
}

static void insertSXWoutputVarsRow(int year, int iter, int MAP_mm, double MAT_C, double AET_cm, double AT_cm, double TotalRelsize, double TotalPR, double TotalTransp) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputVars (Year,Iteration,MAP_mm,MAT_C,AET_cm,AT_cm,TotalRelsize,TotalPR,TotalTransp) VALUES (%d, %d, %d, %5.2f, %5.4f, %5.4f, %.4f, %.4f, %.4f);", year, iter, MAP_mm, MAT_C, AET_cm, AT_cm, TotalRelsize, TotalPR, TotalTransp);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutVars, 1, year);
	sqlite3_bind_int(stmt_OutVars, 2, iter);
	sqlite3_bind_int(stmt_OutVars, 3, MAP_mm);
	sqlite3_bind_double(stmt_OutVars, 4, MAT_C);
	sqlite3_bind_double(stmt_OutVars, 5, AET_cm);
	sqlite3_bind_double(stmt_OutVars, 6, AT_cm);
	sqlite3_bind_double(stmt_OutVars, 7, TotalRelsize);
	sqlite3_bind_double(stmt_OutVars, 8, TotalPR);
	sqlite3_bind_double(stmt_OutVars, 9, TotalTransp);

	sqlite3_step(stmt_OutVars);
	sqlite3_clear_bindings(stmt_OutVars);
	sqlite3_reset(stmt_OutVars);
}

void insertOutputVars(RealF * _resource_cur) {
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;
	int p;
	int t;
	int r;
	double sum = 0;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	ForEachTrPeriod(p)
	{
		//ForEachTranspLayer(t) sum += SXW.transp[Ilp(t,p)];
		for (t = 0; t < SXW.NSoLyrs; t++)
			sum += SXW.transpTotal[Ilp(t, p)];
	}

	ForEachGroup(r) {
			sum1 += RGroup[r]->relsize;
			sum2 += RGroup[r]->pr;
			sum3 += _resource_cur[r];
	}
	beginTransaction();
	insertSXWoutputVarsRow(Year, Iteration, Env.ppt, Env.temp, SXW.aet, sum, sum1,sum2,sum3);
	endTransaction();
}

static void insertSXWoutputRgroupRow(int year, int iter, int RGroupID, double Biomass, double Realsize, double PR, double Transpiration) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputRgroup (Year,Iteration,RGroupID,Realsize,PR,Transpiration) VALUES (%d, %d, %d, %.4f, %.4f, %.4f);", year, iter, RGroupID, Realsize, PR, Transpiration);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutRgroup, 1, year);
	sqlite3_bind_int(stmt_OutRgroup, 2, iter);
	sqlite3_bind_int(stmt_OutRgroup, 3, RGroupID);
	sqlite3_bind_double(stmt_OutRgroup, 4, Biomass);
	sqlite3_bind_double(stmt_OutRgroup, 5, Realsize);
	sqlite3_bind_double(stmt_OutRgroup, 6, PR);
	sqlite3_bind_double(stmt_OutRgroup, 7, Transpiration);

	sqlite3_step(stmt_OutRgroup);
	sqlite3_clear_bindings(stmt_OutRgroup);
	sqlite3_reset(stmt_OutRgroup);
}

void insertRgroupInfo(RealF * _resource_cur) {
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;
	int r;

	beginTransaction();
	ForEachGroup(r) {
		insertSXWoutputRgroupRow(Year, Iteration, r+1, RGroup_GetBiomass(r),RGroup[r]->relsize, RGroup[r]->pr, _resource_cur[r]);
	}
	endTransaction();
}

static void insertSXWoutputProdRow(int year, int iter, int Month, double BMass, double PctLive, double LAIlive, double VegCov, double TotAGB) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputProd (Year,Iteration,Month,BMass,PctLive,LAIlive,VegCov,TotAGB) VALUES (%d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f);", year, iter,Month,BMass,PctLive,LAIlive,VegCov,TotAGB);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
}

void insertOutputProd(SW_VEGPROD *v) {
	int p;
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;

	beginTransaction();
	int doy = 1;
	ForEachMonth(p)
	{
		int days = 31, i;
		double pct_live = 0, lai_live = 0, vegcov = 0, total_agb = 0, biomass = 0;

		if (p == Apr || p == Jun || p == Sep || p == Nov) //all these months have 30 days
			days = 30;
		else if (p == Feb) { //February has either 28 or 29 days
			days = 28;
			if (Is_LeapYear(SW_Model.year))
				days = 29;
		} // all the other months have 31 days

		for (i = doy; i < (doy + days); i++) { //accumulating the monthly values...
			lai_live += (v->tree.lai_live_daily[doy])
					+ (v->shrub.lai_live_daily[doy])
					+ (v->grass.lai_live_daily[doy]);
			vegcov += (v->tree.vegcov_daily[doy]) + (v->shrub.vegcov_daily[doy])
					+ (v->grass.vegcov_daily[doy]);
			total_agb += (v->tree.total_agb_daily[doy])
					+ (v->shrub.total_agb_daily[doy])
					+ (v->grass.total_agb_daily[doy]);
			pct_live += (v->tree.pct_live_daily[doy])
					+ (v->shrub.pct_live_daily[doy])
					+ (v->grass.pct_live_daily[doy])
					+ (v->forb.pct_live_daily[doy]);
			biomass += (v->tree.biomass_daily[doy])
					+ (v->shrub.biomass_daily[doy])
					+ (v->grass.biomass_daily[doy])
					+ (v->forb.biomass_daily[doy]);
		}
		doy += days; //updating the doy

		pct_live /= days;
		biomass /= days;
		lai_live /= days; //getting the monthly averages...
		vegcov /= days;
		total_agb /= days;

		// had to update this commented out code because the VegProds are now daily values instead of monthly ones...
		//fprintf(f,"%4d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		//      p+1,v->biomass[p], v->pct_live[p],v->lai_live[p],
		//    v->vegcov[p], v->total_agb[p]);
		insertSXWoutputProdRow(Year,Iteration,p+1,biomass,pct_live, lai_live, vegcov, total_agb);
	}
	endTransaction();
}

static void insertSXWoutputRootsSumRow(int Year, int Iteration, int Layer, int VegProdType, double January, double February, double March, double April, double May, double June, double July, double August, double September, double October, double November, double December) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputRootsSum (Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutRootsSum, 1, Year);
	sqlite3_bind_int(stmt_OutRootsSum, 2, Iteration);
	sqlite3_bind_int(stmt_OutRootsSum, 3, Layer);
	sqlite3_bind_int(stmt_OutRootsSum, 4, VegProdType);
	sqlite3_bind_double(stmt_OutRootsSum, 5, January);
	sqlite3_bind_double(stmt_OutRootsSum, 6, February);
	sqlite3_bind_double(stmt_OutRootsSum, 7, March);
	sqlite3_bind_double(stmt_OutRootsSum, 8, April);
	sqlite3_bind_double(stmt_OutRootsSum, 9, May);
	sqlite3_bind_double(stmt_OutRootsSum, 10, June);
	sqlite3_bind_double(stmt_OutRootsSum, 11, July);
	sqlite3_bind_double(stmt_OutRootsSum, 12, August);
	sqlite3_bind_double(stmt_OutRootsSum, 13, September);
	sqlite3_bind_double(stmt_OutRootsSum, 14, October);
	sqlite3_bind_double(stmt_OutRootsSum, 15, November);
	sqlite3_bind_double(stmt_OutRootsSum, 16, December);

	sqlite3_step(stmt_OutRootsSum);
	sqlite3_clear_bindings(stmt_OutRootsSum);
	sqlite3_reset(stmt_OutRootsSum);
}

void insertRootsSum(RealD * _roots_active_sum) {
	int l;
	int p;
	int i;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;

	beginTransaction();
	for(i=1; i<=4; i++) {
		for (l = 0; l < SXW.NSoLyrs; l++) {
			for (p = 0; p < 12; p++) {
				m[p] = _roots_active_sum[Itlp(i-1,l, p)];
			}
			insertSXWoutputRootsSumRow(Year, Iteration, l+1, i+1, m[0], m[1],
					m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10],
					m[11]);
		}
	}
	endTransaction();
}

static void insertSXWoutputRootsRelativeRow(int Year, int Iteration, int Layer, int RGroupID, double January, double February, double March, double April, double May, double June, double July, double August, double September, double October, double November, double December) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputRootsRelative (Year,Iteration,Layer,RGroupID,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", Year,Iteration,Layer,RGroupID,January,February,March,April,May,June,July,August,September,October,November,December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutRootsRel, 1, Year);
	sqlite3_bind_int(stmt_OutRootsRel, 2, Iteration);
	sqlite3_bind_int(stmt_OutRootsRel, 3, Layer);
	sqlite3_bind_int(stmt_OutRootsRel, 4, RGroupID);
	sqlite3_bind_double(stmt_OutRootsRel, 5, January);
	sqlite3_bind_double(stmt_OutRootsRel, 6, February);
	sqlite3_bind_double(stmt_OutRootsRel, 7, March);
	sqlite3_bind_double(stmt_OutRootsRel, 8, April);
	sqlite3_bind_double(stmt_OutRootsRel, 9, May);
	sqlite3_bind_double(stmt_OutRootsRel, 10, June);
	sqlite3_bind_double(stmt_OutRootsRel, 11, July);
	sqlite3_bind_double(stmt_OutRootsRel, 12, August);
	sqlite3_bind_double(stmt_OutRootsRel, 13, September);
	sqlite3_bind_double(stmt_OutRootsRel, 14, October);
	sqlite3_bind_double(stmt_OutRootsRel, 15, November);
	sqlite3_bind_double(stmt_OutRootsRel, 16, December);

	sqlite3_step(stmt_OutRootsRel);
	sqlite3_clear_bindings(stmt_OutRootsRel);
	sqlite3_reset(stmt_OutRootsRel);
}

void insertRootsRelative(RealD * _roots_active_rel) {
	int l;
	int p;
	int r;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;

	beginTransaction();
	ForEachGroup(r)
	{
		for (l = 0; l < SXW.NSoLyrs; l++) {
			for (p = 0; p < 12; p++) {
				m[p] = _roots_active_rel[Iglp(r, l, p)];
			}
			insertSXWoutputRootsRelativeRow(Year, Iteration, l+1, r+1, m[0], m[1],
					m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10],
					m[11]);
		}
	}
	endTransaction();
}

static void insertSXWoutputTranspirationRow(int Year, int Iteration, int Layer, int VegProdType, double January, double February, double March, double April, double May, double June, double July, double August, double September, double October, double November, double December) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputTranspiration (Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutTransp, 1, Year);
	sqlite3_bind_int(stmt_OutTransp, 2, Iteration);
	sqlite3_bind_int(stmt_OutTransp, 3, Layer);
	sqlite3_bind_int(stmt_OutTransp, 4, VegProdType);
	sqlite3_bind_double(stmt_OutTransp, 5, January);
	sqlite3_bind_double(stmt_OutTransp, 6, February);
	sqlite3_bind_double(stmt_OutTransp, 7, March);
	sqlite3_bind_double(stmt_OutTransp, 8, April);
	sqlite3_bind_double(stmt_OutTransp, 9, May);
	sqlite3_bind_double(stmt_OutTransp, 10, June);
	sqlite3_bind_double(stmt_OutTransp, 11, July);
	sqlite3_bind_double(stmt_OutTransp, 12, August);
	sqlite3_bind_double(stmt_OutTransp, 13, September);
	sqlite3_bind_double(stmt_OutTransp, 14, October);
	sqlite3_bind_double(stmt_OutTransp, 15, November);
	sqlite3_bind_double(stmt_OutTransp, 16, December);

	sqlite3_step(stmt_OutTransp);
	sqlite3_clear_bindings(stmt_OutTransp);
	sqlite3_reset(stmt_OutTransp);
}

void insertTranspiration() {
	int l;
	int p;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;

	beginTransaction();
	//Total - 0
	for (l = 0; l < SXW.NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW.transpTotal[Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year,Iteration,l+1,0,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}

	//Tree - 1
	for (l = 0; l < SXW.NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW.transpTrees[Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year,Iteration,l+1,1,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}

	//Shrub - 2
	for (l = 0; l < SXW.NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW.transpShrubs[Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year,Iteration,l+1,2,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}

	//Grass - 3
	for (l = 0; l < SXW.NSoLyrs; l++) {
		for (p = 0; p < 12; p++) {
			m[p] = SXW.transpGrasses[Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year, Iteration, l+1, 3, m[0], m[1], m[2],
				m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
	}

	//Forb - 4
	for (l = 0; l < SXW.NSoLyrs; l++) {
		for (p = 0; p < 12; p++) {
			m[p] = SXW.transpForbs[Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year, Iteration, l+1, 4, m[0], m[1], m[2],
				m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
	}
	endTransaction();
}

static void insertSXWoutputSWCBulkRow(int Year, int Iteration, int Layer, double January, double February, double March, double April, double May, double June, double July, double August, double September, double October, double November, double December) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputSWCBulk (Year,Iteration,Layer,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", Year,Iteration,Layer,January,February,March,April,May,June,July,August,September,October,November,December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutSWC, 1, Year);
	sqlite3_bind_int(stmt_OutSWC, 2, Iteration);
	sqlite3_bind_int(stmt_OutSWC, 3, Layer);
	sqlite3_bind_double(stmt_OutSWC, 4, January);
	sqlite3_bind_double(stmt_OutSWC, 5, February);
	sqlite3_bind_double(stmt_OutSWC, 6, March);
	sqlite3_bind_double(stmt_OutSWC, 7, April);
	sqlite3_bind_double(stmt_OutSWC, 8, May);
	sqlite3_bind_double(stmt_OutSWC, 9, June);
	sqlite3_bind_double(stmt_OutSWC, 10, July);
	sqlite3_bind_double(stmt_OutSWC, 11, August);
	sqlite3_bind_double(stmt_OutSWC, 12, September);
	sqlite3_bind_double(stmt_OutSWC, 13, October);
	sqlite3_bind_double(stmt_OutSWC, 14, November);
	sqlite3_bind_double(stmt_OutSWC, 15, December);

	sqlite3_step(stmt_OutSWC);
	sqlite3_clear_bindings(stmt_OutSWC);
	sqlite3_reset(stmt_OutSWC);
}

void insertSWCBulk() {
	int l;
	int p;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals.currIter;

	beginTransaction();
	for (l = 0; l < SXW.NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW.swc[Ilp(l, p)];
		}
		insertSXWoutputSWCBulkRow(Year,Iteration,l+1,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}
	endTransaction();
}

static void prepareStatements() {
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwInputVars (Year,Iteration,FracGrass,FracShrub,FracTree,FracForb,FracBareGround) VALUES (@Year,@Iteration,@FracGrass,@FracShrub,@FracTree,@FracForb,@FracBareGround);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_InputVars, NULL);

	sprintf(sql, "INSERT INTO sxwInputProd (Year,Iteration,VegProdType,Month,Litter,Biomass,PLive,LAI_conv) VALUES (@Year,@Iteration,@VegProdType,@Month,@Litter,@Biomass,@PLive,@LAI_conv);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_InputProd, NULL);

	sprintf(sql, "INSERT INTO sxwInputSoils (Year,Iteration,Layer,Tree_trco,Shrub_trco,Grass_trco,Forb_trco) VALUES (@Year,@Iteration,@Layer,@Tree_trco,@Shrub_trco,@Grass_trco,@Forb_trco);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_InputSoils, NULL);

	sprintf(sql, "INSERT INTO sxwOutputVars (Year,Iteration,MAP_mm,MAT_C,AET_cm,AT_cm,TotalRelsize,TotalPR,TotalTransp) VALUES (@Year,@Iteration,@MAP_mm,@MAT_C,@AET_cm,@AT_cm,@TotalRelsize,@TotalPR,@TotalTransp);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutVars, NULL);

	sprintf(sql, "INSERT INTO sxwOutputRgroup (Year,Iteration,RGroupID,Biomass,Realsize,PR,Transpiration) VALUES (@Year,@Iteration,@RGroupID,@Biomass,@Realsize,@PR,@Transpiration);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutRgroup, NULL);

	sprintf(sql, "INSERT INTO sxwOutputProd (Year,Iteration,Month,BMass,PctLive,LAIlive,VegCov,TotAGB) VALUES (@Year,@Iteration,@Month,@BMass,@PctLive,@LAIlive,@VegCov,@TotAGB);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutProd, NULL);

	sprintf(sql, "INSERT INTO sxwOutputRootsSum (Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (@Year,@Iteration,@Layer,@VegProdType,@January,@February,@March,@April,@May,@June,@July,@August,@September,@October,@November,@December);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutRootsSum, NULL);

	sprintf(sql, "INSERT INTO sxwOutputRootsRelative (Year,Iteration,Layer,RGroupID,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (@Year,@Iteration,@Layer,@RGroupID,@January,@February,@March,@April,@May,@June,@July,@August,@September,@October,@November,@December);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutRootsRel, NULL);

	sprintf(sql, "INSERT INTO sxwOutputTranspiration (Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (@Year,@Iteration,@Layer,@VegProdType,@January,@February,@March,@April,@May,@June,@July,@August,@September,@October,@November,@December);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutTransp, NULL);

	sprintf(sql, "INSERT INTO sxwOutputSWCBulk (Year,Iteration,Layer,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (@Year,@Iteration,@Layer,@January,@February,@March,@April,@May,@June,@July,@August,@September,@October,@November,@December);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutSWC, NULL);
}

static void finalizeStatements() {
	sqlite3_finalize(stmt_InputVars);
	sqlite3_finalize(stmt_InputProd);
	sqlite3_finalize(stmt_InputSoils);
	sqlite3_finalize(stmt_OutVars);
	sqlite3_finalize(stmt_OutRgroup);
	sqlite3_finalize(stmt_OutProd);
	sqlite3_finalize(stmt_OutRootsSum);
	sqlite3_finalize(stmt_OutRootsRel);
	sqlite3_finalize(stmt_OutTransp);
	sqlite3_finalize(stmt_OutSWC);
}

void createTables() {
	int rc;
	char *zErrMsg = 0;

	char *table_PrjInfo = "CREATE TABLE info(StartYear INT, Years INT, Iterations INT, RGroups INT, TranspirationLayers INT, SoilLayers INT, PlotSize REAL);";
	char *table_rgroups = "CREATE TABLE RGroups(ID INT PRIMARY KEY NOT NULL, NAME TEXT NOT NULL, VegProdType INT NOT NULL);";
	char *table_phen = "CREATE TABLE sxwphen(RGroupID INT NOT NULL, Month INT NOT NULL, GrowthPCT REAL NOT NULL, PRIMARY KEY(RGroupID, Month));";
	char *table_prod = "CREATE TABLE sxwprod(Month INT NOT NULL, BMASS REAL, LITTER REAL, PRIMARY KEY(Month));";
	char *table_rootsXphen = "CREATE TABLE rootsXphen(RGroupID INT NOT NULL, Layer INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(RGroupID, Layer));";
	char *table_InputVars =
			"CREATE TABLE sxwInputVars(Year INT NOT NULL, Iteration INT NOT NULL, FracGrass REAL, FracShrub REAL, FracTree REAL, FracForb REAL, FracBareGround REAL, PRIMARY KEY(Year, Iteration));";
	char *table_InputProd =
			"CREATE TABLE sxwInputProd(Year INT NOT NULL, Iteration INT NOT NULL, VegProdType INT NOT NULL, Month INT NOT NULL, Litter REAL, Biomass REAL, PLive REAL, LAI_conv REAL, PRIMARY KEY(Year, Iteration, VegProdType, Month));";
	char *table_InputSoils =
			"CREATE TABLE sxwInputSoils(Year INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, Tree_trco REAL, Shrub_trco REAL, Grass_trco REAL, Forb_trco REAL, PRIMARY KEY(Year, Iteration, Layer));";
	char *table_OutputVars =
			"CREATE TABLE sxwOutputVars(Year INT NOT NULL, Iteration INT NOT NULL, MAP_mm INT, MAT_C REAL, AET_cm REAL, AT_cm REAL, TotalRelsize REAL, TotalPR REAL, TotalTransp REAL, PRIMARY KEY(Year, Iteration));";
	char *table_OutputRgroup =
			"CREATE TABLE sxwOutputRgroup(YEAR INT NOT NULL, Iteration INT NOT NULL, RGroupID INT NOT NULL, Biomass REAL, Realsize REAL, PR REAL, Transpiration REAL, PRIMARY KEY(Year, Iteration, RGroupID));";
	char *table_OutputProd =
			"CREATE TABLE sxwOutputProd(YEAR INT NOT NULL, Iteration INT NOT NULL, Month INT NOT NULL, BMass REAL, PctLive REAL, LAIlive REAL, VegCov REAL, TotAGB REAL, PRIMARY KEY(Year, Iteration, Month));";
	char *table_OutputRootsSum =
			"CREATE TABLE sxwOutputRootsSum(YEAR INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, VegProdType INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(Year, Iteration, Layer, VegProdType));";
	char *table_OutputRootsRelative =
			"CREATE TABLE sxwOutputRootsRelative(YEAR INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, RGroupID INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(Year, Iteration, Layer, RGroupID));";
	char *table_OutputTransp =
			"CREATE TABLE sxwOutputTranspiration(YEAR INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, VegProdType INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(Year, Iteration, Layer, VegProdType));";
	char *table_OutputSWCBulk =
			"CREATE TABLE sxwOutputSWCBulk(YEAR INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(Year, Iteration, Layer));";

	rc = sqlite3_exec(db, table_PrjInfo, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_rgroups, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_phen, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_prod, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_rootsXphen, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_InputVars, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_InputProd, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_InputSoils, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputVars, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputRgroup, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputProd, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputRootsSum, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputRootsRelative, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputTransp, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputSWCBulk, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	insertRgroups();
	prepareStatements();
}
