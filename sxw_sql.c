/**
 * \file sxw_sql.c
 * \brief Outputs \ref SXW information to an SQL database.
 * 
 * \author Ryan J. Murphy
 * \date 15 January 2015
 * \ingroup SQL
 * \ingroup SXW_PRIVATE
 */

#include <string.h>
#include <sqlite3.h>
#include "ST_steppe.h"
#include "ST_globals.h"
#include "sxw_module.h"
#include "sxw.h"

extern SW_MODEL SW_Model;
extern SXW_t* SXW;
extern SW_SITE SW_Site;
extern SW_VEGPROD SW_VegProd;
extern SXW_resourceType* SXWResources;

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
			sprintf(sql, "INSERT INTO sxwphen (RGroupID, Month, GrowthPCT) VALUES (%d, %d, %f);", g+1, m+1,SXWResources->_phen[Igp(g,m)]);
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
	ForEachGroup(g)
	{
	for(m=0;m<12;m++) {
		sql[0] = 0;
		sprintf(sql, "INSERT INTO sxwprod (RGroupID, Month, BMASS, LITTER, PCTLIVE) VALUES (%d, %d, %f, %f, %f);", 
					  g+1, m+1, SXWResources->_prod_bmass[Igp(g,m)], SXWResources->_prod_litter[g][m], SXWResources->_prod_pctlive[Igp(g,m)]);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
		sqlcheck(rc, zErrMsg);
	}
	}
	endTransaction();
}

void insertInfo() {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql, "INSERT INTO info (StartYear, Years, Iterations, RGroups, TranspirationLayers, SoilLayers, PlotSize) VALUES (%d, %d, %d, %d, %d, %d, %f);", 
				  SW_Model.startyr, SuperGlobals.runModelYears, SuperGlobals.runModelIterations, Globals->grpCount, SXW->NTrLyrs, SXW->NSoLyrs, Globals->plotsize);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void insertRootsXphenRow(GrpIndex g, int Layer, double var_Jan, double var_Feb, double var_Mar, double var_Apr, double var_May, double var_Jun, double var_Jul, double var_Aug, double var_Sep, double var_Oct, double var_Nov, double var_Dec) {
	int rc;
	char *zErrMsg = 0;

	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwRootsXphen (RGroupID,Layer,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", g, Layer, var_Jan, var_Feb, var_Mar, var_Apr, var_May, var_Jun, var_Jul, var_Aug, var_Sep, var_Oct, var_Nov, var_Dec);
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
	int Iteration = Globals->currIter;
	SW_VEGPROD *v = &SW_VegProd;

	beginTransaction();
	insertSXWinputVarsRow(Year, Iteration, v->veg[3].cov.fCover, v->veg[1].cov.fCover, v->veg[0].cov.fCover, v->veg[2].cov.fCover, v->bare_cov.fCover);
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
	int Iteration = Globals->currIter;
	int p;
	SW_VEGPROD *v = &SW_VegProd;

	beginTransaction();
	ForEachTrPeriod(p) {
		insertSXWinputProdRow(Year, Iteration, 1, p+1, v->veg[0].litter[p], v->veg[0].biomass[p], v->veg[0].pct_live[p], v->veg[0].lai_conv[p]);
		insertSXWinputProdRow(Year, Iteration, 2, p+1, v->veg[1].litter[p], v->veg[1].biomass[p], v->veg[1].pct_live[p], v->veg[1].lai_conv[p]);
		insertSXWinputProdRow(Year, Iteration, 3, p+1, v->veg[3].litter[p], v->veg[3].biomass[p], v->veg[3].pct_live[p], v->veg[3].lai_conv[p]);
		insertSXWinputProdRow(Year, Iteration, 4, p+1, v->veg[2].litter[p], v->veg[2].biomass[p], v->veg[2].pct_live[p], v->veg[2].lai_conv[p]);
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
	int Iteration = Globals->currIter;
	int l;
	SW_SITE *s = &SW_Site;

	beginTransaction();
	ForEachSoilLayer(l)
	{
		insertSXWinputSoilsRow(Year, Iteration, l+1, s->lyr[l]->transp_coeff[0], s->lyr[l]->transp_coeff[1], s->lyr[l]->transp_coeff[3], s->lyr[l]->transp_coeff[2]);
	}
	endTransaction();
}

static void insertSXWoutputVarsRow(int year, int iter, int MAP_mm, double MAT_C, double AET_cm, double T_cm, double ADT_cm, double AT_cm, double TotalRelsize, double TotalPR, double TotalTransp) {
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
	sqlite3_bind_double(stmt_OutVars, 6, T_cm);
	sqlite3_bind_double(stmt_OutVars, 7, ADT_cm);
	sqlite3_bind_double(stmt_OutVars, 8, AT_cm);
	sqlite3_bind_double(stmt_OutVars, 9, TotalRelsize);
	sqlite3_bind_double(stmt_OutVars, 10, TotalPR);
	sqlite3_bind_double(stmt_OutVars, 11, TotalTransp);

	sqlite3_step(stmt_OutVars);
	sqlite3_clear_bindings(stmt_OutVars);
	sqlite3_reset(stmt_OutVars);
}

void insertOutputVars(RealF * _resource_cur, RealF added_transp) {
	int Year = SW_Model.year;
	int Iteration = Globals->currIter;
	int p;
	int t;
	int r;
	double sum = 0;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	ForEachTrPeriod(p)
	{
		for (t = 0; t < SXW->NSoLyrs; t++)
			sum += SXW->transpTotal[Ilp(t, p)];
	}

	ForEachGroup(r) {
			sum1 += getRGroupRelsize(r);
			sum2 += RGroup[r]->pr;
			sum3 += _resource_cur[r];
	}
	beginTransaction();
	insertSXWoutputVarsRow(Year, Iteration, Env->ppt, Env->temp, SXW->aet, sum, added_transp, sum+added_transp, sum1,sum2,sum3);
	endTransaction();
}

static void insertSXWoutputRgroupRow(int year, int iter, int RGroupID, double Biomass, double Realsize, double PR, double pre_bvt, double resource_cur) {
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
	sqlite3_bind_double(stmt_OutRgroup, 7, pre_bvt);
	sqlite3_bind_double(stmt_OutRgroup, 8, resource_cur);

	sqlite3_step(stmt_OutRgroup);
	sqlite3_clear_bindings(stmt_OutRgroup);
	sqlite3_reset(stmt_OutRgroup);
}

void insertRgroupInfo(RealF * _resource_cur) {
	int Year = SW_Model.year;
	int Iteration = Globals->currIter;
	int r;

	beginTransaction();
	ForEachGroup(r) {
		insertSXWoutputRgroupRow(Year, Iteration, r+1, RGroup_GetBiomass(r),getRGroupRelsize(r), RGroup[r]->pr, SXWResources->_resource_cur[r]/RGroup[r]->_bvt, _resource_cur[r]);
	}
	endTransaction();
}

static void insertSXWoutputProdRow(int year, int iter, int Month, double BMass, double PctLive, double LAIlive, double LAItotal, double TotAGB) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputProd (Year,Iteration,Month,BMass,PctLive,LAIlive,LAItotal,TotAGB) VALUES (%d, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f);", year, iter,Month,BMass,PctLive,LAIlive,LAItotal,TotAGB);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
}

void insertOutputProd(SW_VEGPROD *v) {
	int p;
	int Year = SW_Model.year;
	int Iteration = Globals->currIter;

	beginTransaction();
	int doy = 1;
	ForEachMonth(p)
	{
		int days = 31, i;
		double pct_live = 0, lai_live = 0, bLAI_total = 0, total_agb = 0, biomass = 0;

		if (p == Apr || p == Jun || p == Sep || p == Nov) //all these months have 30 days
			days = 30;
		else if (p == Feb) { //February has either 28 or 29 days
			days = 28;
			if (isleapyear(SW_Model.year))
				days = 29;
		} // all the other months have 31 days

		for (i = doy; i < (doy + days); i++) { //accumulating the monthly values...
			lai_live += (v->veg[0].lai_live_daily[i])
					+ (v->veg[1].lai_live_daily[i])
					+ (v->veg[3].lai_live_daily[i])
					+ (v->veg[2].lai_live_daily[i]);
			bLAI_total += (v->veg[0].bLAI_total_daily[i]) + (v->veg[1].bLAI_total_daily[i])
					+ (v->veg[3].bLAI_total_daily[i]) + (v->veg[2].bLAI_total_daily[i]);
			total_agb += (v->veg[0].total_agb_daily[i])
					+ (v->veg[1].total_agb_daily[i])
					+ (v->veg[3].total_agb_daily[i])
					+ (v->veg[2].total_agb_daily[i]);
			pct_live += (v->veg[0].pct_live_daily[i])
					+ (v->veg[1].pct_live_daily[i])
					+ (v->veg[3].pct_live_daily[i])
					+ (v->veg[2].pct_live_daily[i]);
			biomass += (v->veg[0].biomass_daily[i])
					+ (v->veg[1].biomass_daily[i])
					+ (v->veg[3].biomass_daily[i])
					+ (v->veg[2].biomass_daily[i]);
		}
		doy += days; //updating the doy

		pct_live /= days;
		biomass /= days;
		lai_live /= days; //getting the monthly averages...
		bLAI_total /= days;
		total_agb /= days;

		insertSXWoutputProdRow(Year,Iteration,p+1,biomass,pct_live, lai_live, bLAI_total, total_agb);
	}
	endTransaction();
}

static void insertSXWoutputRootsSumRow(int Year, int Iteration, int Layer, int VegProdType, double var_Jan, double var_Feb, double var_Mar, double var_Apr, double var_May, double var_Jun, double var_Jul, double var_Aug, double var_Sep, double var_Oct, double var_Nov, double var_Dec) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputRootsSum (Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", Year,Iteration,Layer,VegProdType,January,February,var_Mar,April,May,June,July,August,September,October,November,December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutRootsSum, 1, Year);
	sqlite3_bind_int(stmt_OutRootsSum, 2, Iteration);
	sqlite3_bind_int(stmt_OutRootsSum, 3, Layer);
	sqlite3_bind_int(stmt_OutRootsSum, 4, VegProdType);
	sqlite3_bind_double(stmt_OutRootsSum, 5, var_Jan);
	sqlite3_bind_double(stmt_OutRootsSum, 6, var_Feb);
	sqlite3_bind_double(stmt_OutRootsSum, 7, var_Mar);
	sqlite3_bind_double(stmt_OutRootsSum, 8, var_Apr);
	sqlite3_bind_double(stmt_OutRootsSum, 9, var_May);
	sqlite3_bind_double(stmt_OutRootsSum, 10, var_Jun);
	sqlite3_bind_double(stmt_OutRootsSum, 11, var_Jul);
	sqlite3_bind_double(stmt_OutRootsSum, 12, var_Aug);
	sqlite3_bind_double(stmt_OutRootsSum, 13, var_Sep);
	sqlite3_bind_double(stmt_OutRootsSum, 14, var_Oct);
	sqlite3_bind_double(stmt_OutRootsSum, 15, var_Nov);
	sqlite3_bind_double(stmt_OutRootsSum, 16, var_Dec);

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
	int Iteration = Globals->currIter;

	beginTransaction();

  ForEachVegType(i) {
		for (l = 0; l < SXW->NTrLyrs; l++) {
			for (p = 0; p < 12; p++) {
				m[p] = _roots_active_sum[Itlp(i, l, p)];
			}
			insertSXWoutputRootsSumRow(Year, Iteration, l+1, i, m[0], m[1],
					m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10],
					m[11]);
		}
	}
	endTransaction();
}

static void insertSXWoutputRootsRelativeRow(int Year, int Iteration, int Layer, int RGroupID, double var_Jan, double var_Feb, double var_Mar, double var_Apr, double var_May, double var_Jun, double var_Jul, double var_Aug, double var_Sep, double var_Oct, double var_Nov, double var_Dec) {
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
	sqlite3_bind_double(stmt_OutRootsRel, 5, var_Jan);
	sqlite3_bind_double(stmt_OutRootsRel, 6, var_Feb);
	sqlite3_bind_double(stmt_OutRootsRel, 7, var_Mar);
	sqlite3_bind_double(stmt_OutRootsRel, 8, var_Apr);
	sqlite3_bind_double(stmt_OutRootsRel, 9, var_May);
	sqlite3_bind_double(stmt_OutRootsRel, 10, var_Jun);
	sqlite3_bind_double(stmt_OutRootsRel, 11, var_Jul);
	sqlite3_bind_double(stmt_OutRootsRel, 12, var_Aug);
	sqlite3_bind_double(stmt_OutRootsRel, 13, var_Sep);
	sqlite3_bind_double(stmt_OutRootsRel, 14, var_Oct);
	sqlite3_bind_double(stmt_OutRootsRel, 15, var_Nov);
	sqlite3_bind_double(stmt_OutRootsRel, 16, var_Dec);

	sqlite3_step(stmt_OutRootsRel);
	sqlite3_clear_bindings(stmt_OutRootsRel);
	sqlite3_reset(stmt_OutRootsRel);
}

void insertRootsRelative(RealD * _roots_active_rel) {
	int l;
	int p;
	int g;
	int nLyrs;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals->currIter;

	beginTransaction();
	ForEachGroup(g)
	{
		nLyrs = getNTranspLayers(RGroup[g]->veg_prod_type);
		for (l = 0; l < nLyrs; l++) {
			for (p = 0; p < 12; p++) {
				m[p] = _roots_active_rel[Iglp(g, l, p)];
			}
			insertSXWoutputRootsRelativeRow(Year, Iteration, l+1, g+1, m[0], m[1],
					m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10],
					m[11]);
		}
	}
	endTransaction();
}

static void insertSXWoutputTranspirationRow(int Year, int Iteration, int Layer, int VegProdType, double var_Jan, double var_Feb, double var_Mar, double var_Apr, double var_May, double var_Jun, double var_Jul, double var_Aug, double var_Sep, double var_Oct, double var_Nov, double var_Dec) {
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
	sqlite3_bind_double(stmt_OutTransp, 5, var_Jan);
	sqlite3_bind_double(stmt_OutTransp, 6, var_Feb);
	sqlite3_bind_double(stmt_OutTransp, 7, var_Mar);
	sqlite3_bind_double(stmt_OutTransp, 8, var_Apr);
	sqlite3_bind_double(stmt_OutTransp, 9, var_May);
	sqlite3_bind_double(stmt_OutTransp, 10, var_Jun);
	sqlite3_bind_double(stmt_OutTransp, 11, var_Jul);
	sqlite3_bind_double(stmt_OutTransp, 12, var_Aug);
	sqlite3_bind_double(stmt_OutTransp, 13, var_Sep);
	sqlite3_bind_double(stmt_OutTransp, 14, var_Oct);
	sqlite3_bind_double(stmt_OutTransp, 15, var_Nov);
	sqlite3_bind_double(stmt_OutTransp, 16, var_Dec);

	sqlite3_step(stmt_OutTransp);
	sqlite3_clear_bindings(stmt_OutTransp);
	sqlite3_reset(stmt_OutTransp);
}

void insertTranspiration() {
	int l;
	int p;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals->currIter;

	beginTransaction();
	//Total - 0
	for (l = 0; l < SXW->NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW->transpTotal[Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year,Iteration,l+1,0,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}

	//Tree - 1
	for (l = 0; l < SXW->NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW->transpVeg[SW_TREES][Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year,Iteration,l+1,1,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}

	//Shrub - 2
	for (l = 0; l < SXW->NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW->transpVeg[SW_SHRUB][Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year,Iteration,l+1,2,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11]);
	}

	//Grass - 3
	for (l = 0; l < SXW->NSoLyrs; l++) {
		for (p = 0; p < 12; p++) {
			m[p] = SXW->transpVeg[SW_GRASS][Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year, Iteration, l+1, 3, m[0], m[1], m[2],
				m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
	}

	//Forb - 4
	for (l = 0; l < SXW->NSoLyrs; l++) {
		for (p = 0; p < 12; p++) {
			m[p] = SXW->transpVeg[SW_FORBS][Ilp(l, p)];
		}
		insertSXWoutputTranspirationRow(Year, Iteration, l+1, 4, m[0], m[1], m[2],
				m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
	}
	endTransaction();
}

static void insertSXWoutputSWCBulkRow(int Year, int Iteration, int Layer, double var_Jan, double var_Feb, double var_Mar, double var_Apr, double var_May, double var_Jun, double var_Jul, double var_Aug, double var_Sep, double var_Oct, double var_Nov, double var_Dec) {
	/*int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	sprintf(sql, "INSERT INTO sxwOutputSWCBulk (Year,Iteration,Layer,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (%d, %d, %d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f);", Year,Iteration,Layer,January,February,March,April,May,June,July,August,September,October,November,December);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);*/

	sqlite3_bind_int(stmt_OutSWC, 1, Year);
	sqlite3_bind_int(stmt_OutSWC, 2, Iteration);
	sqlite3_bind_int(stmt_OutSWC, 3, Layer);
	sqlite3_bind_double(stmt_OutSWC, 4, var_Jan);
	sqlite3_bind_double(stmt_OutSWC, 5, var_Feb);
	sqlite3_bind_double(stmt_OutSWC, 6, var_Mar);
	sqlite3_bind_double(stmt_OutSWC, 7, var_Apr);
	sqlite3_bind_double(stmt_OutSWC, 8, var_May);
	sqlite3_bind_double(stmt_OutSWC, 9, var_Jun);
	sqlite3_bind_double(stmt_OutSWC, 10, var_Jul);
	sqlite3_bind_double(stmt_OutSWC, 11, var_Aug);
	sqlite3_bind_double(stmt_OutSWC, 12, var_Sep);
	sqlite3_bind_double(stmt_OutSWC, 13, var_Oct);
	sqlite3_bind_double(stmt_OutSWC, 14, var_Nov);
	sqlite3_bind_double(stmt_OutSWC, 15, var_Dec);

	sqlite3_step(stmt_OutSWC);
	sqlite3_clear_bindings(stmt_OutSWC);
	sqlite3_reset(stmt_OutSWC);
}

void insertSWCBulk() {
	int l;
	int p;
	double m[12];
	int Year = SW_Model.year;
	int Iteration = Globals->currIter;

	beginTransaction();
	for (l = 0; l < SXW->NSoLyrs; l++) {
		for(p=0;p<12;p++) {
			m[p] = SXW->swc[Ilp(l, p)];
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

	sprintf(sql, "INSERT INTO sxwOutputVars (Year,Iteration,MAP_mm,MAT_C,AET_cm,T_cm,ADT_cm,AT_cm,TotalRelsize,TotalPR,TotalTransp) VALUES (@Year,@Iteration,@MAP_mm,@MAT_C,@AET_cm,@T_cm,@ADT_cm,@AT_cm,@TotalRelsize,@TotalPR,@TotalTransp);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutVars, NULL);

	sprintf(sql, "INSERT INTO sxwOutputRgroup (Year,Iteration,RGroupID,Biomass,Realsize,PR,resource_cur_preBvt,resource_cur) VALUES (@Year,@Iteration,@RGroupID,@Biomass,@Realsize,@PR,@pre_bvt,@resource_cur);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutRgroup, NULL);

	sprintf(sql, "INSERT INTO sxwOutputProd (Year,Iteration,Month,BMass,PctLive,LAIlive,LAItotal,TotAGB) VALUES (@Year,@Iteration,@Month,@BMass,@PctLive,@LAIlive,@LAItotal,@TotAGB);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutProd, NULL);

	sprintf(sql, "INSERT INTO sxwRootsSum (Year,Iteration,Layer,VegProdType,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (@Year,@Iteration,@Layer,@VegProdType,@January,@February,@March,@April,@May,@June,@July,@August,@September,@October,@November,@December);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_OutRootsSum, NULL);

	sprintf(sql, "INSERT INTO sxwRootsRelative (Year,Iteration,Layer,RGroupID,January,February,March,April,May,June,July,August,September,October,November,December) VALUES (@Year,@Iteration,@Layer,@RGroupID,@January,@February,@March,@April,@May,@June,@July,@August,@September,@October,@November,@December);");
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

	char *table_PrjInfo = "CREATE TABLE info(StartYear INT, Years INT, Iterations INT, RGroups INT, TranspirationLayers INT, SoilLayers INT, PlotSize REAL, BVT REAL);";
	char *table_rgroups = "CREATE TABLE RGroups(ID INT PRIMARY KEY NOT NULL, NAME TEXT NOT NULL, VegProdType INT NOT NULL);";
	char *table_phen = "CREATE TABLE sxwphen(RGroupID INT NOT NULL, Month INT NOT NULL, GrowthPCT REAL NOT NULL, PRIMARY KEY(RGroupID, Month));";
	char *table_prod = "CREATE TABLE sxwprod(RGroupID INT NOT NULL, Month INT NOT NULL, BMASS REAL, LITTER REAL, PCTLIVE REAL, PRIMARY KEY(RGroupID, Month));";
	char *table_rootsXphen = "CREATE TABLE sxwRootsXphen(RGroupID INT NOT NULL, Layer INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(RGroupID, Layer));";
	char *table_rootsSum =
			"CREATE TABLE sxwRootsSum(YEAR INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, VegProdType INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(Year, Iteration, Layer, VegProdType));";
	char *table_rootsRelative =
			"CREATE TABLE sxwRootsRelative(YEAR INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, RGroupID INT NOT NULL, January REAL, February REAL, March REAL, April REAL, May REAL, June REAL, July REAL, August REAL, September REAL, October REAL, November REAL, December REAL, PRIMARY KEY(Year, Iteration, Layer, RGroupID));";

	char *table_InputVars =
			"CREATE TABLE sxwInputVars(Year INT NOT NULL, Iteration INT NOT NULL, FracGrass REAL, FracShrub REAL, FracTree REAL, FracForb REAL, FracBareGround REAL, PRIMARY KEY(Year, Iteration));";
	char *table_InputProd =
			"CREATE TABLE sxwInputProd(Year INT NOT NULL, Iteration INT NOT NULL, VegProdType INT NOT NULL, Month INT NOT NULL, Litter REAL, Biomass REAL, PLive REAL, LAI_conv REAL, PRIMARY KEY(Year, Iteration, VegProdType, Month));";
	char *table_InputSoils =
			"CREATE TABLE sxwInputSoils(Year INT NOT NULL, Iteration INT NOT NULL, Layer INT NOT NULL, Tree_trco REAL, Shrub_trco REAL, Grass_trco REAL, Forb_trco REAL, PRIMARY KEY(Year, Iteration, Layer));";
	char *table_OutputVars =
			"CREATE TABLE sxwOutputVars(Year INT NOT NULL, Iteration INT NOT NULL, MAP_mm INT, MAT_C REAL, AET_cm REAL, T_cm REAL, ADT_cm REAL, AT_cm REAL, TotalRelsize REAL, TotalPR REAL, TotalTransp REAL, PRIMARY KEY(Year, Iteration));";
	char *table_OutputRgroup =
			"CREATE TABLE sxwOutputRgroup(YEAR INT NOT NULL, Iteration INT NOT NULL, RGroupID INT NOT NULL, Biomass REAL, Realsize REAL, PR REAL, resource_cur_preBvt Real,resource_cur REAL, PRIMARY KEY(Year, Iteration, RGroupID));";
	char *table_OutputProd =
			"CREATE TABLE sxwOutputProd(YEAR INT NOT NULL, Iteration INT NOT NULL, Month INT NOT NULL, BMass REAL, PctLive REAL, LAIlive REAL, LAItotal REAL, TotAGB REAL, PRIMARY KEY(Year, Iteration, Month));";
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

	rc = sqlite3_exec(db, table_rootsSum, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_rootsRelative, callback, 0, &zErrMsg);
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

	rc = sqlite3_exec(db, table_OutputTransp, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_OutputSWCBulk, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	insertRgroups();
	prepareStatements();
}
