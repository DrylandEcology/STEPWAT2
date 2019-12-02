/*
 * ST_sql.c
 *
 *  Created on: Feb 18, 2015
 *      Author: Ryan J. Murphy
 */


#include <stdio.h>
#include <string.h>
#include <sqlite3.h>
#include "ST_steppe.h"
#include "ST_globals.h"

static sqlite3 *db;
static char sql[1024];

static void beginTransaction(void);
static void endTransaction(void);
static void prepareStatements(void);
static void finalizeStatements(void);

static sqlite3_stmt * stmt_RGroupsYearInfo;
static sqlite3_stmt * stmt_SpeciesYearInfo;
static sqlite3_stmt * stmt_Indiv;
static sqlite3_stmt * stmt_IndivYearInfo;
static sqlite3_stmt * stmt_IndivKill;

static void createTables(void);

void ST_connect(char *stdbName) {
	int rc;
	char name[256] = { 0 };
	char *zErrMsg = 0;
	strcat(name, stdbName);
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

	createTables();
}

void ST_disconnect() {
	finalizeStatements();
	sqlite3_close(db);
}

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

static void prepareStatements() {
	sql[0] = 0;

	sprintf(sql, "INSERT INTO RGroupsYearInfo (Year, Iteration, RGroupID, Estabs, KillYr, YrsNegPR, mmExtraRes, ResRequired, ResAvail, ResExtra, PR, RelSize, EstSppCount, Extirpated, RegenOk) VALUES (@Year, @Iteration, @RGroupID, @Estabs, @KillYr, @YrsNegPR, @mmExtraRes, @ResRequired, @ResAvail, @ResExtra, @PR, @RelSize, @EstSppCount, @Extirpated, @RegenOk);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_RGroupsYearInfo, NULL);
        
	sprintf(sql, "INSERT INTO SpeciesYearInfo (Year, Iteration, SpeciesID, EstabCount, Estabs, RelSize, ExtraGrowth, ReceivedProb, AllowGrowth, sdSGerm) VALUES (@Year, @Iteration, @SpeciesID, @EstabCount, @Estabs, @RelSize, @ExtraGrowth, @ReceivedProb, @AllowGrowth, @sdSGerm);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_SpeciesYearInfo, NULL);

	sprintf(sql, "INSERT INTO Indiv (IndivID, Iteration, CreatedYear, SpeciesID, RGroupID) VALUES (@IndivID, @Iteration, @CreatedYear, @SpeciesID, @RGroupID);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_Indiv, NULL);

	sprintf(sql, "INSERT INTO IndivYearInfo(Year, IndivID, MortalityTypeID, age, mmExtraRes, SlowYrs, YrsNegPR, Killed, RelSize, GrpResProp, ResRequired, ResAvail, ResExtra, PR, GrowthRate, ProbVeggrow) VALUES (@Year, @IndivID, @MortalityTypeID, @age, @mmExtraRes, @SlowYrs, @YrsNegPR, @Killed, @RelSize, @GrpResProp, @ResRequired, @ResAvail, @ResExtra, @PR, @GrowthRate, @ProbVeggrow);");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_IndivYearInfo, NULL);

	sprintf(sql, "UPDATE Indiv SET KilledYear=?, KillTypeID=? WHERE IndivID=?;");
	sqlite3_prepare_v2(db, sql, 1024, &stmt_IndivKill, NULL);
}

static void finalizeStatements() {
	sqlite3_finalize(stmt_RGroupsYearInfo);
	sqlite3_finalize(stmt_SpeciesYearInfo);
	sqlite3_finalize(stmt_IndivYearInfo);
	sqlite3_finalize(stmt_IndivKill);
}

void insertIndivKill(int IndivID, int KillTypeID) {
	sqlite3_bind_int(stmt_IndivKill, 1, Globals->currYear);
	sqlite3_bind_int(stmt_IndivKill, 2, KillTypeID);
	sqlite3_bind_int(stmt_IndivKill, 3, IndivID);

	sqlite3_step(stmt_IndivKill);
	sqlite3_clear_bindings(stmt_IndivKill);
	sqlite3_reset(stmt_IndivKill);
}

static void insertIndivYearInfoRow(int Year, int IndivID, int MortalityTypeID, int age, int mmExtraRes, int SlowYrs, int YrsNegPR, int Killed, float RelSize, float GrpResProp, float ResRequired, float ResAvail, float ResExtra, float PR, float GrowthRate, float ProbVeggrow) {
	sqlite3_bind_int(stmt_IndivYearInfo, 1, Year);
	sqlite3_bind_int(stmt_IndivYearInfo, 2, IndivID);
	sqlite3_bind_int(stmt_IndivYearInfo, 3, MortalityTypeID);
	sqlite3_bind_int(stmt_IndivYearInfo, 4, age);
	sqlite3_bind_int(stmt_IndivYearInfo, 5, mmExtraRes);
	sqlite3_bind_int(stmt_IndivYearInfo, 6, SlowYrs);
	sqlite3_bind_int(stmt_IndivYearInfo, 7, YrsNegPR);
	sqlite3_bind_int(stmt_IndivYearInfo, 8, Killed);
	sqlite3_bind_double(stmt_IndivYearInfo, 9, RelSize);
	sqlite3_bind_double(stmt_IndivYearInfo, 10, GrpResProp);
	sqlite3_bind_double(stmt_IndivYearInfo, 11, ResRequired);
	sqlite3_bind_double(stmt_IndivYearInfo, 12, ResAvail);
	sqlite3_bind_double(stmt_IndivYearInfo, 13, ResExtra);
	sqlite3_bind_double(stmt_IndivYearInfo, 14, PR);
	sqlite3_bind_double(stmt_IndivYearInfo, 15, GrowthRate);
	sqlite3_bind_double(stmt_IndivYearInfo, 16, ProbVeggrow);
	sqlite3_step(stmt_IndivYearInfo);
	sqlite3_clear_bindings(stmt_IndivYearInfo);
	sqlite3_reset(stmt_IndivYearInfo);
}

void insertIndivYearInfo(IndivType *ind) {
	//Take off of age because this happens after rgroup_IncrAges remove if you put before.
	insertIndivYearInfoRow(Globals->currYear, ind->id, ind->killedby, (ind->age - 1), ind->mm_extra_res, ind->slow_yrs, ind->yrs_neg_pr, ind->killed, ind->relsize, ind->grp_res_prop, ind->res_required, ind->res_avail, ind->res_extra, ind->pr, ind->growthrate, ind->prob_veggrow);
}

static void insertIndivRow(int IndivID, int Iteration, int CreatedYear, int SpeciesID, int RGroupID) {
	sqlite3_bind_int(stmt_Indiv, 1, IndivID);
	sqlite3_bind_int(stmt_Indiv, 2, Iteration);
	sqlite3_bind_int(stmt_Indiv, 3, CreatedYear);
	sqlite3_bind_int(stmt_Indiv, 4, SpeciesID);
	sqlite3_bind_int(stmt_Indiv, 5, RGroupID);
	sqlite3_step(stmt_Indiv);
	sqlite3_clear_bindings(stmt_Indiv);
	sqlite3_reset(stmt_Indiv);
}

void insertIndiv(IndivType *ind) {
	insertIndivRow(ind->id, Globals->currIter, Globals->currYear, ind->myspecies+1, Species[ind->myspecies]->res_grp+1);
}

static void insertSpeciesYearInfoRow(int Year, int Iteration, int SpeciesID, int EstabCount, int Estabs, float RelSize, float ExtraGrowth, float ReceivedProb, int AllowGrowth, int sdSGerm) {
	sqlite3_bind_int(stmt_SpeciesYearInfo, 1, Year);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 2, Iteration);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 3, SpeciesID);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 4, EstabCount);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 5, Estabs);
	sqlite3_bind_double(stmt_SpeciesYearInfo, 6, RelSize);
	sqlite3_bind_double(stmt_SpeciesYearInfo, 7, ExtraGrowth);
	sqlite3_bind_double(stmt_SpeciesYearInfo, 8, ReceivedProb);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 9, AllowGrowth);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 10, sdSGerm);

	sqlite3_step(stmt_SpeciesYearInfo);
	sqlite3_clear_bindings(stmt_SpeciesYearInfo);
	sqlite3_reset(stmt_SpeciesYearInfo);
}

void insertSpecieYearInfo(SppIndex s) {
	SpeciesType *sp = Species[s];
	insertSpeciesYearInfoRow(Globals->currYear, Globals->currIter, s+1, sp->est_count, sp->estabs, getSpeciesRelsize(s), sp->extragrowth, sp->received_prob, sp->allow_growth, sp->sd_sgerm);
}

static void insertRGroupYearInfoRow(int Year, int Iteration, int RGroupID, int Estabs, int KillYr, int YrsNegPR, float mmExtraRes, float ResRequired, float ResAvail, float ResExtra, float PR, float RelSize, int EstSppCount, int Extirpated, int RegenOk) {
	sqlite3_bind_int(stmt_RGroupsYearInfo, 1, Year);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 2, Iteration);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 3, RGroupID);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 4, Estabs);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 5, KillYr);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 6, YrsNegPR);
	sqlite3_bind_double(stmt_RGroupsYearInfo, 7, mmExtraRes);
	sqlite3_bind_double(stmt_RGroupsYearInfo, 8, ResRequired);
	sqlite3_bind_double(stmt_RGroupsYearInfo, 9, ResAvail);
	sqlite3_bind_double(stmt_RGroupsYearInfo, 10, ResExtra);
	sqlite3_bind_double(stmt_RGroupsYearInfo, 11, PR);
	sqlite3_bind_double(stmt_RGroupsYearInfo, 12, RelSize);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 13, EstSppCount);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 14, Extirpated);
	sqlite3_bind_int(stmt_RGroupsYearInfo, 15, RegenOk);

	sqlite3_step(stmt_RGroupsYearInfo);
	sqlite3_clear_bindings(stmt_RGroupsYearInfo);
	sqlite3_reset(stmt_RGroupsYearInfo);
}

void insertRGroupYearInfo(GrpIndex g) {
	GroupType *rg = RGroup[g];
	insertRGroupYearInfoRow(Globals->currYear, Globals->currIter, g+1, rg->estabs, rg->killyr, rg->yrs_neg_pr, rg->mm_extra_res, rg->res_required, rg->res_avail, rg->res_extra, rg->pr, getRGroupRelsize(g), rg->est_count, rg->extirpated, rg->regen_ok);
}

static void insertRgroups(void) {
	int rc;
	GrpIndex g;
	GroupType *rg;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	ForEachGroup(g)
	{
		rg = RGroup[g];
		sql[0] = 0;
		sprintf(sql, "INSERT INTO RGroups (RGroupID, NAME, MaxStretch, MaxSppEstab, MaxSpp, MaxAge, StartYr, KillFreq, Extirp, GrpNum, VegProdType, MinResReq, MaxDensity, MaxPerSqm, MaxBmass, XGrow, SlowRate, PptSlope1, PptSlope2, PptSlope3, PptIntcpt1, PptIntcpt2, PptIntcpt3, Succulent, UseExtraRes, UseMe, UseMort, EstAnnually, DepthClassID) VALUES (%d, '%s', %d, %d, %d, %d, %d, %f, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %d);", g+1, rg->name, rg->max_stretch, rg->max_spp_estab, rg->max_spp, rg->max_age, rg->startyr, rg->killfreq, rg->extirp, rg->grp_num, rg->veg_prod_type, rg->min_res_req, rg->max_density, rg->max_per_sqm, rg->max_bmass, rg->xgrow, rg->slowrate, rg->ppt_slope[0], rg->ppt_slope[1], rg->ppt_slope[2], rg->ppt_intcpt[0], rg->ppt_intcpt[1], rg->ppt_intcpt[2], rg->succulent, rg->use_extra_res, rg->use_me, rg->use_mort, rg->est_annually, rg->depth);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
		sqlcheck(rc, zErrMsg);
	}
	endTransaction();
}

static void insertSpecies(void) {
	int rc;
	GrpIndex s;
	SpeciesType *sp;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	ForEachSpecies(s)
	{
		sp = Species[s];
		sql[0] = 0;
		sprintf(sql,
				"INSERT INTO Species (SpeciesID, RGroupID, NAME, MaxAge, ViableYrs, MaxSeedEstab, MaxVegUnits, MaxSlow, SPnum, MaxRate, IntrinRate, RelSeedlingsSize, SeedlingBiomass, MatureBiomass, SeedlingEstabProbOld, SeedlingEstabProb, AnnMortProb, CohortSurv, ExpDecay, ProbVeggrow1, ProbVeggrow2, ProbVeggrow3, ProbVeggrow4, sdParam1, sdPPTdry, sdPPTwet, sdPmin, sdPmax, sdH, sdVT, TempClassID, DisturbClassID, isClonal, UseTempResponse, UseMe, UseDispersal) VALUES (%d, %d, '%s', %d, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %d);",
				s + 1, sp->res_grp + 1, sp->name, sp->max_age, sp->viable_yrs,
				sp->max_seed_estab, sp->max_vegunits, sp->max_slow, sp->sp_num,
				sp->max_rate, sp->intrin_rate, sp->relseedlingsize,
				sp->seedling_biomass, sp->mature_biomass,
				sp->seedling_estab_prob_old, sp->seedling_estab_prob,
				sp->ann_mort_prob, sp->cohort_surv, sp->exp_decay,
				sp->prob_veggrow[0], sp->prob_veggrow[1], sp->prob_veggrow[2],
				sp->prob_veggrow[3], sp->sd_Param1, sp->sd_PPTdry,
				sp->sd_PPTwet, sp->sd_Pmin, sp->sd_Pmax, sp->sd_H, sp->sd_VT,
				sp->tempclass, sp->disturbclass, sp->isclonal,
				sp->use_temp_response, sp->use_me, sp->use_dispersal);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
		sqlcheck(rc, zErrMsg);
	}
	endTransaction();
}

static void insertInfo(void) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql,
			"INSERT INTO info (Years, Iterations, Seed, RGroups, Species,  PlotSize) VALUES (%d, %d, %d, %d, %d, %f);",
			SuperGlobals.runModelYears, SuperGlobals.runModelIterations, (int) SuperGlobals.randseed,
			Globals->grpCount, Globals->sppCount, Globals->plotsize);
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void insertTempClass(void) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql, "INSERT INTO TempClass (TempClassID, Name) VALUES (0, 'No Season');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO TempClass (TempClassID, Name) VALUES (1, 'Cool Season');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO TempClass (TempClassID, Name) VALUES (2, 'Warm Season');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void insertDisturbClass(void) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql, "INSERT INTO DisturbClass (DisturbClassID, Name) VALUES (0, 'Very Sensitive');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DisturbClass (DisturbClassID, Name) VALUES (1, 'Sensitive');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DisturbClass (DisturbClassID, Name) VALUES (2, 'Insensitive');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DisturbClass (DisturbClassID, Name) VALUES (3, 'Very Insensitive');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void insertDepthClass(void) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql, "INSERT INTO DepthClass (DepthClassID, Name) VALUES (0, 'Depth Non Comp');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DepthClass (DepthClassID, Name) VALUES (1, 'Depth Shallow');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DepthClass (DepthClassID, Name) VALUES (2, 'Depth Medium');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DepthClass (DepthClassID, Name) VALUES (3, 'Depth Deep');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO DepthClass (DepthClassID, Name) VALUES (4, 'Depth Last');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void insertKillTypes(void) {
	int rc;
	char *zErrMsg = 0;
	sql[0] = 0;

	beginTransaction();
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (0, 'Species_Kill New Iteration / Plot Init');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (1, 'Species_Kill Pat');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (2, 'Species_Kill Mound');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (3, 'Species_Kill Burrow');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (4, 'Species_Kill Annuals');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (5, 'mort_EndOfYear rgroup_Extirpate Species_Kill');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (6, 'mort_EndOfYear RGroup_Kill Species_Kill');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (7, '_succulents');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (8, '_slow_growth');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (9, '_age_independent');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (10, '_no_resources');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
	sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (11, '_stretched_clonal');");
	rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);
        sprintf(sql, "INSERT INTO KillTypes (KillTypeID, Name) VALUES (12, '_kill_old_plants');");
        rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
        sqlcheck(rc, zErrMsg);
	endTransaction();
}

static void createTables(void) {
	int rc;
	char *zErrMsg = 0;

	char *table_PrjInfo = "CREATE TABLE info(Years INT, Iterations INT, Seed INT, RGroups INT, Species INT,  PlotSize REAL);";
	char *table_tempClass = "CREATE TABLE TempClass(TempClassID INT PRIMARY KEY NOT NULL, Name TEXT);";
	char *table_disturbClass = "CREATE TABLE DisturbClass(DisturbClassID INT PRIMARY KEY NOT NULL, Name TEXT);";
	char *table_depthClass = "CREATE TABLE DepthClass(DepthClassID INT PRIMARY KEY NOT NULL, Name TEXT);";
	char *table_killTypes = "CREATE TABLE KillTypes(KillTypeID INT PRIMARY KEY, Name TEXT);";

	char *table_rgroups = "CREATE TABLE RGroups(RGroupID INT PRIMARY KEY NOT NULL, NAME TEXT NOT NULL, MaxStretch INT, MaxSppEstab INT, MaxSpp INT, MaxAge INT, StartYr INT, KillFreq REAL, Extirp INT, GrpNum INT, VegProdType INT, MinResReq REAL, MaxDensity REAL, MaxPerSqm REAL, MaxBmass REAL, XGrow REAL, SlowRate REAL, PptSlope1 REAL, PptSlope2 REAL, PptSlope3 REAL, PptIntcpt1 REAL, PptIntcpt2 REAL, PptIntcpt3 REAL, Succulent INT, UseExtraRes INT, UseMe INT, UseMort INT, EstAnnually INT, DepthClassID INT );";
	char *table_rgroupsYearInfo = "CREATE TABLE RGroupsYearInfo(Year INT NOT NULL, Iteration INT NOT NULL, RGroupID INT NOT NULL, Estabs INT, KillYr INT, YrsNegPR INT, mmExtraRes INT, ResRequired REAL, ResAvail REAL, ResExtra REAL, PR REAL, RelSize REAL, EstSppCount REAL, Extirpated INT, RegenOk INT, PRIMARY KEY(Year, Iteration, RGroupID));";

	char *table_species = "CREATE TABLE Species(SpeciesID INT PRIMARY KEY NOT NULL, RGroupID INT NOT NULL, NAME TEXT NOT NULL, MaxAge INT, ViableYrs INT, MaxSeedEstab INT, MaxVegUnits INT, MaxSlow INT, SPnum INT, MaxRate REAL, IntrinRate REAL, RelSeedlingsSize REAL, SeedlingBiomass REAL, MatureBiomass REAL, SeedlingEstabProbOld REAL, SeedlingEstabProb REAL, AnnMortProb REAL, CohortSurv REAL, ExpDecay REAL, ProbVeggrow1 REAL, ProbVeggrow2 REAL, ProbVeggrow3 REAL, ProbVeggrow4 REAL, sdParam1 REAL, sdPPTdry REAL, sdPPTwet REAL, sdPmin REAL, sdPmax REAL, sdH REAL, sdVT REAL, TempClassID INT, DisturbClassID INT, isClonal INT, UseTempResponse INT, UseMe INT, UseDispersal INT);";
	char *table_speciesYearInfo = "CREATE TABLE SpeciesYearInfo(Year INT NOT NULL, Iteration INT NOT NULL, SpeciesID INT NOT NULL, EstabCount INT, Estabs INT, RelSize REAL, ExtraGrowth REAL, ReceivedProb REAL, AllowGrowth INT, sdSGerm INT, PRIMARY KEY(Year, Iteration, SpeciesID));";

	char *table_indiv = "CREATE TABLE Indiv(IndivID INT NOT NULL, Iteration INT NOT NULL, CreatedYear INT NOT NULL, SpeciesID INT NOT NULL, RGroupID INT NOT NULL, KilledYear INT, KillTypeID INT, PRIMARY KEY(IndivID, Iteration, CreatedYear, SpeciesID));";
	char *table_indivYearInfo = "CREATE TABLE IndivYearInfo(Year INT NOT NULL, IndivID INT NOT NULL, MortalityTypeID INT, age INT, mmExtraRes INT, SlowYrs INT, YrsNegPR INT, Killed INT, RelSize REAL, GrpResProp REAL, ResRequired REAL, ResAvail REAL, ResExtra REAL, PR REAL, GrowthRate REAL, ProbVeggrow REAL, PRIMARY KEY(Year, IndivID));";
	//char *table_indivKill = "CREATE TABLE IndivKill(Year INT NOT NULL, IndivID INT NOT NULL, KillTypeID INT NOT NULL, PRIMARY KEY(Year, IndivID));";


	rc = sqlite3_exec(db, table_PrjInfo, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_tempClass, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_disturbClass, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_depthClass, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_rgroups, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_rgroupsYearInfo, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_species, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_speciesYearInfo, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_indiv, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_indivYearInfo, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	//rc = sqlite3_exec(db, table_indivKill, callback, 0, &zErrMsg);
	//sqlcheck(rc, zErrMsg);

	rc = sqlite3_exec(db, table_killTypes, callback, 0, &zErrMsg);
	sqlcheck(rc, zErrMsg);

	insertInfo();
	insertTempClass();
	insertDisturbClass();
	insertDepthClass();
	insertRgroups();
	insertSpecies();
	insertKillTypes();

	prepareStatements();
}
