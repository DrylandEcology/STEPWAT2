/** 
 * \file ST_sql.c
 * \brief A collection of functions that generate an SQL database output file.
 * 
 * Using SQLite in C can be intimidating if you aren't familiar with the API. 
 * You can learn about the API [here](https://www.sqlite.org/cintro.html).
 * 
 * \ingroup SQL
 * \date 18 February 2015
 * \author Ryan J. Murphy
 */

#include <string.h>
#include <sqlite3.h>
#include "ST_steppe.h"
#include "ST_globals.h"

/**
 * \brief A reference to the database object.
 * 
 * The database is opened or created in ST_connect(), and a reference is stored here.
 * 
 * \ingroup SQL
 */
static sqlite3 *db;

/**
 * \brief A buffer to store SQL statements.
 * \ingroup SQL
 */
static char sql[1024];

static void beginTransaction(void);
static void endTransaction(void);
static void prepareStatements(void);
static void finalizeStatements(void);

/**
 * \brief The [resource group](\ref GroupType) year row SQL statement.
 * \ingroup SQL
 */
static sqlite3_stmt * stmt_RGroupsYearInfo;
/**
 * \brief The [Species](\ref SpeciesType) year row SQL statement.
 * \ingroup SQL
 */
static sqlite3_stmt * stmt_SpeciesYearInfo;
/**
 * \brief The [individual](\ref IndivType) SQL statement.
 * \ingroup SQL
 */
static sqlite3_stmt * stmt_Indiv;
/**
 * \brief The [individual](\ref IndivType) year row SQL statement.
 * \ingroup SQL
 */
static sqlite3_stmt * stmt_IndivYearInfo;
/**
 * \brief The [individual](\ref IndivType) mortality SQL statement.
 * \ingroup SQL
 */
static sqlite3_stmt * stmt_IndivKill;

static void createTables(void);

/**
 * \brief Creates a connection to the database named \p stdbName
 * 
 * \param stdbName is the name of the output file.
 * 
 * \sideeffect 
 *      Opens the database for writing. If no database exists named stdbName a
 *      file will be created.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Disconnects from the database opened in ST_connect()
 * 
 * \sideeffect Finalizes the database and closes the reference.
 * 
 * \ingroup SQL
 */
void ST_disconnect() {
	finalizeStatements();
	sqlite3_close(db);
}

/**
 * \brief Prepares the database to recieve information.
 * 
 * This is equivalent to the `BEGIN TRANSACTION` command in SQL.
 * This function is called inside of the insertion functions so there
 * is no need to call in explicitly unless you are writting a new
 * insertion function.
 * 
 * \ingroup SQL
 */
static void beginTransaction() {
	char * sErrMsg = 0;
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);
}

/**
 * \brief Ends a transaction with the database
 * 
 * This is equivalent to the `END TRANSACTION` command in SQL.
 * This function is called inside of the insertion functions so there
 * is no need to call in explicitly unless you are writting a new
 * insertion function.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Determines if an SQL interaction was successfull.
 * 
 * \param rc is the code returned when executing an SQL instruction.
 * \param zErrMsg is the string to print if \p rc indicates that there
 *                was an error
 * 
 * \sideeffect If an error occured \p zErrMsg is printed to stdout.
 * 
 * \ingroup SQL
 */
static void sqlcheck(int rc, char *zErrMsg) {
	if (rc != SQLITE_OK) {
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	} else {
		//fprintf(stdout, "Table created successfully\n");
	}
}

/**
 * \brief Prepares the database to receive \ref RGROUP, \ref SPECIES, \ref INDIVIDUAL and \ref MORTALITY information.
 * 
 * This must be called before the actual data is inserted.
 * 
 * \ingroup SQL
 * 
 * \sa finalizeStatements() which is the sister function to prepareStatements()
 */
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

/**
 * \brief finalizes the insertion of \ref RGROUP, \ref SPECIES, \ref INDIVIDUAL, and \ref MORTALITY data.
 * 
 * This must be called following the insertion of data.
 * 
 * \ingroup SQL
 * 
 * \sa prepareStatements() which is the sister function to finalizeStatements()
 */
static void finalizeStatements() {
	sqlite3_finalize(stmt_RGroupsYearInfo);
	sqlite3_finalize(stmt_SpeciesYearInfo);
	sqlite3_finalize(stmt_IndivYearInfo);
	sqlite3_finalize(stmt_IndivKill);
}

/** 
 * \brief Inserts \ref MORTALITY information pertaining to a single 
 *        \ref INDIVIDUAL into the database.
 * 
 * \param IndivID the unique ID of the [individual](\ref IndivType).
 * \param KillTypeID the \ref MortalityType code that killed the [individual](\ref IndivType).
 * 
 * \sideeffect The database is populated with the given data.
 * 
 * \ingroup SQL
 */
void insertIndivKill(int IndivID, int KillTypeID) {
	sqlite3_bind_int(stmt_IndivKill, 1, Globals->currYear);
	sqlite3_bind_int(stmt_IndivKill, 2, KillTypeID);
	sqlite3_bind_int(stmt_IndivKill, 3, IndivID);

	sqlite3_step(stmt_IndivKill);
	sqlite3_clear_bindings(stmt_IndivKill);
	sqlite3_reset(stmt_IndivKill);
}

/**
 * \brief Insert a row into the database for a given [individual](\ref IndivType).
 * 
 * \param Year is the year that this data was collected.
 * \param IndivID is the unique ID of the [individual](\ref IndivType).
 * \param MortalityTypeID is the \ref MortalityType code that killed the [individual](\ref IndivType).
 * \param age is the age of the [individual](\ref IndivType).
 * \param mmExtraRes is the amount of extra resources, in millimeters, that the [individual](\ref IndivType) recieved.
 * \param SlowYrs are the number of years that the [individual](\ref IndivType) has been experiencing slow growth.
 * \param YrsNegPR are the number of years that the [individual](\ref IndivType)'s PR has been negative.
 * \param Killed is treated as a bool. TRUE if the [individual](\ref IndivType) was killed this year.
 * \param RelSize is the relative size of the [individual](\ref IndivType).
 * \param GrpResProp is the proportion of the [group](\ref RGROUP)'s total resources used by this [individual](\ref IndivType).
 * \param ResRequired is the amount of resources this [individual](\ref IndivType) requires to grow.
 * \param ResAvail is the amount to resources availible to this [individual](\ref IndivType).
 * \param ResExtra is the amount of resources that this [individual](\ref IndivType) recieved which it did not require.
 * \param PR is the PR of this [individual](\ref IndivType).
 * \param GrowthRate is the intrinsic growth rate of this [individual](\ref IndivType).
 * \param ProbVeggrow is the probability that this [individual](\ref IndivType) produces vegetative growth.
 * 
 * \sideeffect The database will be populated with this information.
 * 
 * \ingroup SQL
 * 
 * It would be a pain to pull out this much information from an [individual](\ref IndivType) every time, so instead use
 * insertIndivYearInfo() which takes a pointer to an individual and and pulls the necessary information then calls this 
 * function.
 */
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

/**
 * \brief Inserts an [individual](\ref IndivType)'s information into the database.
 * 
 * \param ind is a pointer to an [individual](\ref IndivType).
 * 
 * \sideeffect The [individual](\ref IndivType)'s information will be inserted into the database.
 * 
 * \sa insertIndivYearInfoRow() which is called in this function to bind the data.
 * 
 * \ingroup SQL
 */
void insertIndivYearInfo(IndivType *ind) {
	//Take off of age because this happens after rgroup_IncrAges remove if you put before.
	insertIndivYearInfoRow(Globals->currYear, ind->id, ind->killedby, (ind->age - 1), ind->mm_extra_res, ind->slow_yrs, ind->yrs_neg_pr, ind->killed, ind->relsize, ind->grp_res_prop, ind->res_required, ind->res_avail, ind->res_extra, ind->pr, ind->growthrate, ind->prob_veggrow);
}

/**
 * \brief Inserts big-picture information about an [individual](\ref IndivType) as a new row.
 * 
 * \param IndivID is the unique ID of the [individual](\ref IndivType).
 * \param CreatedYear is the year in which the [individual](\ref IndivType) established.
 * \param SpeciesID is the [index](\ref SppIndex) in \ref Species of the species to which this [individual](\ref IndivType)
 *                  belongs.
 * \param RGroupID is the [index](\ref GrpIndex) in \ref RGroup of the resource group to which this 
 *                 [individual](\ref IndivType) belongs.
 * \param Iteration is the simulation's current iteration.
 * 
 * \sideeffect The given information is inserted into the database.
 * 
 * Calling a function with this many parameters would be a pain. Instead, call insertIndiv() which takes
 * a pointer to an [individual](\ref IndivType), extracts the parameters, then calls this function.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Inserts big-picture information about an [individual](\ref IndivType) as a row in the database.
 * 
 * \param ind is a pointer to an \ref IndivType.
 * 
 * \sideeffect Information about the [individual](\ref IndivType) will be inserted into the database.
 * 
 * \sa insertIndivRow() which is called by this function to bind the SQL data.
 * 
 * \ingroup SQL
 */
void insertIndiv(IndivType *ind) {
	insertIndivRow(ind->id, Globals->currIter, Globals->currYear, ind->myspecies+1, Species[ind->myspecies]->res_grp+1);
}

/**
 * \brief Insert information on a given [species](\ref SpeciesType).
 * 
 * \param Year is the year in which the data was collected.
 * \param Iteration is the current iteration.
 * \param SpeciesID is the [index](\ref SppIndex) is \ref Species of the given [species](\ref SpeciesType).
 * \param EstabCount is the number of established [individuals](\ref IndivType) in the [species](\ref SpeciesType).
 * \param Estabs is the number of [individuals](\ref IndivType) that established in the current year.
 * \param RelSize is the relative size of the [species](\ref SpeciesType).
 * \param ExtraGrowth is the amount of superfluous growth this [species](\ref SpeciesType) experienced in the current year.
 * \param ReceivedProb is the probability that this [species](\ref SpeciesType) received seeds from an adjacent plot.
 * \param AllowGrowth is treaded as a \ref Bool. TRUE if this [species](\ref SpeciesType) was allowed to grow this year.
 * \param sdSGerm I'm not sure what this is. It appears to be something regarding seed dispersal which is currently
 *                an [open issue](https://github.com/DrylandEcology/STEPWAT2/issues/309) on GitHub.
 * 
 * \sideeffect The provided parameters will be inserted into the database.
 * 
 * Providing this many parameters is a waste of type. Instead call insertSpecieYearInfo() which takes an [index](\ref SppIndex)
 * in \ref Species, exracts the information, then calls this function.
 * 
 * \ingroup SQL
 */
static void insertSpeciesYearInfoRow(int Year, int Iteration, int SpeciesID, int EstabCount, int Estabs, float RelSize, float ExtraGrowth, float ReceivedProb, int sdSGerm) {
	sqlite3_bind_int(stmt_SpeciesYearInfo, 1, Year);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 2, Iteration);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 3, SpeciesID);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 4, EstabCount);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 5, Estabs);
	sqlite3_bind_double(stmt_SpeciesYearInfo, 6, RelSize);
	sqlite3_bind_double(stmt_SpeciesYearInfo, 7, ExtraGrowth);
	sqlite3_bind_double(stmt_SpeciesYearInfo, 8, ReceivedProb);
	sqlite3_bind_int(stmt_SpeciesYearInfo, 9, sdSGerm);

	sqlite3_step(stmt_SpeciesYearInfo);
	sqlite3_clear_bindings(stmt_SpeciesYearInfo);
	sqlite3_reset(stmt_SpeciesYearInfo);
}

/**
 * \brief Insert a [species](\ref SpeciesType)'s information into the database.
 * 
 * \param s is the [index](\ref SppIndex) in \ref Species of the [species](\ref SpeciesType) to insert.
 * 
 * \sideeffect Information about the species will be inserted into the database.
 * 
 * \sa insertSpeciesYearInfoRow() which is responsible for binding the data.
 * 
 * \ingroup SQL
 */
void insertSpecieYearInfo(SppIndex s) {
	SpeciesType *sp = Species[s];
	insertSpeciesYearInfoRow(Globals->currYear, Globals->currIter, s+1, sp->est_count, sp->estabs, getSpeciesRelsize(s), sp->extragrowth, sp->received_prob, sp->seedsPresent);
}

/**
 * \brief Insert a [resource group](\ref GroupType)'s information into the database.
 * 
 * \param Year is the year in the simulation which generated the information.
 * \param Iteration is the iteration in the simulation.
 * \param RGroupID is the [index](\ref GrpIndex) in \ref RGroup of the given [group](\ref GroupType).
 * \param Estabs is the number of established [individuals](\ref IndivType) in the [group](\ref GroupType).
 * \param KillYr is the last year in the simulation that had a disturbance.
 * \param YrsNegPR is the number of years that the [resource group](\ref GroupType) has experienced negative PR.
 * \param mmExtraRes is the millimeters of extra resources the [resource group](\ref GroupType) has received this year.
 * \param ResRequired is the amount of resources required to sustain all individuals in the [resource group](\ref GroupType).
 * \param ResAvail is the amount of resources available to the [resource group](\ref GroupType) in the current year.
 * \param ResExtra is the amount of extra resources available to the [resource group](\ref GroupType) in the current year.
 * \param PR is the ratio of resources required to resources available.
 * \param RelSize is the relative size of the [resource group](\ref GroupType).
 * \param EstSppCount is the number of [species](\ref SpeciesType) in the [resource group](\ref GroupType) with at least
 *        one [individual](\ref IndivType) established.
 * \param Extirpated is TRUE if this [resource group](\ref GroupType) has been killed and prevented from regenerating.
 * \param RegenOk is only for annuals. TRUE if the [resource group](\ref GroupType) is allowed to regenerate.
 * 
 * Filling in all of these parameters would be a pain. Instead call insertRGroupYearInfo() which only takes a group pointer
 * and extracts all of these parameters.
 * 
 * \sideeffect A row of resource group data will be inserted into the database.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Insert a row of [resource group](\ref GroupType) information into the database.
 * 
 * \param g is the [index](\ref GrpIndex) in \ref RGroup of the requested resource group.
 * 
 * \ingroup SQL
 * 
 * \sideeffect A row of resouce group data will be inserted into the database.
 * 
 * \sa insertRGroupYearInfoRow() which is called by this function to bind the data.
 */
void insertRGroupYearInfo(GrpIndex g) {
	GroupType *rg = RGroup[g];
	insertRGroupYearInfoRow(Globals->currYear, Globals->currIter, g+1, rg->estabs, rg->killyr, rg->yrs_neg_pr, rg->mm_extra_res, rg->res_required, rg->res_avail, rg->res_extra, rg->pr, getRGroupRelsize(g), rg->est_count, rg->extirpated, rg->regen_ok);
}

/**
 * \brief Inserts all [resource groups](\ref RGroup) into the database.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Inserts information for all [species](\ref Species) into the database.
 * 
 * \ingroup SQL
 */
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
				"INSERT INTO Species (SpeciesID, RGroupID, NAME, MaxAge, ViableYrs, MaxSeedEstab, MaxVegUnits, MaxSlow, SPnum, MaxRate, IntrinRate, RelSeedlingsSize, SeedlingBiomass, MatureBiomass, SeedlingEstabProbOld, SeedlingEstabProb, AnnMortProb, CohortSurv, ExpDecay, ProbVeggrow1, ProbVeggrow2, ProbVeggrow3, ProbVeggrow4, minReproductiveSize, Height, TempClassID, DisturbClassID, isClonal, UseTempResponse, UseMe, UseDispersal) VALUES (%d, %d, '%s', %d, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %d);",
				s + 1, sp->res_grp + 1, sp->name, sp->max_age, sp->viable_yrs,
				sp->max_seed_estab, sp->max_vegunits, sp->max_slow, sp->sp_num,
				sp->max_rate, sp->intrin_rate, sp->relseedlingsize,
				sp->seedling_biomass, sp->mature_biomass,
				sp->seedling_estab_prob_old, sp->seedling_estab_prob,
				sp->ann_mort_prob, sp->cohort_surv, sp->exp_decay,
				sp->prob_veggrow[0], sp->prob_veggrow[1], sp->prob_veggrow[2],
				sp->prob_veggrow[3], sp->minReproductiveSize,
				getSpeciesHeight(sp), sp->tempclass, sp->disturbclass,
				sp->isclonal, sp->use_temp_response, sp->use_me, 
				sp->use_dispersal);
		rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
		sqlcheck(rc, zErrMsg);
	}
	endTransaction();
}

/**
 * \brief Inserts general information about the simulation.
 * 
 * As of 9/9/19 this information is the number of years, number of iterations, the 
 * random seed, the number of groups, the number of species, and the plot size.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Insert the index values of the temperature classes.
 * 
 * These values are hard coded, and I'm not sure they are correct. This should NOT
 * be used as a reference. Instead see the \ref TempClass enumerator.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Insert the integer values assigned to each \ref DisturbClass.
 * 
 * This function is hard coded and does not actually reference the \ref DisturbClass enum.
 * Therefore, I would not use it as a reference.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Insert the mapping of depth classes to integers.
 * 
 * This value is hard coded which means it could be wrong. Instead see the \ref DepthClass
 * enumerator.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Insert the mapping of \ref MortalityType to integers.
 * 
 * This mapping is hard coded so do not assume it is correct. Instead reference
 * [the enum](\ref MortalityType) directly.
 * 
 * \ingroup SQL
 */
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

/**
 * \brief Creates the various tables that comprise the database. Consider this the main function.
 * 
 * These include the [species](\ref SpeciesType), [resource group](\ref GroupType) and [individual](\ref IndivType) tables plus
 * the enumerator mappings (I wouldn't trust the enumerator mappings because they are hard coded).
 * 
 * This function takes care of inserting everything EXCEPT the [rgroup](\ref GroupType), [species](\ref SpeciesType)
 * and [individual](\ref IndivType) row information.
 * 
 * This function also assumes ST_connect() has already been called to initialize [the database](\ref db)
 * 
 * \sideeffect the tables are created and (with the exception of row information) populated.
 * 
 * \ingroup SQL
 */
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

	char *table_species = "CREATE TABLE Species(SpeciesID INT PRIMARY KEY NOT NULL, RGroupID INT NOT NULL, NAME TEXT NOT NULL, MaxAge INT, ViableYrs INT, MaxSeedEstab INT, MaxVegUnits INT, MaxSlow INT, SPnum INT, MaxRate REAL, IntrinRate REAL, RelSeedlingsSize REAL, SeedlingBiomass REAL, MatureBiomass REAL, SeedlingEstabProbOld REAL, SeedlingEstabProb REAL, AnnMortProb REAL, CohortSurv REAL, ExpDecay REAL, ProbVeggrow1 REAL, ProbVeggrow2 REAL, ProbVeggrow3 REAL, ProbVeggrow4 REAL, minReproductiveSize REAL, sdPPTdry REAL, sdPPTwet REAL, sdPmin REAL, sdPmax REAL, sdH REAL, sdVT REAL, TempClassID INT, DisturbClassID INT, isClonal INT, UseTempResponse INT, UseMe INT, UseDispersal INT);";
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
