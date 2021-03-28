/**
 * \file ST_grid.h
 * \brief All exported functions, variables, ect. from the [grid](\ref GRID) 
 *        module
 * 
 * The gridded mode exports many things including enumerators, functions,
 * variables and structs. 
 * 
 * \author Chandler Haukap
 * \date August 2019. 
 * \ingroup GRID
 */

#ifndef GRID_H
#define GRID_H

/******** These modules are necessary to compile ST_grid.c ********/
#include "ST_stats.h"
#include "sw_src/pcg/pcg_basic.h"
#include "sxw_vars.h"
#include "sw_src/SW_SoilWater.h"
#include "sw_src/SW_VegProd.h"
#include "sw_src/SW_Model.h"
#include "sw_src/SW_Weather.h"
#include "ST_seedDispersal.h"

/*********************** Grid Structures ****************************/

/** 
 * \brief Soil layer information for a single [cell](\ref CellType)
 * 
 * This struct requires dynamic memory allocation, so be carefull when
 * instanciating it.
 * 
 * \ingroup GRID
 */
struct Soil_st
{
	/** \brief Number of soil layers. */
	int num_layers;
	/** Name of the roots file associated with this cell. */
	char rootsFile[20];
	/** \brief Depth in cm of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *depth;
    /** \brief matricd of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *matricd;
	/** \brief Percent gravel for each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *gravel;
	/** \brief evco of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *evco;
	/** \brief trco for grass of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *trco_grass;
	/** \brief trco for shrubs of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *trco_shrub;
	/** \brief trco for trees of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *trco_tree;
	/** \brief trco for forbs of each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *trco_forb;
	/** \brief Percent sand for each soil layer. Size of array defined by 
	 *         num_layers.  */
    RealF *psand;
	/** \brief Percent clay for each soil layer. size of array defined by 
	 *         num_layers.  */
    RealF *pclay;
	/** \brief impermiability of each soil layer. size of array defined by 
	 *         num_layers.  */
    RealF *imperm;
	/** \brief temperature of each soil layer. size of array defined by 
	 *         num_layers.  */
    RealF *soiltemp;
}typedef SoilType;

/** 
 * \brief [Spinup](\ref SPINUP) information for a given [cell](\ref CellType).
 * 
 * \ingroup GRID
 */
struct grid_init_species_st
{
	/** \brief TRUE if at least one species has requested spinup in
	 *         this cell. */
	int useSpinup;
	/** 
	 * \brief Array of boolean values that correspond to this
	 *         [cell](\ref CellType)'s [species](\ref Species) array.
	 * 
	 *  TRUE if Species[sp] should use spinup.
	 */
	int *shouldSpinup;
}typedef Grid_Init_Species_St;

/** 
 * \brief All information specific to a cell.
 * 
 * This struct is basically the cell. It holds all information specific to one
 * cell, and a 2d array of these structs makes up a gridded mode simulation.
 * 
 * \sa gridCells
 * \ingroup GRID
 */
struct grid_cell_st
{
	/** \brief RGroup corresponding to this cell */
	GroupType **myGroup;
	/** \brief Species corresponding to this cell */
	SpeciesType **mySpecies;
	/** \brief Succulents corresponding to this cell */
	SucculentType mySucculent;
	/** \brief This cell's environment. We expect each cell to
	 *         have slightly different weather each year */
	EnvType myEnvironment;
	/** \brief Cell's plot data */
	PlotType myPlot;
	/** \brief Global variables corresponding to this cell */ 
	ModelType myGlobals;
	/** \brief If TRUE this cell should use seed dispersal */
	Bool useSeedDispersal;
	/** \brief TRUE if this cell is in spinup mode */
	Bool DuringSpinup;
	/** \brief Species' [spinup](\ref SPINUP) information. */
	Grid_Init_Species_St mySpeciesInit;
	/** \brief This cell's cheatgrass-wildfire parameters. */
	CheatgrassPrecip* myCheatgrassPrecip;

	/** \brief TRUE if an individual was killed this year. */
	Bool* someKillage;

	/** \brief this cell's version of the \ref UseCheatgrassWildfire variable
	 *         from the \ref MORTALITY module. */
	Bool UseCheatgrassWildfire;
	
	/* ---------------- accumulators -------------------- */

	/** \brief The disturbance event statistics accumulator. */
	StatType *_Dist;
	/** \brief The precipitation statistics accumulator. */
	StatType *_Ppt;
	/** \brief The temperature statistics accumulator. */
	StatType *_Temp;	
	/** \brief The [rgroup](\ref GroupType) biomass statistics accumulator. */
  	StatType *_Grp;
	/** \brief The [rgroup](\ref GroupType) size statistics accumulator. */
	StatType *_Gsize;
	/** \brief The [rgroup](\ref GroupType) PR statistics accumulator. */ 
	StatType *_Gpr;
	/** \brief The [rgroup](\ref GroupType) mortality statistics
	 *         accumulator. */
	StatType *_Gmort; 
	/** \brief The [rgroup](\ref GroupType) establishment statistics 
	 *         accumulator. */ 
	StatType*_Gestab;
	/** \brief The [species](\ref SpeciesType) biomass statistics 
	 *         accumulator. */ 
  	StatType *_Spp;
	/** \brief The [individual](\ref IndivType) statistics accumulator. */ 
	StatType *_Indv;
	/** \brief The [species](\ref SpeciesType) mortality statistics 
	 *         accumulator. */ 
	StatType *_Smort;
	/** \brief The [species](\ref SpeciesType) establishment statistics 
	 *         accumulator. */
	StatType *_Sestab;
	/** \brief The [species](\ref SpeciesType) 
	 *         [seed dispersal](\ref SEED_DISPERSAL) statistics accumulator. */ 
	StatType *_Sreceived;
	/** \brief The [rgroup](\ref GroupType) wildfire statistics accumulator. */
	FireStatsType *_Gwf;

	/**
	 * \brief TRUE if this cell has initialized its statistics accumulators.
	 */
	Bool stats_init;
	/* -------------- end accumulators ------------------ */

	/* -------------------- SXW ------------------------- */
	/**
	 * \brief The transpiration window struct specific to this cell.
	 * 
	 * Intended to be swapped into the [SXW private](\ref SXW_PRIVATE) variable
	 * \ref transp_window.
	 * 
	 * \sa copy_sxw_variables
	 */
	transp_t* myTranspWindow;
	/**
	 * \brief The SXW struct specific to this cell. 
	 * 
	 * Intended to be swapped into the [SXW private](\ref SXW_PRIVATE) variable
	 * \ref SXW.
	 * 
	 * \sa copy_sxw_variables
	 */
	SXW_t* mySXW;
	/**
	 * \brief The SXW resource struct specific to this cell.
	 * 
	 * Intended to be swapped into the [SXW private](\ref SXW_PRIVATE) varable
	 * \ref SXWResources
	 * 
	 * \sa copy_sxw_variables
	 */
	SXW_resourceType* mySXWResources;
	/* ------------------ End SXW ----------------------- */

    /* ----------------- SOILWAT2 output ---------------- */
    /**
     * \brief If TRUE this cell should print SOIWLAT2 output.
     */
    Bool generateSWOutput;
    /* --------------- End SOILWAT2 output -------------- */

	/* ------------------- Soils ------------------------ */
	/** \brief  Soil layer information specific to this cell. */
	SoilType mySoils;
	/* ------------------ End Soils --------------------- */
} typedef CellType;

/**************************** Enumerators *********************************/
/** 
 * \brief Indices for the \ref grid_directories array.
 * \ingroup GRID
 */
typedef enum
{
	/** \brief The index in \ref grid_directories of the STEPWAT directory. */
    GRID_DIRECTORY_STEPWAT_INPUTS,

    /** \brief The number of directories in \ref grid_directories. 
	 * Automatically generate number of directories because enums start at 0 */
    N_GRID_DIRECTORIES
} Directory_Indices;

/** 
 * \brief Indices for the \ref grid_files array.
 * \ingroup GRID
 */
typedef enum
{
	/** \brief Location in \ref grid_files of the logfile name. */
    GRID_FILE_LOGFILE,
	/** \brief Location in \ref grid_files of the setup file name. */
    GRID_FILE_SETUP,
	/** \brief Location in \ref grid_files of the disturbance file name. */
    GRID_FILE_DISTURBANCES,
	/** \brief Location in \ref grid_files of the soil file name. */
    GRID_FILE_SOILS,
	/** \brief Location in \ref grid_files of the species spinup file
	 *         name. */
    GRID_FILE_SPINUP_SPECIES,
	/** \brief Location in \ref grid_files of the colonization file name. */
	GRID_FILE_COLONIZATION,
    /** \brief Location in \ref grid_files of the SOILWAT2 output cells file. */
    GRID_FILE_SOILWAT2_OUTPUT,

	/** \brief Location in \ref grid_files of the STEPWAT2 input file name. */
    GRID_FILE_FILES,
	/** \brief Location in \ref grid_files of the max rgroup and species file
	 *         name. */
    GRID_FILE_MAXRGROUPSPECIES,

	/** \brief Location in \ref grid_files of the biomass output file name. */
    GRID_FILE_PREFIX_BMASSAVG,
	/** \brief Location in \ref grid_files of the mortality output file 
	 *         name. */
    GRID_FILE_PREFIX_MORTAVG,
	/** \brief Location in \ref grid_files of the 
	 *         [seed dispersal](\ref SEED_DISPERSAL) output file name. */
    GRID_FILE_PREFIX_DISPERSALEVENTS,
	/** \brief Location in \ref grid_files of the average biomass output file
	 *         name. */
    GRID_FILE_PREFIX_BMASSCELLAVG,
	/** \brief Location in \ref grid_files of the average mortality file
	 *         name. */
    GRID_FILE_PREFIX_MORTCELLAVG,

    /* Automatically generate number of files since enums start at 0 */
    N_GRID_FILES
} File_Indices;

/**
 * \brief Enumerator for reading soil input parameters.
 * 
 * \ingroup GRID
 */
typedef enum
{
	/** \brief Successfully read soil input parameters. */
	SOIL_READ_SUCCESS = -2,
	/** \brief Failed to read soil input parameters. */
	SOIL_READ_FAILURE = -1,
} Soil_Read_Return_Values;

/************************ Exported Variable Declarations **************************/

/** 
 * \brief The main struct of the gridded mode module.
 * 
 * gridCells[i][j] denotes the cell at position (i,j)
 * 
 * \author Chandler Haukap
 * \ingroup GRID
 */
CellType** gridCells;

/** 
 * \brief Rows in the [grid](\ref gridCells).
 * \ingroup GRID
 */
int grid_Rows;

/** 
 * \brief Columns in the [grid](\ref gridCells).
 * \ingroup GRID
 */
int grid_Cols;

/** 
 * \brief Array of file names. Use the \ref File_Indices enum to pick the 
 *         correct index. 
 * \ingroup GRID
 */
char *grid_files[N_GRID_FILES];

/** \brief Array of directory names. Use the \ref Directory_Indices enum to 
 *         pick the correct index.
 * \ingroup GRID
 */
char *grid_directories[N_GRID_DIRECTORIES];
/** 
 * \brief TRUE if every cell should write its own output file.
 * \ingroup GRID
 */
Bool writeIndividualFiles;

/**
 * \brief Set to TRUE to tell the gridded mode to output SOILWAT2 files.
 * \ingroup GRID
 */
Bool writeSOILWAT2Output;

/**************************** Exported Functions **********************************/
/* See ST_grid.c for documentation of these functions. */
void runGrid(void);
void load_cell(int row, int col);
void unload_cell(void);
void rereadInputs(void);
void free_grid_memory(void);

#endif