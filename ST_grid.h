/******************************************************************/
/* ST_grid.h
    This header defines global functions, variables, structs and 
    enumerators from ST_grid.c. If you are looking for the functions 
    that this module exports check the bottom of the file.

    Initial programming by Chandler Haukap in August 2019. */
/******************************************************************/

#ifndef GRID_H
#define GRID_H

/******** These modules are necessary to compile ST_grid.c ********/
#include "ST_stats.h"
#include "ST_defines.h"
#include "sw_src/pcg/pcg_basic.h"
#include "sxw_vars.h"
#include "SW_Site.h"
#include "SW_SoilWater.h"
#include "SW_VegProd.h"
#include "SW_Model.h"
#include "SW_Weather.h"
#include "ST_seedDispersal.h"

/*********************** Grid Structures ****************************/

// Contains the input data for all the soil layers of a cell
struct Soil_st
{
	// Number of soil layers (size of lyr array)
	int num_layers;
	// Name of the roots file belonging to this cell
	char rootsFile[20];
	// Specific layer's information. The size of these arrays is the
    // same as num_layers.
    RealF *depth;       /* Each Layer's max depth */
    RealF *matricd;     /* Each Layer's matricx */
    RealF *gravel;      /* Each Layer's gravel content */
    RealF *evco;        /* Each Layer's evco */
    RealF *trco_grass;  /* Each Layer's trco for grass */
    RealF *trco_shrub;  /* Each Layer's trco for shrubs */
    RealF *trco_tree;   /* Each Layer's trco for tree */
    RealF *trco_forb;   /* Each Layer's trco for forbs */
    RealF *psand;       /* Each Layer's sand content */
    RealF *pclay;       /* Each Layer's clay content */
    RealF *imperm;      /* Each layer's impermiability */
    RealF *soiltemp;    /* Each Layer's temperature */
}typedef SoilType;

/* Initialization information. */
struct grid_init_species_st
{
	/* TRUE if at least one species has requested initialization */
	int useInitialization;
	/* Array of Boolean values. TRUE if given species
	   should be included in spinup */
	int *shouldBeInitialized;
}typedef Grid_Init_Species_St;

/* Holds all cell-specific parameters */
struct grid_cell_st
{
	/* RGroup corresponding to this cell */
	GroupType **myGroup;
	/* Species corresponding to this cell */
	SpeciesType **mySpecies;
	/* Succulents corresponding to this cell */
	SucculentType mySucculent;
	/* This cell's environment. We expect each cell to
	 * have slightly different weather each year */
	EnvType myEnvironment;
	/* Cell's plot data */
	PlotType myPlot;
	/* Global variables corresponding to this cell */ 
	ModelType myGlobals;
	/* If TRUE this cell should use seed dispersal */
	Bool useSeedDispersal;
	/* TRUE if this cell is in spinup mode */
	Bool DuringInitialization;
	/* species spinup information */
	Grid_Init_Species_St mySpeciesInit;

	Bool* someKillage;
	
	/* ---------------- accumulators -------------------- */
	StatType *_Dist, *_Ppt, *_Temp,
  		*_Grp, *_Gsize, *_Gpr, *_Gmort, *_Gestab,
  		*_Spp, *_Indv, *_Smort, *_Sestab, *_Sreceived;
	FireStatsType *_Gwf;
	Bool stats_init;
	/* -------------- end accumulators ------------------ */

	/* -------------------- SXW ------------------------- */
	transp_t* myTranspWindow;
	SXW_t* mySXW;
	SXW_resourceType* mySXWResources;
	/* ------------------ End SXW ----------------------- */

	/* ------------------- Soils ------------------------ */
	// Soil layer information for this cell.
	SoilType mySoils;
	/* ------------------ End Soils --------------------- */
} typedef CellType;

/**************************** Enumerators *********************************/
/* Indices for grid_directories go here */
typedef enum
{
    GRID_DIRECTORY_STEPWAT_INPUTS,

    /* Automatically generate number of directories since enums start at 0 */
    N_GRID_DIRECTORIES
} Directory_Indices;

/* Indices for grid_files go here */
typedef enum
{
    GRID_FILE_LOGFILE,
    GRID_FILE_SETUP,
    GRID_FILE_DISTURBANCES,
    GRID_FILE_SOILS,
    GRID_FILE_INIT_SPECIES,
    GRID_FILE_FILES,
    GRID_FILE_MAXRGROUPSPECIES,

    GRID_FILE_PREFIX_BMASSAVG,
    GRID_FILE_PREFIX_MORTAVG,
    GRID_FILE_PREFIX_RECEIVEDPROB,
    GRID_FILE_PREFIX_BMASSCELLAVG,
    GRID_FILE_PREFIX_MORTCELLAVG,

    /* Automatically generate number of files since enums start at 0 */
    N_GRID_FILES
} File_Indices;

typedef enum
{
	SOIL_READ_SUCCESS = -2,
	SOIL_READ_FAILURE = -1,
} Soil_Read_Return_Values;

/************************ Exported Variable Declarations **************************/

/* gridCells[i][j] denotes the cell at position (i,j) */
CellType** gridCells;
/* Rows in the grid */
int grid_Rows;
/* Columns in the grid */
int grid_Cols;
/* Array of file names. Use the File_Indices enum to pick the correct index. */
char *grid_files[N_GRID_FILES];
/* Array of directory names. Use the Directory_Indices enum to pick the correct index. */
char *grid_directories[N_GRID_DIRECTORIES];
/* TRUE if every cell should write its own output file. */
Bool writeIndividualFiles;

/**************************** Exported Functions **********************************/

void runGrid(void);
void load_cell(int row, int col);
void unload_cell(void);
void rereadInputs(void);
void free_grid_memory(void);

#endif