/*********************************************************************************
 * ST_seedDispersal.c
 * 
 * Contains definitions for all functions belonging to the seed dispersal module.
 *********************************************************************************/

#include "ST_globals.h"
#include "ST_defines.h"
#include "ST_grid.h"
#include "ST_seedDispersal.h"
#include "rands.h"
#include "myMemory.h"

int _do_bulk_dispersal(SppIndex sp);
void _do_precise_dispersal(int leftoverSeeds, SppIndex sp);
float _cell_dist(int row1, int row2, int col1, int col2, float cellLen);

/* RNG unique to seed dispersal. */
pcg32_random_t dispersal_rng;

/* Derives the probabilities that a given cell will disperse seeds to any other cell. 
   results of this function can be accessed using 
   gridCells[a][b].mySeedDispersal[s].probabilityOfDispersing[c][d] 
   where (a,b) are the coordinates of the sender,
   (b, c) are the coordinates of the reciever,
   and s is the species. */
void initDispersalParameters(void)
{
	int senderRow, senderCol, recieverRow, recieverCol;
	SppIndex sp;
	double MAXD; /* maximum distance that a given species can disperse seeds */
	double maxRate; /* Dispersability of the seeds */
	double plotWidth; /* width of the plots (all plots are equal and square) */
	double distanceBetweenPlots; /* distance between the sender and the reciever */
    CellType* sender;

    RandSeed(SuperGlobals.randseed, &dispersal_rng);

	/* sender denotes that these loops refer to the cell distributing seeds */
	for(senderRow = 0; senderCol < grid_Rows; ++senderRow){
		for(senderCol = 0; senderCol < grid_Cols; ++senderCol){
			/* Cell is loaded to ensure the ForEachSpecies loop works */
			load_cell(senderRow, senderCol);
            sender = &gridCells[senderRow][senderCol];

			/* Allocate seed dispersal information. */
			sender->mySeedDispersal = Mem_Calloc(MAX_SPECIES, sizeof(Grid_SD_St), "initDispersalParameters");

			ForEachSpecies(sp){
				if (!(Species[sp]->use_me && Species[sp]->use_dispersal)){
					continue;
				}

				/* These are the three values we need to calculate the probability of dispersal */
				MAXD = ((Species[sp]->sd_H * Species[sp]->sd_VW) / Species[sp]->sd_VT) / 100.0; // divided by 100 to convert from cm to m.
				maxRate = -(log(0.005) / MAXD);
				plotWidth = sqrt(Globals->plotsize);

				/* Allocate the probabilityOfDispersing 2d array */
				sender->mySeedDispersal[sp].probabilityOfDispersing = Mem_Calloc(grid_Rows, 
							sizeof(double*), "initDispersalParameters: probabilityOfDispersing");
				for(recieverRow = 0; recieverRow < grid_Rows; ++recieverRow){
					sender->mySeedDispersal[sp].probabilityOfDispersing[recieverRow] = 
								Mem_Calloc(grid_Cols, sizeof(double), "initDispersalParameters: probabilityOfDispersing[i]");
				}

				/* Loop through all possible recipients of seeds. */
				for(recieverRow = 0; recieverRow < grid_Rows; ++recieverRow){
					for(recieverCol = 0; recieverCol < grid_Cols; ++recieverCol){
						if(senderRow == recieverRow && senderCol == recieverCol){
							continue; // No need to calculate a probability for dispersal to itself
						}
                        
						distanceBetweenPlots = _cell_dist(senderRow, recieverRow, senderCol, recieverCol, plotWidth);

                        /* The value that we are after should be saved to the sender cell. */
						sender->mySeedDispersal[sp].probabilityOfDispersing[recieverRow][recieverCol]
                                    = (distanceBetweenPlots > MAXD) ? (0.0) : (exp(-maxRate * distanceBetweenPlots));
					}
				}
			}
		}
	}
	unload_cell();
}

/* Perform seed dispersal durring the simulation. This is NOT functional yet. */
void disperseSeeds(void)
{
	/************ TODO: overhaul seed dispersal. This block prevents seed dispersal from running. **************/
	printf("\nSeed dispersal during the simulation is not yet functional.\n"
		    "Check out GitHub for updates on this feature.\n");
	UseSeedDispersal = FALSE;
	return;
	/***********************************************************************************************************/

/* TODO: discuss this implementation:
    SppIndex sp;
    ForEachSpecies(sp){
        int leftovers = _do_bulk_dispersal(sp);
        _do_precise_dispersal(leftovers, sp);
    }
*/

/* 
	float biomass, randomN, LYPPT, presentProb, receivedProb;
	int i, j, germ, sgerm, year, row, col;
	SppIndex s;
	CellType *cell;

    // Load the first cell so we can access seedling_estab_prob and currYear,
    // which are not specific to each cell.
    load_cell(0, 0);

	if (Globals->currYear == 1)
	{ //since we have no previous data to go off of, use the current years...
		for (i = 0; i < grid_Cells; i++)
		{
		    row = i / grid_Cols;
		    col = i % grid_Cols;
		    cell = &gridCells[row][col];

		    load_cell(row, col);

            ForEachSpecies(s)
            {
                if (!(Species[s]->use_me && Species[s]->use_dispersal))
                    continue;
                Species[s]->allow_growth = Species[s]->sd_sgerm =
                        1;// since it's the first year, we have to allow growth...
                if (UseDisturbances)
                    // RGroup[x]->killyr is the same for all x
                    if (1 == RGroup[0]->killyr)
                        Species[s]->allow_growth = 0;
                cell->mySeedDispersal[s].lyppt = Env->ppt;
            }
        }
	}
	else
	{
		// figure out whether or not to allow growth for the current year... based upon whether the species already has plants or germination allowed this year and seeds received last year...
		ForEachSpecies(s)
		{

			if (!(Species[s]->use_me && Species[s]->use_dispersal))
				continue;

			// germination probability
			randomN = RandUni(&grid_rng);
			germ = LE(randomN, Species[s]->seedling_estab_prob);

			year = Globals->currYear - 1;

			for (i = 0; i < grid_Cells; i++)
			{
                row = i / grid_Cols;
                col = i % grid_Cols;
                cell = &gridCells[row][col];

                load_cell(row, col);

				if (Globals->currYear <= SuperGlobals.runInitializationYears)
				{
					cell->mySeedDispersal[s].seeds_present = 1;
				}
				else if (Globals->currYear <= SuperGlobals.runInitializationYears
						 && cell->mySpeciesInit.shouldBeInitialized[s])
				{
					cell->mySeedDispersal[s].seeds_present = 1;
				}

				sgerm = (cell->mySeedDispersal[s].seeds_present
						|| cell->mySeedDispersal[s].seeds_received) && germ; //refers to whether the species has seeds available from the previous year and conditions are correct for germination this year
				Species[s]->allow_growth = FALSE;
				biomass = getSpeciesRelsize(s)
						* Species[s]->mature_biomass;

				if (UseDisturbances)
				{
					if ((sgerm || year < RGroup[0]->killyr
							|| RGroup[0]->killyr <= 0 || GT(biomass, 0.0))
					&& (year != RGroup[0].killyr))
					{
						//commented above one condition as it was causing a bug there, next year of killing year will make
						//allow_growth flag to false as 	year = Globals->currYear - 1 , so for example if killing year= 6 and Globals->currYear=7 then here
						// year variable will be 7-1 =6 that is equal to killing year 6, so this condition (year != RGroup[0].killyr)
						//will fail and allow_growth will not become TRUE, then when Globals.currYear=8 this allow_growth= FALSE will carry forward and there will no call
						// to other functions. Last year size will carry forward so in final output year 7 and year 8 will
						// have same output that is not correct.
						Species[s]->allow_growth = TRUE;
					}

				}
				else if (sgerm || GT(biomass, 0.0))
					Species[s]->allow_growth = TRUE;
				Species[s]->sd_sgerm = sgerm; //based upon whether we have received/produced seeds that germinated
			}
		}

	}

	// calculate whether or not seeds were received/produced this year, this data is used the next time the function is called
	ForEachSpecies(s)
	{
		if (!(Species[s]->use_me && Species[s]->use_dispersal))
			continue;

		IndivType* indiv;

		// figure out which species in each cell produced seeds...
		for (i = 0; i < grid_Cells; i++)
		{
            row = i / grid_Cols;
            col = i % grid_Cols;
            cell = &gridCells[row][col];

            load_cell(row, col);

			cell->mySeedDispersal[s].seeds_present = cell->mySeedDispersal[s].seeds_received =
					Species[s]->received_prob = 0;

			biomass = 0;	//getting the biggest individual in the species...
			ForEachIndiv(indiv, Species[s])
				if (indiv->relsize * Species[s]->mature_biomass
						> biomass)
					biomass = indiv->relsize
							* Species[s]->mature_biomass;

			if (GE(biomass,
					Species[s]->mature_biomass
							* Species[s]->sd_Param1))
			{
				randomN = RandUni(&grid_rng);

				LYPPT = cell->mySeedDispersal[s].lyppt;
				float PPTdry = Species[s]->sd_PPTdry, PPTwet =
						Species[s]->sd_PPTwet;
				float Pmin = Species[s]->sd_Pmin, Pmax =
						Species[s]->sd_Pmax;

				//p3 = Pmin, if LYPPT < PPTdry
				//p3 = 1 - (1-Pmin) * exp(-d * (LYPPT - PPTdry)) with d = - ln((1 - Pmax)/(1 - Pmin)) / (PPTwet - PPTdry), if PPTdry <= LYPPT <= PPTwet
				//p3 = Pmax, if LYPPT > PPTwet

				presentProb = 0.0;
				if (PPTdry <= LYPPT && LYPPT <= PPTwet)
				{
					float d = -log(((1 - Pmax) / (1 - Pmin)))
							/ (PPTwet - PPTdry); //log is the natural log in STD c's math.h
					presentProb = 1 - (1 - Pmin) * exp((-d * (LYPPT - PPTdry)));
				}
				else if (LYPPT < PPTdry)
					presentProb = Pmin;
				else if (LYPPT > PPTwet)
					presentProb = Pmax;

				if (LE(randomN, presentProb))
					cell->mySeedDispersal[s].seeds_present = 1;
			}
		}

		// figure out which species in each cell received seeds...
		for (i = 0; i < grid_Cells; i++)
		{
            row = i / grid_Cols;
            col = i % grid_Cols;
            cell = &gridCells[row][col];

            load_cell(row, col);

			if (cell->mySeedDispersal[s].seeds_present)
				continue;
			receivedProb = 0;

			for (j = 0; j < cell->mySeedDispersal[s].size; j++)
				if (cell->mySeedDispersal[cell->mySeedDispersal[s].cells[j]].seeds_present)
					receivedProb += cell->mySeedDispersal[s].probabilityOfDispersing[j];

			randomN = RandUni(&grid_rng);
			if (LE(randomN, receivedProb) && !ZRO(receivedProb))
				cell->mySeedDispersal[s].seeds_received = 1;
			else
				cell->mySeedDispersal[s].seeds_received = 0;

			Species[s]->received_prob = receivedProb;
		}
	}

	unload_cell();
*/
}

int _do_bulk_dispersal(SppIndex sp){
    return 0;
}

void _do_precise_dispersal(int leftoverSeeds, SppIndex sp){
    return;
}

/* Returns the distance in meters between two cells. 
   To calculate the distance between cells (a,b) and (c,d) with plot width w,
   input _cell_dist(a, c, b, d, w) */
float _cell_dist(int row1, int row2, int col1, int col2, float cellLen)
{
    double rowDist = row1 - row2;
    double colDist = col1 - col2;

	//returns the distance between the two grid cells
	if (row1 == row2)
	{
		return (abs(colDist) * cellLen);
	}
	else if (col1 == col2)
	{
		return (abs(rowDist) * cellLen);
	}
	else
	{ // row1 != row2 && col1 != col2
		//the problem can be thought of in terms of a right triangle...
		//using the pythagorean theorem: c = sqrt(a^2 + b^2)... c (the hypotenuse) represents the distance that we need.  a is the distance between columns and b is the distance between rows.
		return sqrt(pow(abs(colDist)*cellLen, 2.0) + pow(abs(rowDist)*cellLen, 2.0));
	}
}