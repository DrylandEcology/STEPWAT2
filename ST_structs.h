/** 
 * \file ST_structs.h
 * \brief Contains the definitions of STEPWAT structures.
 * 
 * These structures include individuals, species, resource groups, succulents,
 * plots, environments, and global variables.
 * For SXW struct definitions see \ref sxw.h.
 * 
 * A note on annuals:
 * Basically, the number of individuals is determined by the available resources 
 * and a seedbank mechanism, and the growth rate is simply the proportion of max
 * size (1.0) each individual can achieve in the current year, computed by
 * maxbio * irate / PR. This is not the current implementation, but annuals are
 * currently being reevaluated so keep it in mind when modifying the code.
 *  
 * \author Chandler Haukap
 * 
 * \date 21 August 2019
 * 
 * \ingroup STEPPE
 */

#ifndef STEPPE_STRUCT_DEF
#define STEPPE_STRUCT_DEF

#include "sw_src/generic.h"
#include "ST_mortality.h"

/**
 * \brief Holds information on perennial plant individuals. 
 * 
 * Most of the code associated with this struct is found in \ref ST_indivs.c.
 * NOTE: If you add a variable to this struct, make sure you add code to \ref copy_individual()
 * 
 * \sa indiv_ann_st
 * 
 * \ingroup INDIVIDUAL
 */
struct indiv_st {
      /** \brief unique identifier for each individual. */
  IntUS id,
      /** \brief age in years of this individual. */
  	  age,
      /** \brief millimeters of extra resources. */
      mm_extra_res,
      /** \brief number of years this individual has slow growth */
      slow_yrs,
      /** \brief Species[indv->myspecies] is the species this individual belongs to.
       * \sa Species  */
      myspecies,
      /** \brief num yrs pr > 1 (stretched resources) */
      yrs_neg_pr;
      /** \brief What killed this individual.
       * \sa MortalityType */
  MortalityType killedby;
      /** \brief only for clonal plants, means veggrow possible. */
  Bool killed;
      /** \brief Size relative to a full-sized individual. */
  RealF relsize,
       /** \brief Previous year's relsize. Used for proportional recovery calculation.
        * \sa proportion_Recovery() */
       prv_yr_relsize,
       /** \brief proportion of total group relsize belonging to this individual. */
       grp_res_prop,
       /** \brief resources required, biomass of the individual. 
        * \sa rgroup_ResPartIndiv */
       res_required,
       /** \brief Resources available. */
       res_avail,
       /** \brief Resources applied to superfluous growth. */
       res_extra,
       /** \brief Ratio of resources required to amount available. */
       pr,
       /** \brief Actual growth rate .
        * \sa rgroup_Grow() */
       growthrate,
       /** \brief biomass the plant gained this year excluding superfluous biomass. Used for grazing. 
        * \sa grazing_EndOfYear() */
       normal_growth,
       /** \brief set when killed; 0 if not clonal.
        * \sa indiv_Kill_Partial() */
       prob_veggrow;
       /** \brief Allows for a doubly-linked list of individuals. Implemented in \ref Species.
        * \sa Species */
  struct indiv_st *Next, 
       /** \brief Allows for a doubly-linked list of individuals. Implemented in \ref Species.
        * \sa Species */
                  *Prev;
};

/** 
 * \brief Holds information on annual plant individuals.
 * 
 * Most of the code associated with this struct is found in \ref ST_indivs.c.
 * 
 * \sa indiv_st
 * 
 * \ingroup INDIVIDUAL
 */
struct indiv_ann_st {
      /** \brief millimeters of extra resources. */
  IntUS mm_extra_res,
      /** \brief Species[indv->myspecies] is the species this individual belongs to.
       * \sa Species  */
      myspecies;
      /** \brief Size relative to a full-sized individual. */
  RealF relsize,
       /** \brief proportion of total group relsize belonging to this individual. */
       grp_res_prop,
       /** \brief resources required, biomass of the individual. 
        * \sa rgroup_ResPartIndiv */
       res_required,
       /** \brief Resources available. */
       res_avail,
       /** \brief Resources applied to superfluous growth. */
       res_extra,
       /** \brief Ratio of resources required to amt available. */
       pr,
       /** \brief Actual growth rate .
        * \sa rgroup_Grow() */
       growthrate;
        /** \brief Allows for a doubly-linked list of individuals. Implemented in Species.
         * \sa Species */
  struct indiv_ann_st *Next, 
        /** \brief Allows for a doubly-linked list of individuals. Implemented in Species.
         * \sa Species */
                  *Prev;
};

/** 
 * \brief Holds all species-specific information.
 * 
 * This struct contains constants and variables that differ between
 * species. \ref Species is a global array of these structs.
 * 
 * NOTE: If you add a variable to this struct, make sure you add code to \ref copy_species()
 * 
 * \ingroup SPECIES
 */
struct species_st {

  /**** Quantities that can change during model runs *****/

      /** \brief Number of individuals established (growing) */
  SppIndex est_count;
      /** \brief Pointer to array of individuals killed, indexed by age at death. */
  IntUS *kills,
      /** \brief Number of individuals established in current iteration. */
        estabs,
      /** \brief Specific to annuals. Array of seeds produced indexed by year. 
       * Not to be confused with [seed dispersal](\ref SEED_DISPERSAL) */
        *seedprod,
      /** \brief not sure what this does yet.
       * This is NOT a part of [seed dispersal](\ref SEED_DISPERSAL) */
         seedbank,
      /** \brief Average number of seeds produced by annual species per 1g of biomass, per 1m^2 and per year.
       * internally re-calculated as seeds per 1 g biomass per plot and per year. */
         pseed;
      /** \brief relsize from the previous year, used for annual establishment. 
       * \sa rgroup_Establish() */
  RealF lastyear_relsize,
      /** \brief Amount of superfluous growth from extra resources. */
        extragrowth,
      /** \brief Probability that this species received seeds this year. */
	      received_prob,
      /** \brief Alpha parameter for random number draw from beta distribution in establishment of annual species. 
       * \sa _add_annuals() */
        alpha,
      /** \brief Beta parameter for random number draw from beta distribution in establishment of annual species.
       * \sa _add_annuals() */
        beta;
      /** \brief Variance parameter of the beta distribution for establishment of annual species. */
  float var;
      /** \brief Head of a doubly-linked list of all individuals of this species. 
       * \sa indiv_st*/
  struct indiv_st *IndvHead;
      /** \brief Whether seeds where produced/received and germinated this year. 
       * \sa ST_seedDispersal.c 
       * \ingroup SEED_DISPERSAL */
  Bool seedsPresent;

  /**** Quantities that DO NOT change during model runs *****/

      /** \brief 4-letter code (plus \0) for genus and species. */
  char* name;
      /** \brief Max age of mature plant. Also used to determine annual species. */
  IntUS max_age,
      /** \brief Annuals: max years of viability of seeds. */
        viable_yrs,
      /** \brief Max seedlings that can establish in 1 year. */
        max_seed_estab,
      /** \brief Max vegetative regrowth units, for example tillers. */
        max_vegunits,
      /** \brief Years slow growth is allowed before mortality. */
        max_slow,
      /** \brief Index of this species in the Species variable.
       * \sa Species */
        sp_num,
      /** \brief Resource group that this species belongs to. 
       * \sa RGroup */
        res_grp;
      /** \brief Max growth rate. Defined as intrin_rate * proportion. */
  RealF max_rate,
      /** \brief Intrinsic growth rate. */
        intrin_rate,
      /** \brief Starting size of a new individual of this species, relative to an adult plant. */
        relseedlingsize,
      /** \brief Seedling biomass in grams of an individual of this species. */
        seedling_biomass,
      /** \brief Maximum biomass in grams an individual of this species can achieve. */
        mature_biomass,
      /** \brief Holds seedling establishment probability directly from inputs. 
       * \sa seedling_estab_prob */
        seedling_estab_prob_old,
      /** \brief Seedling establishment probability.
       * \sa seedling_estab_prob_old */
        seedling_estab_prob,
      /** \brief Annual mortality probability. This variable is currently unused. */
        ann_mort_prob,
      /** \brief Used in age-independent mortality.
       * \sa _age_independent() */
        cohort_surv,
      /** \brief Annuals-specific exponent for viability decay function. */
        exp_decay,
      /** \brief 1 value for each mortality type, if clonal*/
        prob_veggrow[4],
      /** \brief Minimum relsize at which an individual can produce seeds.
       * \ingroup SEED_DISPERSAL */
  	    minReproductiveSize,
      /** \brief The maximum height of an individual of this species.
       * \ingroup SEED_DISPERSAL */
        maxHeight,
      /** \brief The slope of the allometric relationship between height and 
       *         biomass. Constant read in from inputs.
       * \ingroup SEED_DISPERSAL */
        heightSlope,
      /** \brief The probability that a species disperses seeds 
       *         \ref maxDispersalDistance meters away. 
       *  \ingroup SEED_DISPERSAL */
        maxDispersalProbability;
      /** \brief Temperature class for this species.
       * \sa TempClass */
  TempClass tempclass;
      /** \brief Disturbance class for this species.
       * \sa DisturbClass */
  DisturbClass disturbclass;
      /** \brief TRUE if this species is clonal. */
  Bool isclonal,
      /** \brief TRUE if this species should respond to temperature. Not currently implemented. */
       use_temp_response,
      /** \brief If FALSE do not establish this species. */
       use_me,
      /** \brief TRUE if the user has requested seed dispersal for this species.
       * \sa ST_seedDispersal.c 
       * \ingroup SEED_DISPERSAL */
       use_dispersal;
};

/** 
 * \brief Contains all resource group-specific fields. 
 * 
 * Some fields are constant and some are variable, so be careful what you are modifying. 
 * \ref RGroup is a global array of these structs used in STEPWAT2.
 * 
 * NOTE: If you add a variable to this struct, make sure you add code to \ref copy_rgroup()
 * 
 * \ingroup RGROUP
 */
struct resourcegroup_st {

  /**** Quantities that can change during model runs *****/

      /** \brief individuals in group killed. Index by age at death. */
  IntUS *kills,
      /** \brief Total indivs in group established during current iteration. */
        estabs,
      /** \brief If killyr equals the current year there will be a fire.
       * \sa killfreq
       * \sa mort_EndOfYear() */
        killyr,
      /** \brief Accumulator for consecutive years of low resources. */
        yrs_neg_pr,
      /** \brief Treated as a Bool. If TRUE there was a wildfire in the current year. 
       * \sa mort_EndOfYear() */
        wildfire,
      /** \brief Treated as a Bool. If TRUE there was a prescribed fire in the current year. 
       * \sa mort_EndOfYear() */
        prescribedfire,
      /** \brief Extra resources in millimeters */
        mm_extra_res;
      /** \brief Resources required for current biomass of resgroup. */
  RealF res_required,
      /** \brief Resources available adjusted to account for competition. */
        res_avail,
      /** \brief Resource applied to superficial growth. */
        res_extra,
      /** \brief Ratio of resources required to resources available */
        pr,
      /** \brief Fraction of the total biomass for this vegtype contributed by this resgroup.
       * \sa veg_prod_type */
        rgroupFractionOfVegTypeBiomass,
      /** \brief ratio of biomass/m2 / transp/m2.
       * Used in \ref sxw.c */
        _bvt;
      /** \brief Number of species established in group. */
  SppIndex est_count,
      /** \brief Array of indexes of species established in this resgroup. 
       * \sa Species */
           *est_spp;
      /** \brief If TRUE this group is extirpated, no more regeneration. */
  Bool extirpated,
      /** \brief For annuals: If TRUE this group can regenerate. */
       regen_ok;

  /**** Quantities that DO NOT change during model runs *****/

      /** \brief Number of years resources can be stretched without killing the group. */
  IntUS max_stretch,
      /** \brief Max number of species that can add new plants per year*/
        max_spp_estab,
      /** \brief Number of species in the group*/
        max_spp,
      /** \brief Longest lifespan in resgroup. Used to malloc kills[] */
        max_age,
      /** \brief Don't start trying to grow until this year. */
        startyr,
      /** \brief Start year for kill frequency. 
       * \sa killfreq */
	      killfreq_startyr,
      /** \brief Year in which group is extirpated. Setting to 0 turns off extirpation. */
        extirp,
      /** \brief Index number of this group in the RGroup array.
       * \sa RGroup */
        grp_num,
      /** \brief SOILWAT2 vegetation type:
       * STEPWAT2 inputs via "rgroup.in" are defined as: 1 for tree, 2 for shrub, 3 for grass, 4 for forb;
       * however, the inputs get translated by get_SW2_veg_index() to SOILWAT2 values immediately upon reading the inputs. 
       * \sa SW_Defines.h */
        veg_prod_type,
      /** \brief Start year for grazing frequency.
       * \sa grazingfrq */
        grazingfreq_startyr;
      /** \brief Array of species indexes belonging to this group.
       * \sa Species */
  SppIndex *species;
      /** \brief Input from table. */
  RealF space,
      /** \brief input space from table, rescaled if one or more rgroups is not established. */ 
        min_res_req,
      /** \brief density of mature plants in units of plants / m^2 */
        max_per_sqm,
      /** \brief Perform grazing on this group at this frequency. */
		grazingfrq,
      /** \brief Number of mature plants allowed per plot. */
        max_density,
      /** \brief Max biomass of group. */
        max_bmass,
      /** \brief Kill group at this frequency.
       * Values < 1 result in probablistic fires
       * values > 1 result in fires at a return interval = killfreq. 
       * \sa mort_EndOfYear() */
        killfreq,
      /** \brief Ephemeral growth. Used to calculate superficial growth when precipitation is high. */
        xgrow,
      /** \brief Growthrate that triggers mortality. Defined in inputs. */
        slowrate,
      /** \brief Space equation: slope for wet/dry/normal years */
        ppt_slope[3],
      /** \brief Space equation: intercept for wet/dry/normal years*/
        ppt_intcpt[3],
      /** \brief Proportion of group killed. */
		    proportion_killed,
      /** \brief Proportion recovery after killing year */
		    proportion_recovered,
      /** \brief Proportion of biomass removed during grazing year.
       * \sa grazingfrq */
        proportion_grazing;
      /** \brief TRUE if this group is comprised of succulents. */
  Bool succulent,
      /** \brief IF TRUE this group responds to other group's unused resources. */
       use_extra_res,
      /** \brief If FALSE establish no species of this group. */
       use_me,
      /** \brief If TRUE use age-independent and slowgrowth mortality. */
       use_mort,
      /** \brief Establish this group every year if TRUE. */
       est_annually;
      /** \brief Rooting depth class.
       * \sa DepthClass */
  DepthClass depth;
      /** \brief name of this group, specified in inputs. */
  char *name;
};

/** 
 * \brief Succulent specific constants read from or derived from inputs.
 * 
 * A global instance of this struct is Succulent.
 * 
 * \sa Succulent.
 * 
 * \ingroup SUCCULENTS
 */
struct succulent_st {
      /** \brief Growth modifier parameters for succulents (eqn 10). */
  RealF growth[2],
      /** \brief Mortality modifier parameters for succulents (eqn 16). */
        mort[2],
      /** \brief If not killed, reduce by eqn 10. */
        reduction,
      /** \brief Calculated from eqn 16. */
        prob_death;
};

/**
 * \brief Contains variables used in generating the environment.
 * 
 * This struct is instanciated by the global variable Env.
 * 
 * \sa Env_Generate()
 * 
 * \ingroup STEPPE
 */
struct environs_st {
      /** \brief Is this site wet, dry, or normal?.
       * \sa PPTClass */
  PPTClass wet_dry;
      /** \brief Precipitation for the year in millimeters. */
  IntS ppt,
      /** \brief Precipitation for the previous (last) year in millimeters. */
       lyppt,
      /** \brief Precipitation during the growing season in millimeters. */
       gsppt;
      /** \brief Average daily temp for the year in celsius. */
  RealF temp,
      /** \brief Amount to reduce growth by temperature.
       * (eqns 12,13), one for each tempclass. */
       temp_reduction[2];
};

/**
 * \brief Contains variables used in generating the plot.
 * 
 * plot_st is instanciated by the global variable Plot.
 * 
 * \sa Plot
 * 
 * \ingroup STEPPE
 */
struct plot_st {
  /** \brief Type of disturbance.
   * \sa DisturbanceEvent */
  DisturbEvent disturbance;
  /** \brief If TRUE fecal pats can be removed. */
  Bool pat_removed;
  /** \brief Years remaining before recolonization (ie, new establishments) can begin again.
   * Otherwise, if disturbance is fecalpat, number of years it has been ongoing. Set to 0
   * when the disturbance effect is expired. */
  IntUS disturbed;
};

/**
 * \brief Stores precipitation information.
 * 
 * Used inside the globals_st struct.
 * 
 * \sa globals_st
 * 
 * \ingroup STEPPE
 */
struct ppt_st {
      /** \brief Average precipitation. */
  RealF avg,
      /** \brief Standard deviation of precipitation. */
       std;
      /** \brief Minimum precipitation. */
  IntUS  min,
      /** \brief Maximum precipitation. */
       max,
      /** \brief Number of dry days. */
       dry,
      /** \brief Number of wet days. */
       wet;
};

/**
 * \brief Stores temperature information.
 * 
 * Used in the globals_st struct.
 * 
 * \sa globals_st
 * 
 * \ingroup STEPPE
 */
struct temp_st {
  RealF avg,
        std,
        min,
        max,
        gstemp;
};

/** 
 * \brief Fecal pat information.
 * 
 * Used in the globals_st struct.
 * 
 * \sa globals_st
 * 
 * \ingroup MORTALITY
 */
struct fecalpats_st {
      /** \brief If TRUE fecal pats will be used. */
  Bool use;
      /** \brief Between 0 and 1. How often pats occur. */
  RealF occur,
      /** \brief Between 0 and 1. How likely pats are removed. */
        removal,
      /** \brief 1 elem. for slope and intercept. */
        recol[2];
};

/**
 * \brief Ant mound information.
 * 
 * Used in the globals_st struct.
 * 
 * \sa globals_st
 * 
 * \ingroup MORTALITY
 */
struct antmounds_st {
      /** \brief If TRUE ant mounds will be used. */ 
  Bool use;
      /** \brief The probability of ant mound occurrence. Between 0 and 1. */
  RealF occur;
      /** \brief The minimum number of ant mounds to generate in a year. */
  IntUS minyr,
      /** \brief The maximum number of ant mounds to generate in a year. */
       maxyr;
};

/**
 * \brief Holds information on animal burrows.
 * 
 * Used in the globals_st struct.
 * 
 * \sa globals_st
 * 
 * \ingroup MORTALITY
 */
struct burrows_st {
  /** \brief If TRUE burrows will be used. */
  Bool use;
  /** \brief Value between 0 and 1. The probability that a burrow will occur. */
  RealF occur;
  /** \brief Minimum number of burrows in a year. */
  IntUS minyr;
};

/**
 * \brief Stores file names for output.
 * 
 * The struct contains the name of the yearly and summary output files.
 * Used in the globals_st struct.
 * 
 * \sa globals_st
 * 
 * \ingroup STEPPE
 */
struct outfiles_st {
      /** \brief File handle for yearly output. */
  FILE *fp_year,
      /**\brief File handle for averages output. */
       *fp_sumry;
      /** \brief max width of outfile suffix if printing yearly. */
  IntUS suffixwidth;
};

/**
 * \brief Stores global information that is independent of any species or group.
 * 
 * globals_st deals mainly with structural variables like year and iteration. It also
 * contains information on file names, precipitation, temperature, mounds, pats and burrows.
 * This struct is instanciated by the global variable Globals.
 * 
 * \sa Globals
 * 
 * \ingroup STEPPE
 */
struct globals_st {
  /** \brief Precipitation constants. \sa ppt_st */
  struct ppt_st ppt;
  /** \brief Temperature constants. \sa temp_st */
  struct temp_st temp;
  /** \brief Fecal pat constants. \sa fecalpats_st */
  struct fecalpats_st pat;
  /** \brief Ant mound constants. \sa antmounds_st */
  struct antmounds_st mound;
  /** \brief Animal burrow constants. \sa burrows_st */
  struct burrows_st burrow;

      /** \brief Size of the plot in square meters */
  RealF plotsize,
      /** \brief Proportion of ppt during growing season. */
        gsppt_prop,
      /** \brief warm-season/cool-season species growth modifiers. */
        tempparm[2][3];
      /** \brief Number of years to run the model. */
  IntUS runModelYears,
      /** \brief Oldest plant; same as runModelYears for now. \sa runModelYears. */
      Max_Age,
      /** \brief The year that the simulation is currently in. */
      currYear,
      /** \brief Number of iterations to run the simulation. */
      runModelIterations,
      /** \brief The iteration that the simulation is currently in. */
      currIter,
      /** \brief Number of groups defined. */
      grpCount,
      /** \brief Number of species defined*/
      sppCount,
      /** \brief Number of years for which transpiration data is kept. \sa _transp_contribution_by_group() */
      transp_window,
      /** \brief Number of cells to use in Grid, only applicable if gridded mode is being used. */
      nCells;
      /** \brief Random seed from input file. Used to seed the PCG32 RNGs. */
  IntL randseed;
      /** \brief Maximum resource groups allowed. */
  size_t max_rgroups,
      /** \brief Maximum resource group name length. */
         max_groupnamelen,
      /** \brief Maximum species allowed per resource group. */
         max_spp_per_grp,
      /** \brief Maximum individuals allowed per species. */
         max_indivs_per_spp,
      /** \brief Maximum species name length. */
         max_speciesnamelen;

                    /** \brief Output file names for bmass files. */
  struct outfiles_st bmass, 
                    /** \brief Output file names for mort files. */
                     mort;
};

/**
 * \brief Flags for what biomass files we would like output.
 * 
 * These flags are specified in inputs and used heavily in ST_stats.c.
 * They are instanciated globally by BmassFlags.
 * 
 * \sa BmassFlags
 * 
 * \ingroup STEPPE
 */
struct bmassflags_st {
      /** \brief If FALSE print no biomass output. */
  Bool summary,
      /** \brief If TRUE print biomass output for each year of the simulation. */
       yearly,
      /** \brief If TRUE the output files need a header. */
       header,
      /** \brief If TRUE the output files will contain a year column. */
       yr,
      /** \brief If TRUE Output disturbance information. */
       dist,
      /** \brief If TRUE Output precipitation information. */
       ppt,
      /** \brief If TRUE Output precipitation class information. */
       pclass,
      /** \brief If TRUE output the yearly average temperature. */
       tmp,
      /** \brief If TRUE output the group biomasses. */
       grpb,
      /** \brief If TRUE output relsize. */
       pr,
      /** \brief If TRUE output biomass. */
       size,
      /** \brief If TRUE output species information. */
       sppb,
      /** \brief Print total wild fire count across all iterations. */
       wildfire,
      /** \brief Print prescribed fire count across all the iterations */
       prescribedfire,
      /** \brief print individual information like number of establishments by year. */
       indv;
      /** \brief the character used to separate values in the output files. */
  char sep;
};

/**
 * \brief Defines which mortality output files we would like, and what they should contain.
 * 
 * These flags, all from inputs, determine how many files should be printed as well as what 
 * the columns should contain. This struct is instanciated by the global variable MortFlags.
 * 
 * \sa MortFlags
 * 
 * \ingroup MORTALITY
 */
struct mortflags_st {
      /** \brief If FALSE print no mortality output. */
  Bool summary,
      /** \brief Print individual yearly data as well as a summary. */
       yearly,
      /** \brief Print a header line of names in each file. */
       header,
      /** \brief Print mortality data for each group. */
       group,
      /** \brief Print mortality data for each species. */
       species;
      /** \brief The separator between values in the output files. */
  char sep;
};

struct superglobals_st {
    size_t max_rgroups, /* Maximum resource groups allowed. */
           max_groupnamelen, /* Maximum resource group name length. */
           max_spp_per_grp, /* Maximum species allowed per resource group. */
           max_indivs_per_spp, /* Maximum individuals allowed per species. */
           max_speciesnamelen; /* Maximum species name length. */
    
    IntUS runModelIterations,
          runSpinupYears,
          runModelYears,
          nCells;		/* number of cells to use in Grid, only applicable if grid function is being used */

    IntL randseed;
};

#endif
