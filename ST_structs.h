/********************************************************/
/********************************************************/
/*  Source file: ST_structs.h
 *  Type: header
 *  Application: STEPPE - plant community dynamics simulator
 *  Purpose: This is the most interesting header file where
 *           all of the "objects" are defined as structures.
 *           You can think of this as the object dictionary
 *           insofar as you would want to have this file
 *           handy to refer to when perusing the code.
 *  History
 *  History:
 *     6/15/2000 -- INITIAL CODING - cwb
 *     4-Nov-03 (cwb) Added code to handle annuals.
 *        Annuals are different because they complete all the
 *        dynamics of their life cycle in a single time step. Thus the
 *        mechanism for establishing and growing individuals and
 *        associated variables must be different from perennials. For
 *        example, there are no clonal annuals, no multi-year variables
 *        (eg max_slow), annualized probability of mortality, etc.
 *        Basically, the number of individuals is determined by the
 *        available resources and a seedbank mechanism, and the growth
 *        rate is simply the proportion of max size (1.0) each indiv
 *        can achieve in the current year, computed by maxbio * irate /
 *        PR. And of course they all die at the end of the time step.
 *
 *        The method of establishing annuals is based on the number of
 *        viable seeds produced in the previous year and in all years
 *        previous up to some limit of age-related viability (eg 10
 *        years).  Production of viable seeds is related to this year's
 *        resource availability for the group (recall PR is inverse of
 *        resource availability). Viability decreases with age
 *        (=1/age), so an array is kept with the last X year's seed
 *        production.  Maximum possible establishment is the sum of the
          past years' production, weighted by 1/(seed_age^xdecay). */
/********************************************************/
/********************************************************/


#ifndef STEPPE_STRUCT_DEF
#define STEPPE_STRUCT_DEF

#include "generic.h"

/**
 * \brief Holds information on perenial plant individuals. 
 * 
 * An individual is one plant specimen. Most of the code associated with this
 * struct is found in ST_indivs.c.
 * 
 * \sa indiv_ann_st
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
       /** \brief Resources applied to superficial growth. */
       res_extra,
       /** \brief Ratio of resources required to amt available. */
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
       /** \brief Allows for a doubly-linked list of individuals. Implemented in Species.
        * \sa Species */
  struct indiv_st *Next, 
       /** \brief Allows for a doubly-linked list of individuals. Implemented in Species.
        * \sa Species */
                  *Prev;
};

/** 
 * \brief Holds information on annual plant individuals.
 * 
 * An individual is one plant specimen. Most of the code associated with this
 * struct is found in ST_indivs.c.
 * 
 * \sa indiv_st
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
       /** \brief Resources applied to superficial growth. */
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
 * species. Species is a global array of these structs.
 * 
 * \sa Species
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
       * \sa ST_seedDispersal.c */
        *seedprod,
      /** \brief not sure what this does yet.
       * \sa ST_seedDispersal.c */
         seedbank,
      /** \brief not sure what this does yet.
       * \sa ST_seedDispersal.c */
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
      /** \brief Variance parameter of the beta distribution for establishment. */
  float var;
      /** \brief Head of a doubly-linked list of all individuals of this species. 
       * \sa indiv_st*/
  struct indiv_st *IndvHead;
      /** \brief Seed dispersal only- whether to allow growth for the current year.
       * \sa ST_seedDispersal.c */
  Bool allow_growth,
      /** \brief Whether seeds where produced/received and germinated this year. 
       * \sa ST_seedDispersal.c */
	     sd_sgerm;

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
      /** \brief Starting biomass in grams of a new individual of this species. */
        seedling_biomass,
      /** \brief Biomass in grams of a mature individual. */
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
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
  	    sd_Param1,
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
  	    sd_PPTdry,
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
	      sd_PPTwet,
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
	      sd_Pmin,
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
  	    sd_Pmax,
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
	      sd_H,
      /** \brief Seed dispersal parameter read from inputs.
       * \sa ST_seedDispersal.c */
      	sd_VT;
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
       * \sa ST_seedDispersal.c */
       use_dispersal;
};

/** 
 * \brief Contains all resource group-specific fields. 
 * 
 * Some fields are consant and some are variable, so be careful what you are modifying. 
 * RGroup is a global array of these structs used in STEPWAT2.
 * 
 * \sa RGroup 
 */
struct resourcegroup_st {

  /**** Quantities that can change during model runs *****/

      /** \brief individuals in group killed. Index by age at death. */
  IntUS *kills,
      /** \brief Total indivs in group established during current iteration. */
        estabs,
      /** \brief Treated as a Bool. If TRUE kill the group this year.
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
      /** \brief Xxtra resources in millimeters */
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
        rgroupFractionOfVegTypeBiomass;
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

      /** \brief Number of yers resources can be stretched without killing the group. */
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
      /** \brief Perform grazing on this group at this frequency. */
		    grazingfrq,
      /** \brief Start year for grazing frequency.
       * \sa grazingfrq */
        grazingfreq_startyr;
      /** \brief Array of species indexes belonging to this group.
       * \sa Species */
  SppIndex *species;
      /** \brief Minimum resources required by this group, defined in inputs. */
  RealF min_res_req,
      /** \brief Number of mature plants allowed per plot. */
        max_density,
      /** \brief Max plants per meter squared.
       * \sa max_density */
        max_per_sqm,
      /** \brief Max biomass of group. */
        max_bmass,
      /** \brief Kill group at this frequency.
       * Values < 1 result in probablistic fires
       * values > 1 result in fires ever killfreq years. 
       * \sa mort_EndOfYear() */
        killfreq,
      /** \brief Cheatgrass biomass (g/m2) that triggers potential ignition of a wildfire.
       * If ignition == 0 wildfire is turned off. 
       * \sa mort_EndOfYear() */
        ignition,
      /** \brief Intercept of the cheatgrass biomass-wildfire probability relationship. */
        cheatgrass_coefficient,
      /** \brief Slope of the cheatgrass biomass-wildfire probability relationship. */	
        wild_fire_slope,
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

struct antmounds_st {
  Bool use;
  RealF occur;
  IntUS minyr,
       maxyr;
};
struct burrows_st {
  Bool use;
  RealF occur;
  IntUS minyr;
};

struct outfiles_st {
  FILE *fp_year,  /* file handle for yearly so it can stay open*/
       *fp_sumry; /* file handle for averages output */
  IntUS suffixwidth; /* max width of outfile suffix if printing yearly */
};

struct globals_st {
  struct ppt_st ppt;
  struct temp_st temp;
  struct fecalpats_st pat;
  struct antmounds_st mound;
  struct burrows_st burrow;

  RealF plotsize,   /* size of plot in square meters */
        gsppt_prop, /* proportion of ppt during growing season*/
        tempparm[2][3]; /* three parms for Warm/Cool growth mod*/
  IntUS runModelYears,  /* number of years to run the model*/
      Max_Age,        /* oldest plant; same as runModelYears for now */
      currYear,
      runModelIterations, /* run model this many times for statistics */
      currIter,
      grpCount,     /* number of groups defined*/
      sppCount,     /* number of species defined*/
      transp_window, /* Number of years for which transpiration data is kept*/
      nCells;		/* number of cells to use in Grid, only applicable if grid function is being used */
  IntL randseed;     /* random seed from input file */
  size_t max_rgroups, /* Maximum resource groups allowed. */
         max_groupnamelen, /* Maximum resource group name length. */
         max_spp_per_grp, /* Maximum species allowed per resource group. */
         max_indivs_per_spp, /* Maximum individuals allowed per species. */
         max_speciesnamelen; /* Maximum species name length. */

  struct outfiles_st bmass, mort;
};

struct bmassflags_st {
  Bool summary,  /* if FALSE, print no biomass output */
       yearly, /* print individual yearly runs as well as average */
       header,
       yr,
       dist,
       ppt,
       pclass,
       tmp,
       grpb,
       pr,
       size,
       sppb,
       wildfire,/* print wild fires count during all the iterations */
       prescribedfire,/* print prescribed fires count during all the iterations */
       indv;
  char sep;
};

struct mortflags_st {
  Bool summary,  /* if FALSE, print no mortality output */
       yearly, /* print individual yearly data as well as summary */
       header, /* print a header line of names in each file */
       group,  /* print data summarized by group */
       species; /* print data for species */
  char sep;
};




#endif
