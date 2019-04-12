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

/* structure for indiv of perennial species */
struct indiv_st {
  IntUS id,
  	  age,
      mm_extra_res,
      slow_yrs,  /* number of years this individual has slow growth */
      myspecies,
      yrs_neg_pr;     /* num yrs pr > 1 (stretched resources) */
  MortalityType killedby;
  Bool killed;        /* only for clonal plants, means veggrow possible*/
  RealF relsize,      /* relative to full-sized individual -- 0-1.0) */
       prv_yr_relsize, /*previous year relsize relative to full-sized individual, saving before killing, used for proportional recovery calculation) */
       grp_res_prop,  /* prop'l contribution this indiv makes to group relsize */
       res_required, /* resources required, biomass of the individual */
       res_avail,    /* resource available */
       res_extra,    /* resource applied to superficial growth */
       pr,           /* ratio of resources required to amt available */
       growthrate,   /* actual growth rate*/
       normal_growth, /* biomass the plant gained this year excluding superfluous biomass. Used for grazing.*/
       prob_veggrow; /* set when killed; 0 if not clonal*/
  struct indiv_st *Next, *Prev;  /* facility for linked list 8/3/01 */
};

/* structure for indiv of annual species */
struct indiv_ann_st {
	 // (DLM) - 6/5/2013 - this struct was defined when I got steppe, but now that I've been analyzing the individuals file (and it's linked list) I don't actually see it being used anywhere... it seems like it's just a relic of some older version of steppe that was never removed.  Will leave it in the code for now though in case someone else knows.
  IntUS mm_extra_res,
      myspecies;
  RealF relsize,      /* relative to full-sized individual -- 0-1.0) */
       grp_res_prop,  /* prop'l contribution this indiv makes to group relsize */
       res_required, /* resources required, biomass of the individual */
       res_avail,    /* resource available */
       res_extra,    /* resource applied to superficial growth */
       pr,           /* ratio of resources required to amt available */
       growthrate;   /* actual growth rate*/
  struct indiv_ann_st *Next, *Prev;  /* facility for linked list 8/3/01 */
};

/* structure for species with some annuals-only variables thrown in */
struct species_st {

  /**** Quantities that can change during model runs *****/

  SppIndex est_count;  /* number of individuals established (growing)*/
  IntUS *kills,        /* ptr to array of # indivs killed by age. index=age. */
        estabs,        /* number of individuals established in iter */
        *seedprod,    /* annuals: array of previous years' seed production (size = viable_yrs)*/
         seedbank,
         pseed;
  RealF lastyear_relsize,    /* relsize from the previous year, used for annual establishment */
        extragrowth,   /* amt of superfluous growth from extra resources */
	received_prob,	//the chance that this species received seeds this year... only applicable if using seed dispersal and gridded option
        alpha,         /* species-specific alpha parameter used for random number draw from beta distribution in establishment of annual species */
        beta;           /* species-specific beta parameter used for random number draw from beta distribution in establishment of annual species */
  float       var;    /* the variance parameter of the beta distribution for establishment..*/
  struct indiv_st *IndvHead;    /* facility for linked list 8/3/01; top of list */
  Bool allow_growth, //whether to allow growth this year... only applicable if using seed dispersal and gridded option
	sd_sgerm;	//whether seeds where produced/received and germinated this year... only applicable if using seed dispersal and gridded option

  /**** Quantities that DO NOT change during model runs *****/

    /* 4-letter code (plus \0) for genus and species*/
  char name[MAX_SPECIESNAMELEN+1];
  IntUS max_age,         /* max age of mature plant, also flag for annual */
        viable_yrs,      /* annuals: max years of viability of seeds */
        max_seed_estab,  /* max seedlings that can estab in 1 yr*/
        max_vegunits,    /* max vegetative regrowth units (eg tillers)*/
        max_slow,        /* years slow growth allowed before mortality */
        sp_num,          /* index number of this species */
        res_grp;         /* this sp. belongs to this res_grp*/
  RealF max_rate,      /* intrin_rate * proportion*/
        intrin_rate,
        relseedlingsize,
        seedling_biomass,
        mature_biomass,    /* biomass of mature_size individual */
        seedling_estab_prob_old,  /* supports Extirpate() and Kill() */
        seedling_estab_prob,
        ann_mort_prob,
        cohort_surv,
        exp_decay,        /* annuals: exponent for viability decay function */
        prob_veggrow[4],  /* 1 value for each mortality type, if clonal*/
  	sd_Param1,	  /* for seed dispersal */
  	sd_PPTdry,
	sd_PPTwet,
	sd_Pmin,
	sd_Pmax,
	sd_H,
	sd_VT;
  TempClass tempclass;
  DisturbClass disturbclass;
  Bool isclonal,
       use_temp_response,
       use_me,           /* do not establish if this is false */
       use_dispersal;	//whether to use seed dispersal... only applicable if using gridded option
};

struct resourcegroup_st {

  /**** Quantities that can change during model runs *****/

  IntUS *kills,        /* indivs in group killed. index by age killed. */
        estabs,         /* total indivs in group established during iter */
        killyr,         /* kill the group in this year; if 0, don't kill, but see killfreq */
        yrs_neg_pr,     /* counter for consecutive years low resources */
        wildfire, 	/* number of wildfires during all iterations in this year */
        prescribedfire, /* number of prescribed fires during all iterations in this year */
        mm_extra_res;   /* extra resource converted back to mm */
  RealF res_required, /* resource required for current size */
        res_avail,    /* resource available from environment X competition */
        res_extra,    /* resource applied to superficial growth */
        pr,           /* resources required / resources available */
        rgroupFractionOfVegTypeBiomass; /*proportional biomass of the STEPPE functional group out of the SOILWAT2 functional type biomass */
  SppIndex est_count, /* number of species actually established in group*/
           est_spp[MAX_SPP_PER_GRP]; /*list of spp actually estab in grp*/
  Bool extirpated,    /* group extirpated, no more regen */
       regen_ok;      /* annuals: startyr; TM this is a flag for annuals.  When
                       * you set the start year for the group this will flag all
                       * species within that group that establishment will start.
                       * EG start year = 96 I will have 5 flags for one annual species
                       * for year 96,97,98,99,100. BUT extirpate will not un-flag this*/

  /**** Quantities that DO NOT change during model runs *****/

  IntUS max_stretch,    /* num yrs resources can be stretched w/o killing*/
        max_spp_estab,  /* max # species that can add new plants per year*/
        max_spp,        /* number of species in the group*/
        max_age,        /* longest lifespan in group. used to malloc kills[] */
        startyr,        /* don't start trying to grow until this year */
	killfreq_startyr,/* start year for kill frequency*/

        extirp,         /* year in which group is extirpated (0==ignore) */
        grp_num,        /* index number of this group */
        veg_prod_type,  /* type of VegProd.  1 for tree, 2 for shrub, 3 for grass, 4 for forb */
		grazingfrq,     /* grazing effect on group at this frequency: <1=prob, >1=# years */
        grazingfreq_startyr;/* start year for grazing frequency*/
  SppIndex species[MAX_SPP_PER_GRP]; /*list of spp belonging to this grp*/
  RealF min_res_req,  /* input from table */
        max_density,  /* number of mature plants per plot allowed */
        max_per_sqm,  /* convert density and plotsize to max plants/m^2 */
        max_bmass,    /* sum of mature biomass for all species in group */
        killfreq,       /* kill group at this frequency: <1=prob, >1=# years */
        ignition,       /* cheatgrass biomass (g/m2) that triggers potential ignition of a wildfire */
        cheatgrass_coefficient,   /* intercept of the cheatgrass biomass-wildfire probability relationship */
        wild_fire_slope,  /* slope of the cheatgrass biomass-wildfire probability relationship */	
        xgrow,        /* ephemeral growth = mm extra ppt * xgrow */
        slowrate,     /* user-defined growthrate that triggers mortality */
        ppt_slope[3], /* res. space eqn: slope for wet/dry/norm yrs*/
        ppt_intcpt[3],/* res. space eqn: intercept for "" ""*/
		proportion_killed,      /* proportion killing  */
		proportion_recovered,   /* proportion recovery after killing year */
        proportion_grazing;     /* proportion grazing on grazing year */
  Bool succulent,
       use_extra_res, /* responds to other groups' unused resources */
       use_me,        /* establish no species of this group if false; TMartyn5.26.15 - this
                       * is 0 of the onoff in species.in and of the on in rgroup.in */
       use_mort,      /* use age-independent+slowgrowth mortality?  */
       est_annually;  /* establish this group every year if true  */
  DepthClass depth;  /* rooting depth class */
  char name[MAX_GROUPNAMELEN +1];
};

struct succulent_st {
  RealF growth[2], /* growth modifier eqn parms for succulents (eqn 10)*/
        mort[2],   /* mortality eqn parms for succulents (eqn 16)*/
        reduction, /* if not killed, reduce by eqn 10*/
        prob_death;  /* calculated from eqn 16*/
};

struct environs_st {
  PPTClass wet_dry;
  IntS ppt,  /* precip for the year (mm)*/
       lyppt, /* precip for the previous (last) year (mm)*/
       gsppt;  /* precip during growing season (mm)*/
  RealF temp,  /* average daily temp for the year (C)*/
       temp_reduction[2]; /* amt to reduce growth by temp */
                           /*(eqns 12,13), one for each tempclass*/
};

struct plot_st {  /* plot-level things */
  DisturbEvent disturbance;
  Bool pat_removed;  /* fecalpats can be removed which only*/
                     /* kills seedlings and not other plants*/
  IntUS disturbed;  /* years remaining before recolonization*/
                  /* (ie, new establishments) can begin again,*/
                  /* or, if disturbance is fecalpat, number of*/
                  /* years it has been ongoing.  Set to 0*/
                  /* when the disturbance effect is expired.*/
};


struct ppt_st {
  RealF avg,
       std;
  IntUS  min,
       max,
       dry,
       wet;
};
struct temp_st {
  RealF avg,
        std,
        min,
        max;
};
struct fecalpats_st {
  Bool use;
  RealF occur,
        removal,
        recol[2]; /* 1 elem. for slope and intercept*/
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
  char header_line[1024]; /* contains output header line, used to print */
                          /* yearly values and statistics  */

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
