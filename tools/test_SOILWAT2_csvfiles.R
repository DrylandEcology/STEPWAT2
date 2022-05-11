#!/usr/bin/env Rscript

#-------------------------------------------
# R script to check that the `SOILWAT2`-based output for each `STEPWAT2`
# iteration agrees with the across-iteration aggregated output of means and
# standard deviations

# Consider turning on csv-output files for each of the four output time periods:
# --> set `TIMESTEP dy wk mo yr` in file `outsetup.in`

# Run this R code from within the `STEPWAT2/` directory after a
# `STEPWAT2` simulation has been run in the testing directory with
# ```
#   cd testing.sagebrush.master/Stepwat_Inputs/
#   ./stepwat -f files.in -s -i -o
#   cd ../..
#   Rscript tools/test_SOILWAT2_csvfiles.R
# ```

# Adjust variable `reps` below to reflect number of iterations from `STEPWAT2`
# run (as specified as `niter` in file `model.in`).
#-------------------------------------------

dir <- file.path(
  "testing.sagebrush.master",
  "Stepwat_Inputs",
  "Output",
  "sw_output"
)

report_by_var <- FALSE
report_by_period <- TRUE

reps <- 1
OutPeriods <- c("daily", "weekly", "monthly", "yearly")
agg_type <- c("", "_slyrs")

for (atype in agg_type) {
  x <- read.csv(file.path(dir, paste0("sw2_yearly", atype, "_rep1.csv")))

  vars <- colnames(x)[-c(1, ncol(x))]

  for (k in OutPeriods) {
    dev_sum <- c(mean = 0, SD = 0)

    xd <- lapply(seq_len(reps), function(i)
      try(read.csv(file.path(dir, paste0("sw2_", k, atype, "_rep", i, ".csv")))))

    xdm <- try(read.csv(file.path(dir, paste0("sw2_", k, atype, "_agg.csv"))))

    has_failed <- any(
      sapply(xd, function(x) inherits(x, "try-error")),
      inherits(xdm, "try-error"))

    if (has_failed) {
      next
    }

    for (v in vars) {
      reagg <- matrix(NA, nrow = nrow(xdm), ncol = reps)
      for (i in seq_len(reps)) {
        reagg[, i] <- xd[[i]][, v]
      }

      reagg_mean <- rowMeans(reagg)
      dmean <- sum(abs(xdm[, paste0(v, "_Mean")] - reagg_mean))
      dev_sum["mean"] <- dev_sum["mean"] + dmean

      reagg_sd <- apply(reagg, 1, sd)
      dsd <- sum(abs(xdm[, paste0(v, "_SD")] - reagg_sd))
      dev_sum["SD"] <- dev_sum["SD"] + dsd

      if (report_by_var) {
        print(paste0(atype, k, ": ", shQuote(v),
          " has sum of mean-deviation = ", round(dmean, 3)))
        print(paste0(atype, k, ": ", shQuote(v),
          " has sum of sd-deviation = ", round(dsd, 3)))
      }
    }

    if (report_by_period) {
      print(paste0(atype, k, ": has total sum of mean-deviation = ",
        round(dev_sum["mean"], 3)))
      print(paste0(atype, k, ": has total sum of sd-deviation = ",
        round(dev_sum["mean"], 3)))
    }
  }
}
