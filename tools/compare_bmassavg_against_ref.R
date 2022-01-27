#!/usr/bin/env Rscript

#-------------------------------------------
# R script to compare biomass output of a non-gridded testing run
# (e.g., on a development branch)
# against a previously executed "reference" run (e.g., on the release branch)


# Run this R code from within the directory of `STEPWAT2/` after
# two non-gridded `STEPWAT2` simulations have been run in the testing directory
#   cd testing.sagebrush.master/Stepwat_Inputs/
#   ./stepwat -f files.in -o
#   cp -r Output/ Output_ref/
#   ./stepwat -f files.in -o

#-------------------------------------------

#--- Define inputs ------
dir <- "testing.sagebrush.master/Stepwat_Inputs/Output"
dir_ref <- "testing.sagebrush.master/Stepwat_Inputs/Output_ref"
fname_bmass <- "bmassavg.csv"

nyrs_meanbmass <- 200 # Used by Palmquist et al. 2021 GCB

# rgroups from "Input/rgoup.in"
rgroups <- c(
  "sagebrush",
  "a.cool.forb",
  "a.warm.forb",
  "p.cool.forb",
  "p.warm.forb",
  "a.cool.grass",
  "p.cool.grass",
  "p.warm.grass",
  "shrub",
  "succulents"
)


#--- Load values ------
x <- read.csv(file.path(dir, fname_bmass))
x_ref <- read.csv(file.path(dir_ref, fname_bmass))

ids_yrs <- max(1, nrow(x) - nyrs_meanbmass + 1):nrow(x)

ids <-
  apply(x[ids_yrs, rgroups], 2, max) > 0 |
  apply(x_ref[ids_yrs, rgroups], 2, max) > 0
rgroups_used <- rgroups[ids]


#--- Comparisons ------
nrgu <- length(rgroups_used)

colors_rg <- grDevices::hcl.colors(nrgu, palette = "viridis")
colors_rg <- grDevices::palette.colors(nrgu, palette = "Okabe-Ito")
lty_rg <- rep_len(1:5, nrgu)

meanx <- colMeans(x[ids_yrs, rgroups_used])
meanxref <- colMeans(x_ref[ids_yrs, rgroups_used])


#--- * Time-series comparisons ------
n_panels <- c(3, 2)

fname_tsc <- file.path(
  dirname(dir),
  paste0("Fig_TimeSeriesComparison_", format(Sys.time(), "%Y%m%d-%H%M"), ".pdf")
)

if (!file.exists(fname_tsc)) {
  pdf(
    file = fname_tsc,
    height = n_panels[1] * 3,
    width = n_panels[2] * 5
  )

  par_prev <- par(
    mfrow = n_panels,
    mar = c(2.5, 3, 3, 1),
    mgp = c(1.5, 0, 0),
    tcl = 0.3
  )

  matplot(
    x[, rgroups_used],
    type = "l",
    xlab = "Simulation time [Years]",
    ylab = "Biomass of test run [g/m2]",
    col = colors_rg,
    lty = lty_rg,
    main = "Biomass time-series"
  )
  points(rep(nrow(x) + 5, nrgu), meanx, col = colors_rg, pch = 4)
  legend(
    "topleft",
    legend = rgroups_used,
    col = colors_rg,
    lwd = 3,
    cex = 0.65
  )
  matplot(
    x_ref[, rgroups_used],
    type = "l",
    xlab = "Simulation time [Years]",
    ylab = "Biomass of reference run [g/m2]",
    col = colors_rg,
    lty = lty_rg
  )
  points(rep(nrow(x_ref) + 5, nrgu), meanxref, col = colors_rg, pch = 4)

  matplot(
    x[, rgroups_used],
    type = "l",
    log = "y",
    xlab = "Simulation time [Years]",
    ylab = "Biomass of test run [g/m2]",
    col = colors_rg,
    lty = lty_rg,
    main = "Biomass time-series on log-scale"
  )
  points(rep(nrow(x) + 5, nrgu), meanx, col = colors_rg, pch = 4)
  matplot(
    x_ref[, rgroups_used],
    type = "l",
    log = "y",
    xlab = "Simulation time [Years]",
    ylab = "Biomass of reference run [g/m2]",
    col = colors_rg,
    lty = lty_rg
  )
  points(rep(nrow(x_ref) + 5, nrgu), meanxref, col = colors_rg, pch = 4)

  matplot(
    x[, rgroups_used],
    type = "l",
    ylim = c(0, 200),
    xlab = "Simulation time [Years]",
    ylab = "Biomass of test run [g/m2]",
    col = colors_rg,
    lty = lty_rg,
    main = "Biomass time-series on trimmed scale"
  )
  points(rep(nrow(x) + 5, nrgu), meanx, col = colors_rg, pch = 4)
  matplot(
    x_ref[, rgroups_used],
    type = "l",
    ylim = c(0, 200),
    xlab = "Simulation time [Years]",
    ylab = "Biomass of reference run [g/m2]",
    col = colors_rg,
    lty = lty_rg
  )
  points(rep(nrow(x_ref) + 5, nrgu), meanxref, col = colors_rg, pch = 4)

  par(par_prev)
  dev.off()
}


#--- * 1:1 comparisons ------
n_panels <- c(nrgu, 2)

fname_scc <- file.path(
  dirname(dir),
  paste0("Fig_ScatterComparison_", format(Sys.time(), "%Y%m%d-%H%M"), ".pdf")
)

if (!file.exists(fname_scc)) {
  pdf(
    file = fname_scc,
    height = n_panels[1] * 3,
    width = n_panels[2] * 3
  )

  par_prev <- par(
    mfrow = n_panels,
    mar = c(3.5, 4.5, 3, 1),
    mgp = c(2, 0, 0),
    tcl = 0.3
  )

  for (k in seq_len(nrgu)) {
    tmpx <- x[ids_yrs, rgroups_used[k]]
    tmpxref <- x_ref[ids_yrs, rgroups_used[k]]

    vlim <- c(0, max(tmpx, tmpxref))
    plot(
      tmpxref,
      tmpx,
      col = colors_rg[k],
      xlim = vlim,
      ylim = vlim,
      xlab = paste(
        "\nBiomass of reference run [g/m2]\nmean =",
        round(meanxref[k], 2)
      ),
      ylab = paste(
        "Biomass of test run [g/m2]\nmean =",
        round(meanx[k], 2)
      ),
      main = rgroups_used[k]
    )
    abline(0, 1, col = "red", lty = 2)
    abline(v = meanxref[k], h = meanx[k], col = "orange", lwd = 2)

    diffs <- tmpxref - tmpx
    tmp <- range(diffs)
    vlim <- c(if (tmp[1] < 0) 1.05 else 0.95, 1.05) * tmp
    hist(
      diffs,
      col = colors_rg[k],
      xlim = vlim,
      xlab = paste(
        "\nBiomass of reference - test [g/m2]\nmean =",
        round(mean(diffs), 2)
      ),
      main = rgroups_used[k]
    )
    abline(v = 0, col = "red", lty = 2)
    abline(v = mean(diffs), col = "orange", lwd = 2)
  }

  par(par_prev)
  dev.off()
}
