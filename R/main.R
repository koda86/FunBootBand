# ******************************************************************************
# *****************    Main script prediction bands    *************************
# ******************************************************************************
#
# Author: Daniel Koska
# Date: November 2022
# R version: 4.0.5
# Platform: x86_64-pc-linux-gnu
# ------------------------------------------------------------------------------
# This script accompanies the publication:

# Please cite: ...
#
#
# Values of 4 variables should be set/initialized before running the script:
# 1. The data and script directories ('dir.data', 'dir.script')
# 2. The example data set (in example_data())
# 3. The number of bootstrap iterations (n.boot)
# 4. The number of repeated calculations for the uncertainty estimation (n.rep)
#
#
# IMPORTANT: Currently, the script is designed for balanced data sets.
# Unbalanced designs (unequal number of observations) may lead to errors!
# ******************************************************************************

rm(list = ls())

# ------------------------------------------------------------------------------
# Set up the R environment (working, directory, packages, scripts)
# ------------------------------------------------------------------------------
# Please specify the correct paths
dir.script <- "~/floa/R"        # All R scripts of the package are stored here
dir.data <- "~/floa/R/examples" # Directory in which the data are stored

setwd(dir.script)

source("pick_curves.R")
source("floa_boot.R")
source("floa_boot_rep.R")
source("floa_point.R")
source("floa_roislien.R")
source("plot_loa.R")
source("points_within_limits.R")
source("coverage_loocv.R")
source("coverage_curves.R")
source("estimate_uncertainty_kfold_rep.R")

# Load installed packages and install missing packages (automatically)
# Packages needed for floa
packages = c("tidyverse",
             "reshape2",
             "matlab")

package.setup <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# ------------------------------------------------------------------------------
# Data sets
# ------------------------------------------------------------------------------
# Specify the name of your data object (data sets in the paper were created in
# 'simulate_data.R'). Own data sets need to have the same structure as those
# (long format data.frame with 6 columns: device, subjectID, strideID, value, frame).
#
# - Smooth wave data (normal error):       "smooth_realistic"
# - Data with non-gaussian error:          "non_gaussian"
# - Phase shifted data (x-axis direction): "shift"
# - Real world empirical validation data:  "imu_mc"
data <- readRDS(paste0(dir.data, "/", "smooth_realistic.rds"))


# ------------------------------------------------------------------------------
# Plot prediction bands from four methods
# ------------------------------------------------------------------------------
n.boot <- 1000

floa.point     <- floa_point(data)
floa.roislien  <- floa_roislien(data)
floa.boot.rep  <- floa_boot_rep(data,
                            k.coef = 50,
                            n.boot = n.boot,
                            band = "prediction",
                            cp.begin = 0,
                            alpha = 0.05)
floa.boot.iid  <- floa_boot(data,
                            k.coef = 50,
                            n.boot = n.boot,
                            band = "prediction",
                            cp.begin = 0,
                            alpha = 0.05,
                            iid = TRUE) # Draw only a single curve per subject

# Plot prediction bands and difference curves
# POINT (dotted), RØISLIEN (pink), BOOTrep (blue), BOOTiid (yellow)
plot_loa(data,
         floa.point, floa.roislien, floa.boot.rep, floa.boot.iid,
         ylim = c(-5, 5)) # values data set imu_mc: c(-20, 20)


# ------------------------------------------------------------------------------
# Leave-one out cross validation to estimate coverage properties
# ------------------------------------------------------------------------------
# Percentage of covered points
# 100% means that curve is fully contained within the band limits
# Order of output: POINT, ROISLIEN, BOOTrep, BOOTiid
n.strides <- length(unique(data$strideID))
cover.cross.fraction.singlecurve <- coverage_loocv(data, n.boot)

# Proportion of bands that contain p.cover percent of the T=101 curve points
apply(cover.cross.fraction.singlecurve, 2, coverage_curves, p.cover = 1)
apply(cover.cross.fraction.singlecurve, 2, coverage_curves, p.cover = 0.95)
apply(cover.cross.fraction.singlecurve, 2, coverage_curves, p.cover = 0.9)
apply(cover.cross.fraction.singlecurve, 2, coverage_curves, p.cover = 0.5)


# ------------------------------------------------------------------------------
# Uncertainty estimation using repeated k-Fold cross validation (k subjects)
# ------------------------------------------------------------------------------
n.rep <- 300

# Order of output: POINT, ROISLIEN, BOOT.all, BOOT.iid
uncertainty.area.kfold <- estimate_uncertainty_kfold_rep(data,
                                                         n.rep,
                                                         n.boot,
                                                         plot.au = TRUE,
                                                         plot.ylimits = c(-6, 6) # c(-12, 12)
)


