#' @title FunBootBand
#'
#' @description Creates Functional Bootstraped (statistical) Bands
#'
#' @param data Data set
#' @param type Band type (type = c("confidence", "prediction", "tolerance"))
#' @param B Number of bootstrap iterations (e.g., B = 1000)
#' @param iid Assume independent and identically distributed (iid) curves or not
#' (iid = c(TRUE, FALSE))
#'
#' @return A data frame object that contains upper and lower band boundaries
#' @examples
#' data(imu_mc)
#' band.limits <- band(data = imu_mc, type = "prediction", B = 1000, iid = TRUE)
#' @export
#' @importFrom dplyr "%>%"

# R version: 4.0.5
# Platform: x86_64-pc-linux-gnu
#
# IMPORTANT: Currently, the script is designed for balanced data sets.
# Unbalanced designs (unequal number of observations) may lead to errors!
# ******************************************************************************

# ------------------------------------------------------------------------------
# Set up the R environment (working, directory, packages, scripts)
# ------------------------------------------------------------------------------
# Please specify the correct paths
# dir.script <- "~/floa/R"        # All R scripts of the package are stored here
# dir.data <- "~/floa/R/examples" # Directory in which the data are stored

# setwd(dir.script)

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

data <- load("curvesample.RData")

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

