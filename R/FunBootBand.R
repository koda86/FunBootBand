#' @title FunBootBand
#'
#' @description Creates Functional Bootstraped (statistical) Bands
#'
#' @param data Data set
#' @param type Band type (type = c("confidence", "prediction", "tolerance"))
#' @param B Number of bootstrap iterations (e.g., B = 1000)
#' @param iid Assume independent and identically distributed (iid) curves or not
#' (iid = c(TRUE, FALSE))
#' @param alpha Type I error probability
#'
#' @return A data frame object that contains upper and lower band boundaries
#' @examples
#' load("curvesample.RData")
#' band.limits <- band(data = curves, type = "prediction", B = 1000, iid = TRUE)
#' @export
#' @import tidyverse, reshape2, matlab

# ******************************************************************************
# IMPORTANT: Currently, the script is designed for balanced data sets.
# Unbalanced designs (unequal number of observations) may lead to errors!
# ******************************************************************************

# # Load installed packages and install missing packages (automatically)
# # Packages needed for floa
# packages = c("tidyverse",
#              "reshape2",
#              "matlab")
#
# package.setup <- lapply(
#   packages,
#   FUN = function(x) {
#     if (!require(x, character.only = TRUE)) {
#       install.packages(x, dependencies = TRUE)
#       library(x, character.only = TRUE)
#     }
#   }
# )

data <- get(load("curvesample.RData"))

# floa.boot.rep  <- floa_boot_rep(data,
#                                 k.coef = 50,
#                                 B = B,
#                                 band = "prediction",
#                                 cp.begin = 0,
#                                 alpha = 0.05)
# floa.boot.iid  <- floa_boot(data,
#                             k.coef = 50,
#                             B = B,
#                             band = "prediction",
#                             cp.begin = 0,
#                             alpha = 0.05,
#                             iid = TRUE) # Draw only a single curve per subject

# Actual function
band <- function(data, k.coef = 50, B = 50, type = "prediction", cp.begin = 0,
                 alpha = 0.05, iid = TRUE) {
  # ****************************************************************************
  # In this script, functional prediction bands are calculated by adapting the
  # method described in Lenhoff et al. (1999).
  # ****************************************************************************
  # Dimensions
  n.time    <- dim(data)[1] # unique(data$frame)
  n.curves  <- dim(data)[2] # unique(data$strideID)
  time      <- seq(0, (n.time-1))

  # ----------------------------------------------------------------------------
  # Approximate time series (differences) using Fourier functions
  # ----------------------------------------------------------------------------
  fourier.koeffi    <- matlab::zeros(c(k.coef*2 + 1, n.curves))
  fourier.real      <- matlab::zeros(n.time, n.curves)
  fourier.mean      <- matlab::zeros(k.coef*2 + 1)
  fourier.real_mw   <- matlab::zeros(n.time, 1)
  fourier.std1      <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, n.curves)
  fourier.kovarianz <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1)
  fourier.std_all   <- matlab::zeros(n.time, n.time)
  fourier.std       <- matlab::zeros(n.time, 1)

  # Set up a Fourier series
  # General: f(t) = mu + sum(alpha cos(2pi*k*t/T) + beta sin(2pi*k*t/T))
  fourier.s = rep(1, times = n.time)
  for (k in seq(1, k.coef*2, 2)) {
    fourier.s <- cbind(fourier.s, cos(2*pi*(k/2)*time/(n.time-1)))
    fourier.s <- cbind(fourier.s, sin(2*pi*(k/2)*time/(n.time-1)))
  }

  for (i in 1:n.curves) {
    # Least squares Regression
    fourier.koeffi[, i] = pracma::mldivide(fourier.s, data[, i])
    # Fourier curve
    fourier.real[, i] = fourier.s %*% fourier.koeffi[, i]
  }

  # Mean Fourier curve
  fourier.mean[, 1] = rowMeans(fourier.koeffi)
  fourier.real_mw[, 1] = fourier.s %*% fourier.mean[, 1]

  # Standard deviation of the Fourier curve
  for (i in 1:n.curves) {
    # variance-covariance matrix
    fourier.std1[, , i] <- (fourier.koeffi[, i] - fourier.mean[, 1]) %*% t(fourier.koeffi[, i] - fourier.mean[, 1])
  }

  fourier.kovarianz <- apply(fourier.std1, c(1, 2), mean)
  # Lenhoff, Appendix A, Eq. (0.5)
  fourier.std_all <- suppressWarnings(sqrt(fourier.s %*% fourier.kovarianz %*% t(fourier.s)))

  for (i in 1:n.time) {
    # Values are on the diagonal of the square matrix fourier.std_all
    fourier.std[i, 1] = fourier.std_all[i, i]
  }


  # ----------------------------------------------------------------------------
  # Bootstrap
  # ----------------------------------------------------------------------------
  bootstrap_sample        <- matlab::zeros(n.time, 4)
  bootstrap.mean          <- matlab::zeros(k.coef*2 + 1, B)
  bootstrap.real_mw       <- matlab::zeros(n.time, B)
  bootstrap.zz            <- matlab::zeros(n.curves, B)
  bootstrap.pseudo_koeffi <- matlab::zeros(k.coef*2 + 1, n.curves, B)
  bootstrap.real          <- matlab::zeros(n.time, n.curves, B)
  bootstrap.std1          <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, n.curves)
  bootstrap.kovarianz     <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, B)
  bootstrap.std_all       <- matlab::zeros(n.time, n.time, B)
  bootstrap.std           <- matlab::zeros(n.time, B)

  for (i in 1:B) {
    for (k in 1:n.curves) {
      bootstrap.zz[k, i] = sample(n.curves, size=1)
      bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
      bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
    }

    # Mean bootstrap curve and standard deviation
    bootstrap.mean[, i] <- rowMeans(bootstrap.pseudo_koeffi[, , 1])
    bootstrap.real_mw[, i] <- fourier.s %*% bootstrap.mean[, i]

    for (k in 1:n.curves) {
      bootstrap.std1[, , k] <- (bootstrap.pseudo_koeffi[, k, i] - bootstrap.mean[, i]) %*% t(bootstrap.pseudo_koeffi[, k, i] - bootstrap.mean[, i])
    }

    bootstrap.kovarianz[, , i] <- apply(bootstrap.std1, c(1, 2), mean)
    bootstrap.std_all[, , i] <- suppressWarnings(sqrt(fourier.s %*% bootstrap.kovarianz[, , i] %*% t(fourier.s)))

    for (k in 1:n.time) {
      bootstrap.std[k, i] <- bootstrap.std_all[k, k, i]
    }

    # Dummy plots (Appendix) (which curves are drawn per bootstrap iteration?)
    # condition <- ifelse(iid == TRUE, "iid", "rep")
    # filename.dummy <- paste0("~/Nextcloud/project-fab-forschung/Publikationen/FLOA/paper_JournalBiomechanics/Review/dummy_plots/",
    #                          condition,
    #                          i, # bootstrap iteration
    #                          ".png")
    # png(filename.dummy)
    # plot(bootstrap.real[, 1, 1],
    #      type = "l",
    #      ylim = c(-2, 2),
    #      main = paste0("BOOTiid (iteration ", i, "), SD at 87%: ", round(bootstrap.std[, i][87], 2)),
    #      xlab = "Time-normalized signal [%]",
    #      ylab = "Difference")
    # for (ii in 1:dim(bootstrap.real)[2]) {
    #   lines(bootstrap.real[, ii, i])
    # }
    # dev.off()
  }



  # ----------------------------------------------------------------------------
  # Construct prediction or confidence bands
  # ----------------------------------------------------------------------------
  floa.boot.mean <- rowMeans(bootstrap.real_mw)
  floa.boot.sd   <- rowMeans(bootstrap.std)

  if (type == "prediction") {
    cp.data   <- matlab::zeros(n.curves, B)
    cp.data_i <- matlab::zeros(n.curves, B)

    cp.mean <- 0
    cp.bound <- cp.begin
    while (cp.mean < (1-alpha)) {
      for (i in 1:B) {
        for (k in 1:n.curves) {
          # Lenhoff et al., Appendix A, Eq. (0.6)
          cp.data[k, i] <- max(abs(fourier.real[, k] - bootstrap.real_mw[, i]) / bootstrap.std[, i])
          cp.data_i[k, i] <- cp.data[k, i] < cp.bound
        }
      }
      cp.mean <- mean(cp.data_i)
      cp.bound <- cp.bound + 0.05
    }
    cp_out <- cp.bound

    # Construct bands
    # ------------------------------------------------------------------------
    floa.boot <- rbind(floa.boot.mean + cp.bound * floa.boot.sd,
                       floa.boot.mean,
                       floa.boot.mean - cp.bound * floa.boot.sd
    )
  } else { # confidence bands
    # cc.data <- matlab::zeros(n.curves, B)

    for (i in 1:B) {
      for (k in 1:n.curves) {
        # Lenhoff, Appendix A, Eq. (0.6)
        cc.data[k, i] <- max(abs(fourier.real_mw[, k] - bootstrap.real_mw[, i]) / bootstrap.std[, i])
      }
    }
    cc <- quantile(cc.data, probs = 1-alpha)

    floa.boot <- rbind(floa.boot.mean + cc * floa.boot.sd,
                       floa.boot.mean,
                       floa.boot.mean - cc * floa.boot.sd
    )
  }

  row.names(floa.boot) <- c("upper.loa", "mean", "lower.loa")

  return(floa.boot)
}

