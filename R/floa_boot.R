floa_boot <- function(data, k.coef, n.boot, band, cp.begin, alpha, iid) {

  # ****************************************************************************
  # In this script, functional prediction bands are calculated by adapting the
  # method described in Lenhoff et al. (1999).
  # ----------------------------------------------------------------------------
  #
  # The script is an adapted an slightly modified version of the MATLAB script
  # by Doris Oriwol, TU Chemnitz, 12.05.2010
  #
  # Function arguments:
  # k.coef   : Number of (Fourier) coefficients
  # n.boot   : Number of bootstrap iterations
  # band     : Type of interval (prediction or confidence)
  # cp.begin : Initial value quantile
  # alpha    : Significance level
  # iid      : If iid==TRUE, only one curve per subject is drawn, otherwise rep-
  #            ated measurement (several curves per subject) are allowed
  #            If TRUE, only one curve per subject is selected
  # ****************************************************************************

  # ----------------------------------------------------------------------------
  # Initial step: Get difference curves
  # ----------------------------------------------------------------------------
  if (iid == TRUE) {
    # Pick only one curve per subject to satisfy the iid assumption
    data.diff <- pick_curves(data, iid = TRUE)

  } else if (iid == FALSE) {
    # Pick 'n=number of subjects' strides (several strides per subject allowed)
    data.diff <- pick_curves(data, iid = FALSE)
  }

  # Dimensions
  n.time    <- dim(data.diff)[1]
  n.curves  <- dim(data.diff)[2]
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
    fourier.koeffi[, i] = pracma::mldivide(fourier.s, data.diff[, i])
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
  bootstrap.mean          <- matlab::zeros(k.coef*2 + 1, n.boot)
  bootstrap.real_mw       <- matlab::zeros(n.time, n.boot)
  bootstrap.zz            <- matlab::zeros(n.curves, n.boot)
  bootstrap.pseudo_koeffi <- matlab::zeros(k.coef*2 + 1, n.curves, n.boot)
  bootstrap.real          <- matlab::zeros(n.time, n.curves, n.boot)
  bootstrap.std1          <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, n.curves)
  bootstrap.kovarianz     <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, n.boot)
  bootstrap.std_all       <- matlab::zeros(n.time, n.time, n.boot)
  bootstrap.std           <- matlab::zeros(n.time, n.boot)

  for (i in 1:n.boot) {
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

  if (band == "prediction") {
    cp.data   <- matlab::zeros(n.curves, n.boot)
    cp.data_i <- matlab::zeros(n.curves, n.boot)

    cp.mean <- 0
    cp.bound <- cp.begin
    while (cp.mean < (1-alpha)) {
      for (i in 1:n.boot) {
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
      cc.data <- matlab::zeros(n.curves, n.boot)

      for (i in 1:n.boot) {
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

