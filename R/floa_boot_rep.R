floa_boot_rep <- function(data, k.coef, n.boot, band, cp.begin, alpha) {
  
  # ****************************************************************************
  # In this script, functional prediction bands are calculated by adapting the
  # method described in Lenhoff et al. (1999).
  #
  # # The variant implemented here corresponds to a two-stage bootstrap as
  # described in (Davison and Hinkely, 1997).
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
  # ****************************************************************************
  
  # ----------------------------------------------------------------------------
  # Initial step: Get difference curves
  # ----------------------------------------------------------------------------
  devices <- unique(data$device)

  device.1      <- subset(data, device  == devices[1])$value
  device.1.wide <- matrix(device.1, ncol = length(device.1) / 101)
  device.2      <- subset(data, device == devices[2])$value
  device.2.wide <- matrix(device.2, ncol = length(device.2) / 101)

  data.diff <- device.1.wide - device.2.wide
  
  # Dimensions
  n.time             <- dim(data.diff)[1]
  n.subj             <- length(unique(data$subjectID))
  n.curves           <- dim(data.diff)[2]
  curves.per.subject <- n.curves / n.subj
  n.strides          <- length(unique(data$strideID))
  time               <- seq(0, (n.time-1))
  
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
  bootstrap.zz            <- matlab::zeros(curves.per.subject, n.boot)
  bootstrap.pseudo_koeffi <- matlab::zeros(k.coef*2 + 1, curves.per.subject, n.boot)
  bootstrap.real          <- matlab::zeros(n.time, curves.per.subject, n.boot)
  bootstrap.std1          <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, curves.per.subject)
  bootstrap.kovarianz     <- matlab::zeros(k.coef*2 + 1, k.coef*2 + 1, n.boot)
  bootstrap.std_all       <- matlab::zeros(n.time, n.time, n.boot)
  bootstrap.std           <- matlab::zeros(n.time, n.boot)
    
  for (i in 1:n.boot) {
    # The variant implemented here corresponds to a two-stage bootstrap as
    # described in (Davison and Hinkely, 1997):
    # 
    # There are two simple strategies, for both of which the first stage is to
    # randomly sample groups with replacement. At the second stage we randomly
    # sample within the groups selected at the first stage, either without
    # replacement (Strategy 1) or with replacement (Strategy 2).
    # 
    # Note that Strategy 1 keeps selected groups intact.
    # ...
    # On balance, therefore, Strategy 1 more closely mimics the variation
    # properties of the data, and so is the preferable strategy.
    #
    # Therefore the variant implemented here is:
    # 1. stage: Randomly sample groups WITH replacement
    # 2. stage: Sample within these groups (stage 1) WITHOUT replacement
    
    for (k in 1:curves.per.subject) {
      # STAGE 1: Randomly sample /subjects WITH replacement
      stage.1.idx <- sample(1:n.subj, size = n.subj, replace = TRUE)
      
      # STAGE 2: Sample within stage 1 groups (subjects) WITHOUT replacement
      strides <- c()
      for (stride.idx in stage.1.idx) {
        stride.numbers.stage.1 <- seq(from = stride.idx*curves.per.subject - curves.per.subject + 1, to = stride.idx*curves.per.subject)
        tmp <- sample(stride.numbers.stage.1, size = 1, replace = FALSE)
        while (tmp %in% strides) { # Assure drawing WITHOUT replacement
          tmp <- sample(stride.numbers.stage.1, size = 1)
        }
        strides <- c(strides, tmp)
      }
      
      bootstrap.zz[k, i] = strides[k]
      bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
      bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
    }
  
    # Mean bootstrap curve and standard deviation
    bootstrap.mean[, i] <- rowMeans(bootstrap.pseudo_koeffi[, , 1])
    bootstrap.real_mw[, i] <- fourier.s %*% bootstrap.mean[, i]
    
    for (k in 1:curves.per.subject) {
      bootstrap.std1[, , k] <- (bootstrap.pseudo_koeffi[, k, i] - bootstrap.mean[, i]) %*% t(bootstrap.pseudo_koeffi[, k, i] - bootstrap.mean[, i])
    }
    
    bootstrap.kovarianz[, , i] <- apply(bootstrap.std1, c(1, 2), mean)
    bootstrap.std_all[, , i] <- suppressWarnings(sqrt(fourier.s %*% bootstrap.kovarianz[, , i] %*% t(fourier.s)))
    
    for (k in 1:n.time) {
      bootstrap.std[k, i] <- bootstrap.std_all[k, k, i]
    }
    
    # Dummy plots (Appendix) (which curves are drawn per bootstrap iteration?)
    # condition <- "rep_new"
    # filename.dummy <- paste0("~/Nextcloud/project-fab-forschung/Publikationen/FLOA/paper_JournalBiomechanics/Review/dummy_plots/",
    #                          condition,
    #                          i, # bootstrap iteration
    #                          ".png")
    # png(filename.dummy)
    # plot(bootstrap.real[, 1, 1],
    #      type = "l",
    #      ylim = c(-2, 2),
    #      main = paste0("BOOTrep (iteration ", i, "), SD at 87%: ", round(bootstrap.std[, i][87], 2)),
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
    # n <- 1
    while (cp.mean < (1-alpha)) {
      for (i in 1:n.boot) {
        for (k in 1:n.curves) {
          # Lenhoff, Appendix A, Eq. (0.6)
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
    cc.data   <- matlab::zeros(n.curves, n.boot)
    cc.data_i <- matlab::zeros(n.curves, n.boot)
    
    # Determine the quantile for alpha
    cc.mean <- 0
    cc.bound <- cp.begin
    # n <- 1
    
    while (cc.mean < (1 - alpha)) {
      for (i in 1:n.boot) {
        for (k in 1:n.curves) {
          # Lenhoff, Appendix A, Eq. (0.7)
          cc.data[k, i] <- max(abs(bootstrap.real_mw[, i] - fourier.real_mw) / bootstrap.std[, i])
          cc.data_i[k, i] <- cc.data[k, i] < cc.bound
        }
      }
      cc.mean <- mean(cc.data_i)
      cc.bound <- cc.bound + 0.05
    }
    
    cc_out <- cc.bound
    
    # Construct bands
    # ----------------------------------------------------------------------
    floa.boot <- rbind(floa.boot.mean + cc.bound * floa.boot.sd,
                       floa.boot.mean,
                       floa.boot.mean - cc.bound * floa.boot.sd
    )
    }

  row.names(floa.boot) <- c("upper.loa", "mean", "lower.loa")
  
  return(floa.boot)
}

