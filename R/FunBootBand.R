#' @title FunBootBand
#'
#' @description Creates Functional Bootstraped (statistical) Bands.
#' IMPORTANT NOTE: Currently, the script is designed for balanced data sets.
#' Unbalanced designs (unequal number of curves) may lead to errors!
#'
#' @usage band(data, k.coef = 50, B = 10, type = "prediction", alpha = 0.05,
# iid = TRUE)
#'
#' @param data A data set consisting of n curves of length t. Needs to be a
#' numerical matrix of dimensions [t, n], e.g., data[1:101, 1:50] represents 50
#' curves of length 101 points.
#' @param type Band type (c("confidence", "prediction")).
#' @param B Number of bootstrap iterations (e.g., B = 1000).
#' @param iid Assume independent and identically distributed (iid) curves or not
#' (iid = c(TRUE, FALSE)). Setting iid=TRUE runs an ordinary (naive) bootstrap.
#' When setting iid=FALSE, a two-stage bootstrap is run, where clusters
#' (comprising all of their curves) are resampled with replacement in the
#' initial stage, and one curve per cluster is sampled without replacement in
#' the second stage. If iid is set to FALSE, curves are assumed to be nested in
#' curve clusters. These curve clusters need to be indicated in an additional
#' header line (see 'Format').
#' @param alpha Desired type I error probability.
#'
#' @return A data frame object that contains upper and lower band boundaries.
#' @examples
#' load("curvesample.RData")
#' band.limits <- band(data = curves, type = "prediction", B = 1000, iid = TRUE)
#' @export
#'

# TODO:
# - Funktion in C++ entwickeln
# - Vignette schreiben
# - Testen

# The header line needs to consist of letters.
# Technically requires stationary curves.

invisible(get(load("~/FunBootBand/data/curvesample.RData")))

band <- function(data, type, alpha, iid = TRUE, k.coef = 50, B = 400) {

  # Initial check if input arguments match the desired format
  if (class(type) != "character") {
    stop("'type' must be a variable of type 'character'.")
  } else if (!(type %in% c("confidence", "prediction"))) {
    stop("'type' must be either 'confidence' or 'prediction'.")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha > 1) {
    stop("'alpha' must be a numeric value between 0 and 1.")
  }
  if (!is.logical(iid)) {
    stop("'iid' must be a logical value (TRUE or FALSE).")
  }
  if (!is.numeric(k.coef) || k.coef <= 0) {
    stop("'k.coef' must be a positive integer.")
  }
  if (!is.numeric(B) || B <= 0) {
    stop("'B' must be a positive integer.")
  }

  if(all(is.na(suppressWarnings(as.numeric(data[1, ]))))) { # If header exists
    header <- as.character(data[1, ])
    data <- data[-1, ]
    n.time <- dim(data)[1]
    data <- matrix(as.numeric(data), nrow = n.time) # Reduce object size
  } else {
    n.time <- dim(data)[1]
    data <- matrix(as.numeric(data), nrow = n.time)
  }

  time <- seq(0, (n.time-1))
  n.curves  <- dim(data)[2]
  if (iid == FALSE) {
    n.cluster <- length(unique(header))
    curves.per.cluster <- n.curves / n.cluster
  }

  # Approximate curves using Fourier functions ---------------------------------
  fourier.koeffi    <- array(data = 0, dim = c(k.coef*2 + 1, n.curves))
  fourier.real      <- array(data = 0, dim = c(n.time, n.curves))
  fourier.mean      <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1))
  fourier.real_mw   <- array(data = 0, c(n.time, 1))
  fourier.std1      <- array(data = 0, c(k.coef*2 + 1, k.coef*2 + 1, n.curves))
  fourier.cov <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1))
  fourier.std_all   <- array(data = 0, c(n.time, n.time))
  fourier.std       <- array(data = 0, dim = c(n.time, 1))

  # Set up a Fourier series
  # General: f(t) = mu + sum(alpha cos(2pi*k*t/T) + beta sin(2pi*k*t/T))
  fourier.s = rep(1, times = n.time)
  for (k in seq(1, k.coef*2, 2)) {
    fourier.s <- cbind(fourier.s, cos(2*pi*(k/2)*time / (n.time-1)))
    fourier.s <- cbind(fourier.s, sin(2*pi*(k/2)*time / (n.time-1)))
  }

  for (i in 1:n.curves) {
    # Least squares Regression
    # fourier.koeffi[, i] = pracma::mldivide(fourier.s, data[, i]) # old
    fourier.koeffi[, i] = pinv(t(fourier.s) %*% fourier.s) %*% t(fourier.s) %*% data[, i]
    # Fourier curve
    fourier.real[, i] = fourier.s %*% fourier.koeffi[, i]
  }

  # Mean Fourier curve
  fourier.mean[, 1] = rowMeans(fourier.koeffi)
  fourier.real_mw[, 1] = fourier.s %*% fourier.mean[, 1]

  # Standard deviation of the Fourier curve
  for (i in 1:n.curves) {
    # Variance-covariance matrix
    fourier.std1[, , i] <- (fourier.koeffi[, i] - fourier.mean[, 1]) %*%
                           t(fourier.koeffi[, i] - fourier.mean[, 1])
  }

  fourier.cov <- apply(fourier.std1, c(1, 2), mean)
  # Lenhoff, Appendix A, Eq. (0.5)
  fourier.std_all <- suppressWarnings(sqrt(fourier.s %*% fourier.cov %*%
                     t(fourier.s))
                     )

  for (i in 1:n.time) {
    # Values are on the diagonal of the square matrix fourier.std_all
    fourier.std[i, 1] = fourier.std_all[i, i]
  }

  # Bootstrap ------------------------------------------------------------------
  bootstrap_sample        <- array(data = 0, dim = c(n.time, 4))
  bootstrap.mean          <- array(data = 0, dim = c(k.coef*2 + 1, B))
  bootstrap.real_mw       <- array(data = 0, dim = c(n.time, B))
  bootstrap.zz            <- array(data = 0, dim = c(n.curves, B))
  bootstrap.pseudo_koeffi <- array(data = 0, dim = c(k.coef*2 + 1, n.curves, B))
  bootstrap.real          <- array(data = 0, dim = c(n.time, n.curves, B))
  bootstrap.std1          <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1, n.curves))
  bootstrap.cov     <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1, B))
  bootstrap.std_all       <- array(data = 0, dim = c(n.time, n.time, B))
  bootstrap.std           <- array(data = 0, dim = c(n.time, B))

  for (i in 1:B) {
    if (iid == FALSE) {
      for (k in 1:curves.per.cluster) {
        # STAGE 1: Sample curve clusters with replacement
        stage.1.idx <- sample(1:n.cluster, size = n.cluster, replace = TRUE)
        # STAGE 2: Sample within stage clusters without replacement
        curves <- c()
        for (curve.idx in stage.1.idx) {
          curve.numbers.stage.1 <- seq(from = curve.idx*curves.per.cluster -
                                          curves.per.cluster + 1,
                                        to = curve.idx*curves.per.cluster)
          tmp <- sample(curve.numbers.stage.1, size = 1, replace = FALSE)
          while (tmp %in% curves) { # Assure drawing without replacement
            tmp <- sample(curve.numbers.stage.1, size = 1)
          }
          curves <- c(curves, tmp)
        }
        bootstrap.zz[k, i] = curves[k]
        bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
        bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
    }
    } else {
      for (k in 1:n.curves) {
        bootstrap.zz[k, i] = sample(n.curves, size=1)
        bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
        bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
      }
    }

    # Mean bootstrap curve and standard deviation
    bootstrap.mean[, i] <- rowMeans(bootstrap.pseudo_koeffi[, , 1])
    bootstrap.real_mw[, i] <- fourier.s %*% bootstrap.mean[, i]

    for (k in 1:n.curves) {
      bootstrap.std1[, , k] <- (bootstrap.pseudo_koeffi[, k, i] -
                                  bootstrap.mean[, i]) %*%
                        t(bootstrap.pseudo_koeffi[, k, i] - bootstrap.mean[, i])
    }

    bootstrap.cov[, , i] <- apply(bootstrap.std1, c(1, 2), mean)
    bootstrap.std_all[, , i] <- suppressWarnings(sqrt(fourier.s %*%
                                bootstrap.cov[, , i] %*% t(fourier.s))
                                )

    for (k in 1:n.time) {
      bootstrap.std[k, i] <- bootstrap.std_all[k, k, i]
    }
  }

  # Construct bands ------------------------------------------------------------
  band.mean <- rowMeans(bootstrap.real_mw)
  band.sd   <- rowMeans(bootstrap.std)

  if (type == "prediction") {
    cp.data   <- array(data = 0, dim = c(n.curves, B))
    cp.data_i <- array(data = 0, dim = c(n.curves, B))

    cp.mean <- 0
    cp.bound <- 0 # cp.begin
    while (cp.mean < (1-alpha)) {
      for (i in 1:B) {
        for (k in 1:n.curves) {
          # Lenhoff et al., Appendix A, Eq. (0.6)
          cp.data[k, i] <- max(abs(fourier.real[, k] - bootstrap.real_mw[, i]) /
                                 bootstrap.std[, i])
          cp.data_i[k, i] <- cp.data[k, i] < cp.bound
        }
      }
      cp.mean <- mean(cp.data_i)
      cp.bound <- cp.bound + 0.05
    }
    cp_out <- cp.bound

    band.boot <- rbind(band.mean + cp.bound * band.sd,
                       band.mean,
                       band.mean - cp.bound * band.sd
    )
  } else if (type == "confidence") {
    cc.data <- array(data = 0, dim = c(n.curves, B))

    for (i in 1:B) {
      for (k in 1:n.curves) {
        # Lenhoff, Appendix A, Eq. (0.8)
        cc.data[k, i] <- max(abs(bootstrap.real_mw[, i] - fourier.real_mw) /
                               bootstrap.std[, i])
      }
    }
    cc <- quantile(cc.data, probs = 1-alpha)

    band.boot <- rbind(band.mean + cc * band.sd,
                       band.mean,
                       band.mean - cc * band.sd
    )
  }

  row.names(band.boot) <- c("upper", "mean", "lower")

  return(band.boot)
}

