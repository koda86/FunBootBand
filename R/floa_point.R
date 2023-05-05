floa_point <- function(data) {

  # ----------------------------------------------------------------------------
  # Pointwise continuous Limits of Agreement
  #
  # Calculation described in Bland & Altman (2007) and Bland & Altman (1999)
  # Design currently implemented for balanced data
  # ----------------------------------------------------------------------------

  n.frames <- length(unique(data$frame))
  n.subjects <- length(unique(data$subjectID))
  n.strides <- length(unique(data$strideID))
  n.curves.per.subject <- n.strides / n.subjects

  mean.diff <- vector(mode = "list", length = n.frames)
  total.sd <- vector(mode = "list", length = n.frames)

  for (frame.idx in 0:100){
    data.by.frame <- subset(data, frame == frame.idx)

    # Get variance from ANOVA (see Bland & Altman, 2007) -----------------------
    # Prepare data for ANOVA
    data.by.frame.wide <- reshape2::dcast(data = data.by.frame, subjectID + strideID ~ device, value.var = "value")
    data.by.frame.wide$device.diff <- data.by.frame.wide$TWO - data.by.frame.wide$ONE
    data.by.frame.wide$subjectID <- as.factor(data.by.frame.wide$subjectID)

    LMEM.2 <- aov(device.diff ~ subjectID, data = data.by.frame.wide)
    S <- summary(LMEM.2)
    # Example from Bland & Altman (2007): total.sd <- sqrt((S[[1]]$`Mean Sq`[1] - S[[1]]$`Mean Sq`[2]) / ((60^2 - 312)/(11*60)) + S[[1]]$`Mean Sq`[2] )
    # ((60^2 - 312)/(11*60)) reduces to m (here: 10) because all subjects have the same number of observations
    # total.sd[[frame.idx + 1]] <- sqrt((S[[1]]$`Mean Sq`[1] - S[[1]]$`Mean Sq`[2]) / 10 + S[[1]]$`Mean Sq`[2])
    sigma.squared.within <- S[[1]]$`Mean Sq`[2]
    sigma.squared.between <- (S[[1]]$`Mean Sq`[1] - S[[1]]$`Mean Sq`[2]) / n.curves.per.subject
    total.sd[[frame.idx + 1]] <- sqrt(sigma.squared.between + sigma.squared.within)
    mean.diff[[frame.idx + 1]] <- mean(data.by.frame.wide$device.diff)
  }

  floa.point <- data.frame(unlist(mean.diff), unlist(total.sd))
  names(floa.point) <- c("mean.diff", "total.sd")


  # Calculate mean, upper and lower limits of agreement to get the same
  # structure as returned by the other methods ---------------------------------
  z0.975 <- qnorm(0.975, mean = 0, sd = 1)

  floa.point.loa <- rbind(floa.point$mean.diff + z0.975*floa.point$total.sd,
                          floa.point$mean.diff,
                          floa.point$mean.diff - z0.975*floa.point$total.sd
                          )


  # PRECISION OF ESTIMATED LIMITS OF AGREEMENT ---------------------------------
  # Confidence intervals for the 95% limits of agreement (Bland & Altman, 1999)
  n.curves <- length(unique(data$strideID))
  # Formula from Francq et al. (2020)
  standard.error.d2s <- sqrt(1/n.curves + (z0.975^2 / (2*(n.curves-1)))) * floa.point$total.sd

  # t-score of the 95th quantile of the Student t distribution with df = (n - 1)
  t <- qt(.975, df = n.curves - 1)

  # Upper and lower uncertainty limits -----------------------------------------
  upper.ci.upper <- floa.point$mean.diff + z0.975*floa.point$total.sd + t*standard.error.d2s
  upper.ci.lower <- floa.point$mean.diff + z0.975*floa.point$total.sd - t*standard.error.d2s
  lower.ci.upper <- floa.point$mean.diff - z0.975*floa.point$total.sd + t*standard.error.d2s
  lower.ci.lower <- floa.point$mean.diff - z0.975*floa.point$total.sd - t*standard.error.d2s

  floa.point.summary <- rbind(floa.point.loa, upper.ci.upper, upper.ci.lower, lower.ci.upper, lower.ci.lower)
  row.names(floa.point.summary) <- c("upper.loa", "mean", "lower.loa",
                                     "upper.ci.upper", "upper.ci.lower", "lower.ci.upper","lower.ci.lower")

  return(floa.point.summary)
}
