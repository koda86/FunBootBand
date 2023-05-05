estimate_uncertainty_kfold_rep <- function (data, n.rep, n.boot, plot.au, plot.ylimits) {

  # ****************************************************************************
  # Estimate uncertainty in different methods (percentile intervals)
  # ----------------------------------------------------------------------------
  #
  # Function arguments:
  # data        : data
  # n.rep       : Number of repeated calculations
  # n.boot      : Number of bootstrap iterations
  # plot.au     : TRUE includes a plot of the uncertainty ribbons
  # plot.ylimits: (manually set) y-axis limits for plots in different data sets
  # ****************************************************************************

  n.subj <- length(unique(data$subjectID))

  # ----------------------------------------------------------------------------
  # Calculate uncertainty area across repeated calculations (n.rep times)
  # ----------------------------------------------------------------------------
  floa.roislien.upper <- vector(mode = "list", length = n.boot)
  floa.roislien.lower <- vector(mode = "list", length = n.boot)
  floa.boot.rep.upper <- vector(mode = "list", length = n.boot)
  floa.boot.rep.lower <- vector(mode = "list", length = n.boot)
  floa.boot.iid.upper <- vector(mode = "list", length = n.boot)
  floa.boot.iid.lower <- vector(mode = "list", length = n.boot)

  rep.idx <- 1


  while (rep.idx < n.rep) {

    print(rep.idx)

    for (i in 1:n.subj) {

      print(i)

      # Cross validation: In each repetition, one (random) subject is left out
      data.leaveout <- data[data$subjectID != i, ]

      floa.point     <- floa_point(data.leaveout)
      floa.roislien  <- floa_roislien(data.leaveout)
      floa.boot.rep  <- floa_boot_rep(data.leaveout,
                                      k.coef = 50,
                                      n.boot = n.boot,
                                      band = "prediction",
                                      cp.begin = 0,
                                      alpha = 0.05)
      floa.boot.iid  <- floa_boot(data.leaveout,
                                  k.coef = 50,
                                  n.boot = n.boot,
                                  band = "prediction",
                                  cp.begin = 0,
                                  alpha = 0.05,
                                  iid = TRUE) # Draw only a single curve per subject

      # Use plot to check if the method works correctly
      # if (rep.idx == 1 & i == 1) {
      #   plot(floa.point["upper.loa", ], type = "l", ylim = c(-1, 1))
      # } else {
      #   lines(floa.point["upper.loa", ])
      # }

      floa.roislien.upper[[i]] <- floa.roislien["upper.loa", ]
      floa.roislien.lower[[i]] <- floa.roislien["lower.loa", ]
      floa.boot.rep.upper[[i]]  <- floa.boot.rep["upper.loa", ]
      floa.boot.rep.lower[[i]]  <- floa.boot.rep["lower.loa", ]
      floa.boot.iid.upper[[i]]  <- floa.boot.iid["upper.loa", ]
      floa.boot.iid.lower[[i]]  <- floa.boot.iid["lower.loa", ]
    }
    rep.idx <- rep.idx + 1
  }

  floa.roislien.upper.unlist <- matrix(unlist(floa.roislien.upper), 101)
  floa.roislien.lower.unlist <- matrix(unlist(floa.roislien.lower), 101)
  floa.boot.rep.upper.unlist <- matrix(unlist(floa.boot.rep.upper), 101)
  floa.boot.rep.lower.unlist <- matrix(unlist(floa.boot.rep.lower), 101)
  floa.boot.iid.upper.unlist <- matrix(unlist(floa.boot.iid.upper), 101)
  floa.boot.iid.lower.unlist <- matrix(unlist(floa.boot.iid.lower), 101)

  # Pointwise LoA return the same estimate every time (all curves are included)
  floa.point <- data.frame(t(floa_point(data.leaveout)))

  # Pointwise 95% percentile intervals for the distribution of the limits
  # The percentiles form the boundaries of the undertainty area
  pi95 <- rbind(
    apply(floa.roislien.upper.unlist, 1, quantile, probs = c(0.025, 0.975)),
    apply(floa.roislien.lower.unlist, 1, quantile, probs = c(0.025, 0.975)),
    apply(floa.boot.rep.upper.unlist, 1, quantile, probs = c(0.025, 0.975)),
    apply(floa.boot.rep.lower.unlist, 1, quantile, probs = c(0.025, 0.975)),
    apply(floa.boot.iid.upper.unlist, 1, quantile, probs = c(0.025, 0.975)),
    apply(floa.boot.iid.lower.unlist, 1, quantile, probs = c(0.025, 0.975)),
    # Pointwise limits are calculated differently (using CI, see floa_point())
    # but used in the same way here to facilitate plotting (plausibility control)
    floa.point$upper.ci.lower,
    floa.point$upper.ci.upper,
    floa.point$lower.ci.lower,
    floa.point$lower.ci.upper
  )

  pi95 <- data.frame(t(pi95))
  colnames(pi95) <- c(
    "floa.roislien.upper.2.5", "floa.roislien.upper.97.5", "floa.roislien.lower.2.5", "floa.roislien.lower.97.5",
    "floa.boot.rep.upper.2.5", "floa.boot.rep.upper.97.5", "floa.boot.rep.lower.2.5", "floa.boot.rep.lower.97.5",
    "floa.boot.iid.upper.2.5", "floa.boot.iid.upper.97.5", "floa.boot.iid.lower.2.5", "floa.boot.iid.lower.97.5",
    "floa.point.upper.2.5", "floa.point.upper.97.5", "floa.point.lower.2.5", "floa.point.lower.97.5"
  )



  # ----------------------------------------------------------------------------
  # Plot uncertainty areas
  # ----------------------------------------------------------------------------
  if (plot.au == TRUE) {

    # Prepare data for ggploting -------------------------------------------------
    device1 <- data.frame(subset(data.leaveout, device == "TWO")$value)
    device2 <- data.frame(subset(data.leaveout, device == "ONE")$value)
    device.diff <- device1 - device2
    colnames(device.diff)[1] <- "value"

    device.diff$frame <- seq(0, 100)
    n.frames <- length(unique(data.leaveout$frame))
    n.strides <- length(unique(data.leaveout$strideID))
    n.subjects <- length(unique(data.leaveout$subjectID))
    strides.per.subject <- length(unique(data.leaveout$strideID)) / n.subjects
    device.diff$strideID <- as.factor(rep(1:n.strides, each = 101))
    device.diff$subjectID <- as.factor(rep(1:n.subjects, each = strides.per.subject * n.frames))

    # ROISLIEN
    device.diff$floa.roislien.upper.2.5  <- rep(pi95$floa.roislien.upper.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.roislien.upper.97.5 <- rep(pi95$floa.roislien.upper.97.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.roislien.lower.2.5  <- rep(pi95$floa.roislien.lower.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.roislien.lower.97.5 <- rep(pi95$floa.roislien.lower.97.5,
                                                length.out = nrow(device.diff))
    # BOOTrep
    device.diff$floa.boot.rep.upper.2.5  <- rep(pi95$floa.boot.rep.upper.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.boot.rep.upper.97.5 <- rep(pi95$floa.boot.rep.upper.97.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.boot.rep.lower.2.5  <- rep(pi95$floa.boot.rep.lower.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.boot.rep.lower.97.5 <- rep(pi95$floa.boot.rep.lower.97.5,
                                                length.out = nrow(device.diff))
    # BOOTiid
    device.diff$floa.boot.iid.upper.2.5  <- rep(pi95$floa.boot.iid.upper.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.boot.iid.upper.97.5 <- rep(pi95$floa.boot.iid.upper.97.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.boot.iid.lower.2.5  <- rep(pi95$floa.boot.iid.lower.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.boot.iid.lower.97.5 <- rep(pi95$floa.boot.iid.lower.97.5,
                                                length.out = nrow(device.diff))

    # POINT (based on confidence intervals around the band limits)
    device.diff$floa.point.upper.2.5     <- rep(pi95$floa.point.upper.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.point.upper.97.5    <- rep(pi95$floa.point.upper.97.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.point.lower.2.5     <- rep(pi95$floa.point.lower.2.5,
                                                length.out = nrow(device.diff))
    device.diff$floa.point.lower.97.5    <- rep(pi95$floa.point.lower.97.5,
                                                length.out = nrow(device.diff))

    # Plot Matrix ----------------------------------------------------------------

    # For line graphs, the data points must be grouped so that it knows which points to connect.
    # In this case, it is simple -- all points should be connected, so group=1.
    # When more variables are used and multiple lines are drawn, the grouping for
    # lines is usually done by variable.
    PLOT.BOOTiid <- ggplot(data = device.diff, aes(x = frame, y = value)) +
      geom_ribbon(aes(ymin = floa.boot.iid.upper.2.5,
                      ymax = floa.boot.iid.upper.97.5),
                  fill = "#E69F00",
                  alpha = 0.5) +
      geom_ribbon(aes(ymin = floa.boot.iid.lower.2.5,
                      ymax = floa.boot.iid.lower.97.5),
                  fill = "#E69F00",
                  alpha = 0.5) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.15) +
      theme_minimal() +
      ylim(plot.ylimits) +
      labs(x = "Time-normalized signal [%]", y = "Difference")


    PLOT.BOOTrep <- ggplot(data = device.diff, aes(x = frame, y = value)) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.05) +
      geom_ribbon(aes(ymin = floa.boot.rep.upper.2.5,
                      ymax = floa.boot.rep.upper.97.5),
                  fill = "royalblue1",
                  alpha = 0.5) +
      geom_ribbon(aes(ymin = floa.boot.rep.lower.2.5,
                      ymax = floa.boot.rep.lower.97.5),
                  fill = "royalblue1",
                  alpha = 0.5) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.15) +
      theme_minimal() +
      ylim(plot.ylimits) +
      labs(x = "Time-normalized signal [%]", y = "Difference")


    PLOT.POINT <- ggplot(data = device.diff, aes(x = frame, y = value)) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.05) +
      geom_ribbon(aes(ymin = floa.point.upper.2.5,
                      ymax = floa.point.upper.97.5),
                  fill = "grey50",
                  alpha = 0.5) +
      geom_ribbon(aes(ymin = floa.point.lower.2.5,
                      ymax = floa.point.lower.97.5),
                  fill = "grey50",
                  alpha = 0.5) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.15) +
      theme_minimal() +
      ylim(plot.ylimits) +
      labs(x = "Time-normalized signal [%]", y = "Difference")

    # Plot version with ribbon
    # PLOT.POINT <- ggplot(data = device.diff, aes(x = frame, y = value)) +
    #   geom_line(data = device.diff,
    #             aes(group = strideID),
    #             alpha = 0.05) +
    #   geom_line(data = floa.point,
    #             aes(x = seq(0, 100), y = upper.loa),
    #             linetype = "dotted",
    #             size = 1) +
    #   geom_line(data = floa.point,
    #             aes(x = seq(0, 100), y = lower.loa),
    #             linetype = "dotted",
    #             size = 1) +
    #   geom_line(data = device.diff,
    #             aes(group = strideID),
    #             alpha = 0.15) +
    #   theme_minimal() +
    #   theme(axis.text.x = element_text(size = 15),
    #         axis.title.x = element_text(size = 17),
    #         axis.text.y = element_text(size = 15),
    #         axis.title.y = element_text(size = 17),
    #         legend.position = "none") +
    #   ylim(-5, 5) +
    #   labs(x = "Time-normalized signal [%]", y = "Difference")


    PLOT.ROISLIEN <- ggplot(data = device.diff, aes(x = frame, y = value)) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.05) +
      geom_ribbon(aes(ymin = floa.roislien.upper.2.5,
                      ymax = floa.roislien.upper.97.5),
                  fill = "deeppink",
                  alpha = 0.5) +
      geom_ribbon(aes(ymin = floa.roislien.lower.2.5,
                      ymax = floa.roislien.lower.97.5),
                  fill = "deeppink",
                  alpha = 0.5) +
      geom_line(data = device.diff,
                aes(group = strideID),
                alpha = 0.15) +
      theme_minimal() +
      ylim(plot.ylimits) +
      labs(x = "Time-normalized signal [%]", y = "Difference")


    PLOT <- ggpubr::ggarrange(PLOT.BOOTiid, PLOT.BOOTrep, PLOT.POINT, PLOT.ROISLIEN,
                              labels = c("A", "B", "C", "D"),
                              ncol = 2,
                              nrow = 2)

    # ggsave(plot = PLOT,
    #        filename = "~/Nextcloud/project-fab-forschung/Publikationen/FLOA/paper_JournalBiomechanics/Review/tex/Grafiken/au_kfold.png",
    #        bg = "white",
    #        dpi = 300)

    # # # Plot all in one ------------------------------------------------------------
    # # For line graphs, the data points must be grouped so that it knows which points to connect.
    # # In this case, it is simple -- all points should be connected, so group=1.
    # # When more variables are used and multiple lines are drawn, the grouping for lines is usually done by variable.
    # PLOT.DIFF <- ggplot(data = device.diff, aes(x = frame, y = value)) +
    #   # geom_line(data = device.diff,
    #   #           aes(group = strideID),
    #   #           alpha = 0.05) +
    #   geom_ribbon(aes(ymin = floa.boot.iid.upper.2.5,
    #                   ymax = floa.boot.iid.upper.97.5),
    #               fill = "#E69F00",
    #               alpha = 0.5) +
    #   geom_ribbon(aes(ymin = floa.boot.iid.lower.2.5,
    #                   ymax = floa.boot.iid.lower.97.5),
    #               fill = "#E69F00",
    #               alpha = 0.5) +
    #   geom_ribbon(aes(ymin = floa.boot.rep.upper.2.5,
    #                 ymax = floa.boot.rep.upper.97.5),
    #             fill = "royalblue1",
    #             alpha = 0.5) +
    #   geom_ribbon(aes(ymin = floa.boot.rep.lower.2.5,
    #                   ymax = floa.boot.rep.lower.97.5),
    #               fill = "royalblue1",
    #               alpha = 0.5) +
    #   geom_line(data = floa.point,
    #             aes(x = seq(0, 100), y = upper.loa),
    #             linetype = "dotted",
    #             size = 1) +
    #   geom_line(data = floa.point,
    #             aes(x = seq(0, 100), y = lower.loa),
    #             linetype = "dotted",
    #             size = 1) +
    #   geom_ribbon(aes(ymin = floa.roislien.upper.2.5,
    #                   ymax = floa.roislien.upper.97.5),
    #               fill = "deeppink",
    #               alpha = 0.5) +
    #   geom_ribbon(aes(ymin = floa.roislien.lower.2.5,
    #                   ymax = floa.roislien.lower.97.5),
    #               fill = "deeppink",
    #               alpha = 0.5) +
    #   theme_minimal() +
    #   theme(axis.text.x = element_text(size = 15),
    #         axis.title.x = element_text(size = 17),
    #         axis.text.y = element_text(size = 15),
    #         axis.title.y = element_text(size = 17),
    #         legend.position = "none") +
    #   ylim(-7, 7) +
    #   labs(x = "Time-normalized signal [%]", y = "Difference")
    #
    # PLOT.DIFF
  }



  # ----------------------------------------------------------------------------
  # Area between 2.5 and 97.5 percentile borders (Cumulated Uncertainty Area)
  # ----------------------------------------------------------------------------
  # POINT values are identical for the upper and lower limits
  CUA.point.upper     <- round(sum(pi95$floa.point.upper.97.5 - pi95$floa.point.upper.2.5), 2)
  CUA.point.lower     <- round(sum(pi95$floa.point.lower.97.5 - pi95$floa.point.lower.2.5), 2)
  CUA.roislien.upper  <- round(sum(pi95$floa.roislien.upper.97.5 - pi95$floa.roislien.upper.2.5), 2)
  CUA.roislien.lower  <- round(sum(pi95$floa.roislien.lower.97.5 - pi95$floa.roislien.lower.2.5), 2)
  CUA.boot.rep.upper  <- round(sum(pi95$floa.boot.rep.upper.97.5 - pi95$floa.boot.rep.upper.2.5), 2)
  CUA.boot.rep.lower  <- round(sum(pi95$floa.boot.rep.lower.97.5 - pi95$floa.boot.rep.lower.2.5), 2)
  CUA.boot.iid.upper  <- round(sum(pi95$floa.boot.iid.upper.97.5 - pi95$floa.boot.iid.upper.2.5), 2)
  CUA.boot.iid.lower  <- round(sum(pi95$floa.boot.iid.lower.97.5 - pi95$floa.boot.iid.lower.2.5), 2)

  uncert.area <- c(CUA.point.upper + CUA.point.lower,
                   CUA.roislien.upper + CUA.roislien.upper,
                   CUA.boot.rep.upper + CUA.boot.rep.lower,
                   CUA.boot.iid.upper + CUA.boot.iid.lower
  )

  names(uncert.area) <- c("POINT", "ROISLIEN", "BOOT.rep", "BOOT.iid")

  return(uncert.area)
}

