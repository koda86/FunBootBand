plot_loa <- function (data, floa.point, floa.roislien, floa.boot.rep, floa.boot.iid, ylim) {

  # ****************************************************************************
  # Plots the limits across the respective methods against all difference curves
  # in the data set
  # ----------------------------------------------------------------------------
  #
  # Function arguments:
  # - floa.point: point-wise limits of agreement (Bland & Altman, 2007)
  # - floa.roislien: limits of agreement as in Roislien et al. (2012)
  # - floa.boot: limits of agreement as in Lenhoff et al. (1999)
  # - ylim: y-axis limits (2 element numeric ob  ject)
  # ****************************************************************************

  # Prepare data for (gg)plotting
  floa.point    <- data.frame(t(floa.point))
  floa.roislien <- data.frame(t(floa.roislien))
  floa.boot.rep <- data.frame(t(floa.boot.rep))
  floa.boot.iid <- data.frame(t(floa.boot.iid))

  n.frames <- length(unique(data$frame))
  n.strides <- length(unique(data$strideID))
  n.subjects <- length(unique(data$subjectID))
  strides.per.subject <- length(unique(data$strideID)) / n.subjects

  # ----------------------------------------------------------------------------
  # Difference curves
  # ----------------------------------------------------------------------------

  device1 <- subset(data, device == "TWO", select = value)
  device2 <- subset(data, device == "ONE", select = value)

  device.diff <- device1 - device2

  colnames(device.diff)[1] <- "value"

  device.diff$frame <- seq(0, 100)
  device.diff$strideID <- as.factor(rep(1:n.strides,
                                        each = n.frames))
  device.diff$subjectID <- as.factor(rep(1:n.subjects,
                                         each = strides.per.subject * n.frames))

  # ----------------------------------------------------------------------------
  # Plot data
  # ----------------------------------------------------------------------------

  # For line graphs, the data points must be grouped so that it knows which points to connect.
  # In this case, it is simple -- all points should be connected, so group=1.
  # When more variables are used and multiple lines are drawn,
  # the grouping for lines is usually done by variable.
  PLOT.DIFF <- ggplot(data = device.diff, aes(x = frame, y = value)) +
    geom_line(aes(group = strideID),
              alpha = 0.1) +
    scale_color_grey(start = 0.8, end = 0.2) +
    # Plot lower limit of agreement
    geom_line(data = floa.boot.iid,
              aes(x = seq(0, 100), y = upper.loa, col = "red", group = 1),
              linetype = "solid",
              size = 1.2,
              colour = "#E69F00",
              alpha = 0.8) +
    geom_line(data = floa.boot.iid,
              aes(x = seq(0, 100), y = lower.loa, col = "red", group = 1),
              linetype = "solid",
              size = 1.2,
              colour = "#E69F00",
              alpha = 0.9) +
    geom_line(data = floa.boot.rep,
              aes(x = seq(0, 100), y = upper.loa, col = "red", group = 1),
              linetype = "solid",
              size = 1.2,
              colour = "royalblue1",
              alpha = 0.8) +
    geom_line(data = floa.boot.rep,
              aes(x = seq(0, 100), y = lower.loa, col = "red", group = 1),
              linetype = "solid",
              size = 1.2,
              colour = "royalblue1",
              alpha = 0.9) +
    geom_line(data = floa.roislien,
              aes(x = seq(0, 100), y = upper.loa),
              linetype = "solid",
              size = 1.2,
              colour = "deeppink",
              alpha = 1) +
    geom_line(data = floa.roislien,
              aes(x = seq(0, 100), y = lower.loa),
              linetype = "solid",
              size = 1.2,
              colour = "deeppink",
              alpha = 1) +
    geom_line(data = floa.point,
              aes(x = seq(0, 100), y = upper.loa),
              linetype = "dotted",
              size = 1.2,
              colour = "grey10",
              alpha = 0.8) +
    geom_line(data = floa.point,
              aes(x = seq(0, 100), y = lower.loa),
              linetype = "dotted",
              size = 1.2,
              colour = "grey10",
              alpha = 0.8) +
    scale_y_continuous(limits = c(ylim[1], ylim[2])) +
    labs(x = "Time-normalized signal [%]", y = "Difference") +
    theme_minimal()
    # theme(axis.text.x = element_text(size = 20),
    #       axis.title.x = element_text(size = 22),
    #       axis.text.y = element_text(size = 20),
    #       axis.title.y = element_text(size = 22),
    #       legend.position = "none")

  PLOT.DIFF
}
