points_within_limits <- function (curve, floa) {

  # ****************************************************************************
  # Calculate coverage (points of a single curve within the band limits)
  # ----------------------------------------------------------------------------
  #
  # 'curve' is a (single curve) numerical vector
  # 'floa' is a "matrix" "array" containing "lower" and "upper" limits
  # ****************************************************************************

  lwr.limit <- floa["lower.loa", ]
  upr.limit <- floa["upper.loa", ]

  outside <- 0

  # Compare difference curves with upper and lower boundaries
  below.thresh <- curve < lwr.limit
  above.thresh <- curve > upr.limit

  points.outside <- sum(above.thresh) + sum(below.thresh)

  coverage <- 1 - (points.outside / length(curve))
  coverage <- round(coverage, digits = 2)
}
