coverage_curves <- function(vector.cover, p.cover) {

  # ****************************************************************************
  # What is the proportion of bands (across cross validation iterations) that
  # contain at least p.cover of the T = 101 curve points?
  # ----------------------------------------------------------------------------
  #
  # p.cover (0 < p.cover < 1) is the proportion of covered points per curve
  # vector.cover is a column vector with coverage percentages for each curve
  # ****************************************************************************

  proportion.covered.points <- sum(vector.cover >= p.cover) / length(vector.cover)
  round(proportion.covered.points, 2)
}
