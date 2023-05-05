floa_roislien <- function(data) {

  # ****************************************************************************
  # Functional limits of agreement according to Roislien et al. (2012)
  #
  # Currently no functional data methods are included.
  # ****************************************************************************

  # Pick only one curve per subject to satisfy iid assumption
  diff.curves <- pick_curves(data, iid = TRUE)

  # Calculate Limits of Agreement
  func.mean <- apply(diff.curves, 1, mean)
  func.sd <- apply(diff.curves, 1, sd)

  # Get the same structure as returned by the other methods
  floa.roislien <- rbind(func.mean + 1.96 * func.sd,
                         func.mean,
                         func.mean - 1.96 * func.sd)

  rownames(floa.roislien) <- c("upper.loa", "mean", "lower.loa")

  return(floa.roislien)
}
