pick_curves <- function (data, iid) {

  n.subj <- length(unique(data$subjectID))
  n.strides <- length(unique(data$strideID))

  if (iid == TRUE) {
    # Pick a single UNIQUE random stride from every subject
    curve.idx <- c()
    for (subj.idx in 1:n.subj) {

      if (subj.idx %in% unique(data$subjectID)) {
        tmp <- data[data$subjectID == subj.idx, ]

        curve.idx.subj <- as.numeric(sample(tmp$strideID, size = 1))
        curve.idx <- c(curve.idx, curve.idx.subj)
      }
    }
  } else {
    # Pick 'n=number of subjects' strides (several strides per subject allowed)
    curve.idx <- sample(n.strides, size = n.subj)
  }

  curve <- data[data$strideID %in% curve.idx, ]
  curve0 <- subset(curve, device  == "TWO")$value
  curve0 <- matrix(curve0, ncol = length(curve0) / 100)
  curve1 <- subset(curve, device  == "ONE")$value
  curve1 <- matrix(curve1, ncol = length(curve1) / 100)

  diff.curves <- curve0 - curve1
}
