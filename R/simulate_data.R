################################################################################
# Simulate data
################################################################################
#
# Simulated data sets:
#
# 1. Realistic smooth wave data (normal error, constant variance, no trend)
# 2. Same model as in (1), but almost no variation in one of the devices: "smooth"
# 3. Heteroscedastic error
# 4. Non-stationary data (trend, no bias) data:"non_stationary"
# 5. Non-gaussian (Weibull distributed) error (no trend)
# 6. Data with shock peaks (no bias, no trend): "shock"
# 7. Shift data (data phase shifted (x-axis direction))
################################################################################

n.subj <- 11
n.strides <- 10
n.devices <- 2
n.frames <- 101

# ------------------------------------------------------------------------------
# 1. Realistic smooth wave data (constant variance, no trend)
# ------------------------------------------------------------------------------
set.seed(3)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.5, max = 0.5)

  a.sd <- runif(1, min = 0.05, max = 0.15)
  b.sd <- runif(1, min = 0.0001, max = 0.002)

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(1, mean = 3, sd = b.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd)
    c <- 2

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.1 * t)

    one <- rnorm(1, offset.mean, 0.05) + sine.1 + sine.2
    two <- rnorm(1, offset.mean, 0.05) + sine.3 + sine.4

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "smooth_realistic.rds"))
# write.table(data, "~/floa/R/examples/smooth_realistic.txt", sep="\t")



# ------------------------------------------------------------------------------
# 2. Smooth wave data (constant variance, no trend) with subjectwise bias
#
# - Same model as in (1), but almost no variation in one of the devices
# ------------------------------------------------------------------------------
set.seed(3)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.5, max = 0.5)

  a.sd <- runif(1, min = 0.05, max = 0.15)
  b.sd <- runif(1, min = 0.0001, max = 0.002)
  b.sd.1 <- runif(1, min = 0.005, max = 0.01)

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(1, mean = 3, sd = a.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd)
    c <- 2

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.2 * t)

    offset <- rnorm(1, offset.mean, 0.05)

    one <- sine.1 + sine.2
    two <- offset + sine.3 + sine.4

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "smooth.rds"))
# write.table(data, "~/floa/R/examples/smooth.txt", sep="\t")



# ------------------------------------------------------------------------------
# 3. Heteroscedasticity (no trend)
# ------------------------------------------------------------------------------
set.seed(2)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.5, max = 0.5)

  a.sd <- runif(1, min = 0.05, max = 0.15)
  b.sd <- runif(1, min = 0.0001, max = 0.002)
  b.sd.1 <- runif(1, min = 0.005, max = 0.01)

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(1, mean = 3, sd = a.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd.1)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd.1)
    c <- 2

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.2 * t)

    offset <- rnorm(1, offset.mean, 0.05)

    one <- sine.1 + sine.2
    two <- offset + sine.3 + sine.4

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "heteroscedastic.rds"))
# write.table(data, "~/floa/R/examples/heteroscedastic.txt", sep="\t")



# ------------------------------------------------------------------------------
# 4. Smooth wave data with nonlinear trend
# ------------------------------------------------------------------------------
set.seed(3)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.5, max = 0.5)

  a.sd <- runif(1, min = 0.05, max = 0.15)
  b.sd <- runif(1, min = 0.0001, max = 0.002)

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(1, mean = 3, sd = a.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd)
    c <- 2

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.1 * t)

    # Create nonlinear trend
    trend <- (1/100000) * seq(0.5, 50.5, 0.5) ^ 3

    one <- rnorm(1, offset.mean, 0.05) + sine.1 + sine.2
    two <- rnorm(1, offset.mean, 0.05) + sine.3 + sine.4 + trend

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "smooth_trend.rds"))
# write.table(data, "~/floa/R/examples/smooth_trend.txt", sep="\t")



# ------------------------------------------------------------------------------
# 5. Non-gaussian (Weibull distributed) error (no trend)
# ------------------------------------------------------------------------------
set.seed(3)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.5, max = 0.5)

  a.sd <- runif(1, min = 0.05, max = 0.15)
  b.sd <- runif(1, min = 0.0001, max = 0.002)

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(n = 1,
                  # factorial() used to center around 0
                  mean = rweibull(1, shape = 1.5, scale=1) - factorial(1/1.5),
                  sd = a.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd)
    c <- 2

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.1 * t)

    offset <- rnorm(1, offset.mean, 0.05)

    one <- rnorm(1, offset.mean, 0.05) + sine.1 + sine.2
    two <- rnorm(1, offset.mean, 0.05) + sine.3 + sine.4

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "non_gaussian.rds"))
# write.table(data, "~/floa/R/examples/non_gaussian.txt", sep="\t")



# ------------------------------------------------------------------------------
# 6. Shock (spike) data
# ------------------------------------------------------------------------------
set.seed(6)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.005, max = 0.005)

  a.sd <- runif(1, min = 0.05, max = 0.15) # .1
  b.sd <- runif(1, min = 0.0001, max = 0.002) # .001

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(1, mean = 3, sd = a.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd)
    c <- 2

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.1 * t)

    offset <- rnorm(1, offset.mean, 0.01)

    one <- sine.1 + sine.2
    two <- sine.3 + sine.4 + offset

    # Add shock to the curves of a subject
    spike.idx <- round(runif(1, min = 6, max = 7))
    # spike.idx.one <- round(runif(1, min = 5, max = 7))

    two[spike.idx] <- abs(two[spike.idx] + abs(rnorm(1, mean = 0, sd = 0.1)))
    one[spike.idx] <- abs(one[spike.idx] + abs(rnorm(1, mean = 0.02, sd = 0.7)))

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "shock.rds"))
# write.table(data, "~/floa/R/examples/shock.txt", sep="\t")



# ------------------------------------------------------------------------------
# 7. Shift data
# ------------------------------------------------------------------------------
# Inspired by the toy example from:
# https://mjskay.github.io/ggdist/articles/lineribbon.html
set.seed(3)

device <- c()
value <- c()
subjectID <- c()

for (subj.idx in 1:n.subj) {
  value.subj <- c()
  device.subj <- c()
  subjectID.subj <- c()

  t <- seq(0, 100)

  # Subjectwise parameters
  offset.mean <- runif(1, min = -0.5, max = 0.5)

  a.sd <- runif(1, min = 0.05, max = 0.15)
  b.sd <- runif(1, min = 0.0001, max = 0.002)

  for (stride.idx in 1:(n.strides)) {
    a1.1 <- rnorm(1, mean = 3, sd = a.sd)
    a1.2 <- rnorm(1, mean = 3, sd = a.sd)
    a2.1 <- 0.08
    a2.2 <- 0.08
    b1.1 <- rnorm(1, mean = 0.06, sd = b.sd)
    b1.2 <- rnorm(1, mean = 0.06, sd = b.sd)
    b2.1 <- rnorm(1, mean = 0.58, sd = b.sd)
    b2.2 <- rnorm(1, mean = 0.58, sd = b.sd)
    c <- 2

    x.offset <- rnorm(1, offset.mean, 0.05)

    sine.1 <- a1.1 * sin(b1.1 * t) ^ (c + 3)
    sine.2 <- a2.1 * sin(b2.1 * t)
    sine.3 <- a1.2 * sin(b1.2 * t + x.offset) ^ (c + 3)
    sine.4 <- a2.2 * sin(b2.1 * t)

    one <- x.offset + sine.1 + sine.2
    two <- x.offset + sine.3 + sine.4

    value.subj <- c(value.subj, c(two, one))

    device.two <- rep("TWO", n.frames)
    device.one <- rep("ONE", n.frames)
    device.subj <- c(device.subj, c(device.two, device.one))

    subjectID.subj <- c(subjectID.subj, rep(subj.idx, 2 * n.frames))
  }

  # Data is stored for long format output
  device <- c(device, device.subj)
  value <- c(value, value.subj)
  subjectID <- c(subjectID, subjectID.subj)
}

strideID <- rep(1:(n.strides * n.subj), each = n.devices * n.frames)
frame <- rep(0:100, times = n.strides * n.devices * n.subj)

data <- data.frame(device, subjectID, strideID, value, frame)

saveRDS(data, file = paste0("~/floa/R/examples/", "shift.rds"))
# write.table(data, "~/floa/R/examples/shift.txt", sep="\t")






