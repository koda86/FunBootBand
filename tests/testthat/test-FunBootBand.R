library(FunBootBand)

# Test if the function runs without error (assuming valid inputs)
expect_no_error(band(data,
                     type = "prediction",
                     alpha = 0.05,
                     iid = FALSE,
                     k.coef = 50,
                     B = 5))

# Does the function output have a valid data type?
test_that("band() has a valid return type", {
  expect_type(object = band(data,
                            type = "prediction",
                            alpha = 0.05,
                            iid = FALSE,
                            k.coef = 50,
                            B = 5),
              type = "double")
})

# Does the function output have the correct length?
expect_length(object = band(data,
                            type = "prediction",
                            alpha = 0.05,
                            iid = FALSE,
                            k.coef = 50,
                            B = 5),
              n = 3*dim(data)[1]) # n = expected length

# Test the function behavior with invalid inputs
test_that("Invalid inputs throw appropriate errors", {
  # Invalid 'type' (a)
  expect_error(band(matrix(1:10), type = "invalid", alpha = 0.05, iid = TRUE, k.coef = 50, B = 5),
               "must be either 'confidence' or 'prediction'.")
  # Invalid 'type' (b)
  expect_error(band(matrix(1:10), type = 3, alpha = 0.05, iid = TRUE, k.coef = 50, B = 5),
               "'type' must be a variable of type 'character'.")
  # Invalid 'alpha' (a)
  expect_error(band(matrix(1:10), type = "confidence", alpha = 1.5, iid = TRUE, k.coef = 50, B = 5),
               "'alpha' must be a numeric value between 0 and 1.")
  # Invalid 'alpha' (b)
  expect_error(band(matrix(1:10), type = "confidence", alpha = "0.05", iid = TRUE, k.coef = 50, B = 5),
               "'alpha' must be a numeric value between 0 and 1.")
})

# Test if the header is specified correctly
data.without.names <- data.frame(matrix(unlist(data), nrow = ncol(data)))
test_that("Invalid header", {
  expect_error(band(data.without.names, type = "prediction", alpha = 0.05, iid = FALSE, k.coef = 50, B = 5),
               "The header does not\n
        indicate a nested structure even though 'iid' is set to 'FALSE'.")
})

# Test handling of NA's
data.na <- data
data.na[1,1] <- NA
expect_error(band(data.na,
                     type = "prediction",
                     alpha = 0.05,
                     iid = FALSE,
                     k.coef = 50,
                     B = 5))

# Further 'expectations': https://testthat.r-lib.org/reference/
# expect_named() # Does code return a vector with (given) names?
# expect_true() expect_false() Does code return TRUE or FALSE?


# # Weiterer Test
# load("~/FunBootBand/data/curvesample.RData")
#
# # iid case
# dat <- data[-1, ]
#
# test_that("band() returns something", {
#   tmp <- band(dat,
#               type = "prediction",
#               alpha = 0.05,
#               iid = TRUE,
#               k.coef = 50,
#               B = 5)
#   expect_type(tmp, "double")
# })

# Assume iid when data is not iid or
# ...


# Wrong data input format (iid)





# 4. Data Structure Test

