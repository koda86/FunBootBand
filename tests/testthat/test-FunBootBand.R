library(testthat)
library(FunBootBand)
# devtools::load_all()

# 1. Test the function behavior with valid inputs
test_that("band() returns something", {
  tmp <- band(data,
              type = "prediction",
              alpha = 0.05,
              iid = FALSE,
              k.coef = 50,
              B = 5)
  expect_type(tmp, "double")
})

# 2. Does the code return the expected value?
load("~/FunBootBand/data/curvesample.RData")

# iid case
dat <- data[-1, ]

test_that("band() returns something", {
  tmp <- band(dat,
              type = "prediction",
              alpha = 0.05,
              iid = TRUE,
              k.coef = 50,
              B = 5)
  expect_type(tmp, "double")
})

# Assume iid when data is not iid or
# ...

# Check for NAs
# Currently, NAs are not allowed.

# Wrong data input format (iid)

# 3. Test the function behavior with invalid inputs
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

# 4. Data Structure Test
# expect_length(): Checks the length of an object.

# Equality Tests:
#   expect_equal(): Checks if two objects are equal.
# expect_identical(): Checks if two objects are identical (same object in memory).
# expect_equivalent(): Checks if two objects are equivalent, accounting for differences in attributes or class.
#
# Comparison Tests:
#   expect_true(): Checks if an expression is TRUE.
# expect_false(): Checks if an expression is FALSE.
# expect_less_than(): Checks if a value is less than another value.
# expect_greater_than(): Checks if a value is greater than another value.
# expect_match(): Checks if a character string matches a regular expression.
#
# Error Handling:
#   expect_error(): Checks if an expression throws an error.
# expect_warning(): Checks if an expression generates a warning.
# expect_message(): Checks if an expression generates a specific message.
#
# Data Structure Tests:
#   expect_length(): Checks the length of an object.
# expect_is(): Checks if an object is of a specific class.
# expect_inherits(): Checks if an object inherits from a specific class.
# expect_named(): Checks if the elements of an object have specific names.
#
# Test Context:
#   skip_if(): Skips a test if a certain condition is met.
# skip(): Skips a test unconditionally.
# context(): Defines a new context or grouping for related tests.

# # In reality, our function is more complex and aggregates your input if you have duplicates in your id-time units -- this is why the following two tests were essential for us
# ## Test whether the output contains the right number of rows
# test_that("overview_tab() returns a dataframe with correct number of rows", {
#   output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#   expect_equal(nrow(output_table), length(unique(toydata$ccode)))
# })

# ## Test whether the function works on a data frame that has no duplicates in id-time
# test_that("overview_tab() works on a dataframe that is already in the correct
#           format",
#           {
#             df_com <- data.frame(
#               # Countries
#               ccode  = c(
#                 rep("RWA", 4),
#                 rep("AGO", 8),
#                 rep("BEN", 2),
#                 rep("GBR", 5),
#                 rep("FRA", 3)
#               ),
#               # Time frame
#               year =
#                 c(
#                   seq(1990, 1995),
#                   seq(1990, 1992),
#                   seq(1995, 1999),
#                   seq(1991, 1999, by = 2),
#                   seq(1993, 1999, by = 3)
#                 )
#             )
#             output_table <-
#               overview_tab(dat = df_com, id = ccode, time = year)
#             expect_equal(nrow(output_table), 5)
#           })
