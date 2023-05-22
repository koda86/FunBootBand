context("check-output")  # Call "test-check_output.R"
library(testthat)
library(FunBootBand)
# devtools::load_all()

# Test whether ...
test_that("band() returns something", {
  tmp <- band(data, type = "prediction", B = 5, alpha = 0.05)
  expect_type(tmp, "double")
})
tmp <- band(data, type = "prediction", B = 5, alpha = 0.05)
plot(tmp[1,],
     xlim = c(0, 101),
     ylim = range(tmp),
     type = "l")
lines(tmp[2, ])
lines(tmp[3, ])

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
