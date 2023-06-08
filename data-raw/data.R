#' Example data set containing time series (i.e. curve) data
#'
#' The data set is the same as described in Koska et al. (2023)
#' @source  DOI: (\href{https://doi.org/10.1016/j.jbiomech.2023.111506}{Journal of Biomechanics})
#'
#' #' The data set consists of 110 curves and has a hierarchical/nested
#' structure: It contains 10 curves each from 11  subjects. Curves from the same
#' subject all have the same capital letter. Each column represents a curve. All
#' curves are of equal length (n = 101 curve points aka rows).
#'
#' @docType data
#' A subset of data from the World Health Organization Global Tuberculosis
#' Report ...
#'
#' @format ## `data`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{country}{Country name}
#'   \item{iso2, iso3}{2 & 3 letter ISO country codes}
#'   \item{year}{Year}
#'   ...
#' }

## Code to prepare the `data` dataset
data <- readRDS("~/floa/R/examples/non_gaussian.rds")

data.sub <- subset(data, device == "TWO", select = c(subjectID, strideID, value, frame))

n.time <- 101
data.num.long <- data.sub$value
data.num.wide <- matrix(data.num.long, ncol = length(data.num.long) / n.time)
data <- data.num.wide |> data.frame()

cluster.idx <- as.character(rep(LETTERS[1:11], each = 10))
colnames(data) <- cluster.idx

usethis::use_data(data, overwrite = TRUE)

# Convenience function from package 'usethis' that creates 'data_raw/' and adds
# the folder to .Rbuildignore
# usethis::use_data_raw("data")
