#' Standardize a Numeric Vector
#'
#' This function standardizes a numeric vector by subtracting the mean and dividing by the standard deviation
#' It can also be applied to standardize each column (variable) of a matrix
#' This function is particularly useful for preprocessing omics data, where standardization ensures that all variables
#' (e.g., genes, proteins, metabolites) are brought to a common scale before performing GGM estimation
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector, where each value is centered by subtracting the mean and scaled by the standard deviation
#'
#' @examples
#' # Standardize a numeric vector
#' x <- c(1, 2, 3, 4, 5)
#' standardize(x)
#'
#' # Standardize each column of a matrix
#' mat <- matrix(c(1:15), nrow=5, ncol=3)
#' apply(mat, 2, standardize)
#'
#' @export

standardize <- function(x) {
  return((x - mean(x)) / sd(x))
}

