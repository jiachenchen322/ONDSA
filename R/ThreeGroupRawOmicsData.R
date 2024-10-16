#' Three-Group Raw Omics Data
#'
#' This dataset contains synthetic raw data for three groups, stored as `group1_rawdata`, `group2_rawdata`, and `group3_rawdata`.
#' The data is structured for general use in multi-group analyses, especially for omics studies requiring raw measurements.
#' Each group is measured on the same set of proteins (p is the same across groups)
#'
#' @format A `.rda` file containing three data frames:
#' \describe{
#'   \item{group1_rawdata}{Raw data matrix (n1 * p) for the first group.}
#'   \item{group2_rawdata}{Raw data matrix (n2 * p) for the second group.}
#'   \item{group3_rawdata}{Raw data matrix (n3 * p) for the third group.}
#' }
#'
#' @examples
#' # Load the synthetic raw data
#' data(ThreeGroupRawOmicsData)
#'
#' # Inspect the data
#' str(group1_rawdata)
#' str(group2_rawdata)
#' str(group3_rawdata)

