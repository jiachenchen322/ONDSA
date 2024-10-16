#' ONDSA: Omics Network Differential and Similarity Analysis
#'
#' This function performs a two-step differential and similarity analysis on multiple omics networks based on GGMs.
#' It identifies both differential omics network structures (edges that differ across groups) and similar structures (edges that are similar across groups).
#' The procedure conducts rigorous statistical tests with FDR control at both steps.
#' ONDSA requires that each group is measured on the same set of variables, ensuring that the number of variables p remains the same across all groups
#'
#' @param Omega A list of estimated precision matrices, each representing the omics network for a specific group
#' @param p An integer, the number of variables (nodes) in each network
#' @param K An integer, the number of groups (or networks)
#' @param N A numeric vector containing the sample sizes for each group
#' @param alpha A numeric value specifying the FDR control level for differential and similar structures at both steps
#' @param varname A character vector of names for the variables (nodes) in the networks
#'
#' @return A list containing two data frames, differential omics network structures and similar omics network structures:
#' \describe{
#'   \item{differential_structures}{A data frame of edges that differ significantly across the groups}
#'   \item{similar_structures}{A data frame of edges that are similar across the groups}
#' }
#'
#' @details
#' The ONDSA framework works in two main steps:
#' \itemize{
#'   \item Step 1: Identification of differential structures based on a Chi-square statistics with FDR control (based on the test statistics in T1matrix).
#'   \item Step 2: Identification of similar structures using a z-value with FDR control (based on the test statistics in T2matrix).
#' }
#' The function first computes test statistics for differential structures in the precision matrices using \code{Diff.Chisq.omega},
#' followed by FDR control with \code{FDR_control_t1}. In the second step, it identifies edges that are similar across groups using
#' \code{Diff.Z.omega} and applies FDR control with \code{FDR_control_t2}.
#'
#' @importFrom FastGGM FastGGM_Parallel
#' @import Rcpp
#' @import RcppParallel
#'
#' @examples
#' # Load the raw data for the three groups
#' data(ThreeGroupRawOmicsData)
#'
#' # Standardize the raw data for each group
#' Data_1Std <- apply(group1_rawdata, 2, standardize)
#' Data_2Std <- apply(group2_rawdata, 2, standardize)
#' Data_3Std <- apply(group3_rawdata, 2, standardize)
#'
#' # Set parameters
#' alpha <- 0.05
#' p <- 500
#' K <- 3
#' N <- c(nrow(group1_rawdata), nrow(group2_rawdata), nrow(group3_rawdata))
#' varname <- colnames(group1_rawdata)
#'
#' # Estimate group-specific precision matrices using FastGGM
#' n_total <- sum(N)
#' lambda_value <- sqrt(2 * log(p / sqrt(n_total)) / n_total)
#' fastggm_1 <- FastGGM_Parallel(Data_1Std, lambda=lambda_value)
#' fastggm_2 <- FastGGM_Parallel(Data_2Std, lambda=lambda_value)
#' fastggm_3 <- FastGGM_Parallel(Data_3Std, lambda=lambda_value)
#'
#' # Create list of precision matrices
#' Omega <- list(fastggm_1$precision, fastggm_2$precision, fastggm_3$precision)
#'
#' # Run ONDSA with the estimated precision matrices
#' result <- ONDSA(Omega, p, K, N, alpha, varname)
#' result$differential_structures
#' result$similar_structures
#'
#' @export
ONDSA <- function(Omega, p, K, N, alpha, varname) {

  # Step 1: Compute the test statistics for differential structures
  Diff.Chisq = Diff.Chisq.omega(K, p, Omega, N)

  # Extract the test statistics matrix from the output of Diff.Chisq.omega
  T1matrix = Diff.Chisq[[2]]

  # Apply FDR control to select differential edges based on the test statistics in T1matrix
  selected.list = FDR_control_t1(T1matrix, p, alpha)[[1]]

  # Get indices of edges in the upper triangle of the precision matrix
  upperindex = which(upper.tri(Omega[[1]]))

  # Obtain the indices of  non-differential edges
  selected.listT2 = FDR_control_t1(T1matrix, p, alpha)[[2]]

  # Create an empty matrix to mark differential edges
  show = matrix(0, p, p)

  # Mark selected differential edges in the matrix
  show[upperindex[selected.list]] = 1

  # Set column and row names to variable names for better interpretation
  colnames(show) = varname
  rownames(show) = varname

  # Extract the indices of non-zero entries, which indicate the differential edges
  mat_nonzero <- which(show != 0, arr.ind = TRUE, useNames = TRUE)

  # Convert matrix indices into a data frame for easier handling
  mat_nonzero = as.data.frame(mat_nonzero)

  # Add column names to indicate which nodes are represented by each index
  mat_nonzero$colname = varname[mat_nonzero[, 2]]

  # Create a data frame to store the differential edges with both node indices and names
  differential_structures = data.frame(
    node1index = mat_nonzero$row,
    node2index = mat_nonzero$col,
    node1name = varname[mat_nonzero[, 1]],
    node2name = mat_nonzero$colname
  )

  # Step 2: Compute z-values to identify similar structures across groups
  T2matrix = Diff.Z.omega(K, p, selected.listT2, Omega, N)

  # Apply FDR control to select similar edges across all groups
  commonedge = FDR_control_t2(T2matrix, p, selected.listT2, alpha)

  # Reinitialize the empty matrix to mark similar edges
  show = matrix(0, p, p)

  # Mark selected similar edges in the matrix
  show[upperindex[selected.listT2[commonedge]]] = 1

  # Set column and row names to variable names for easier interpretation
  colnames(show) = varname
  rownames(show) = varname

  # Extract the indices of non-zero entries, which indicate the similar edges
  mat_nonzero <- which(show != 0, arr.ind = TRUE, useNames = TRUE)

  # Convert matrix indices into a data frame for easier handling
  mat_nonzero = as.data.frame(mat_nonzero)

  # Add column names to indicate which nodes are represented by each index
  mat_nonzero$colname = varname[mat_nonzero[, 2]]

  # Create a data frame to store the similar edges with both node indices and names
  similar_structures = data.frame(
    node1index = mat_nonzero$row,
    node2index = mat_nonzero$col,
    node1name = varname[mat_nonzero[, 1]],
    node2name = mat_nonzero$colname
  )

  # Step 4: Return the results containing both differential and similar structures
  res = list(
    differential_structures = differential_structures,
    similar_structures = similar_structures
  )

  return(res)
}
