#' Entry-wise Cochran's Q Statistics and Weighted Mean for Across-Group Comparison
#'
#' Internal function to compute Cochran's Q statistics and weighted mean d for comparing precision matrix entries across groups
#' @param arraymean An array of estimates of a specific precision matrix entry across the K groups
#' @param arraysd2 An array of variances for the estimator of a specific precision matrix entry across K groups
#'
#' @keywords internal
#' @noRd
CochranQ <- function(arraymean, arraysd2) {
  weight = 1 / arraysd2
  d = sum(arraymean * weight) / sum(weight)
  Q = sum(weight * (arraymean - d)^2)
  return(c(Q = Q, d = d))
}

#' Test Statistics for Differential Precision Matrix Entries
#'
#' Internal function that transforms the Chi-square statistic into a z-value for assessing differential omics network structures in precision matrices
#'
#' @param K The number of groups
#' @param p The number of variables (nodes) in each network
#' @param Omega A list of precision matrices estimated for each group
#' @param N A numeric vector with the sample sizes for each group
#'
#' @keywords internal
#' @noRd
Diff.Chisq.omega<- function(K,p,Omega,N){
  Qmatrix=matrix(0,p,p)
  T1matrix=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      omega_ij=NULL
      sd2omega_ij=NULL
      for(k in 1:K){
        omega_ij[k]=Omega[[k]][i,j]
        sd2omega_ij[k]=(1/N[k])*(Omega[[k]][i,i]*Omega[[k]][j,j]+(Omega[[k]][i,j])^2)
      }
      Qmatrix[i,j]=CochranQ(arraymean=omega_ij,arraysd2=sd2omega_ij)[1]
      T1matrix[i,j]=qnorm(pchisq(Qmatrix[i,j], df = K-1))
    }
  }

  return(list(Qmatrix=Qmatrix,T1matrix=T1matrix))
}

#' Test Statistics for Similar Precision Matrix Entries
#'
#' Internal function to compute the z-value for assessing similar omics network structures in precision matrices across groups
#'
#' @param K The number of groups
#' @param p The number of variables (nodes) in each network
#' @param selected.listT2 A list of non-differential edges selected from the first FDR control step
#' @param Omega A list of precision matrices estimated for each group
#' @param N A numeric vector with the sample sizes for each group
#'
#' @keywords internal
#' @noRd
Diff.Z.omega<- function(K,p,selected.listT2,Omega,N){
  sT2matrix=matrix(0,p,p)
  T2matrix=matrix(0,p,p)

  for(i in 1:p){
    for(j in 1:p){
      omega_ij=NULL
      sd2omega_ij=NULL
      for(k in 1:K){
        omega_ij[k]=Omega[[k]][i,j]
        sd2omega_ij[k]=(1/N[k])*(Omega[[k]][i,i]*Omega[[k]][j,j]+(Omega[[k]][i,j])^2)
      }

      sT2matrix[i,j]=CochranQ(arraymean=omega_ij,arraysd2=sd2omega_ij)[2]*sqrt(sum(1/sd2omega_ij))
      T2matrix[i,j]=qnorm(2*pnorm(abs(sT2matrix[i,j]))-1)
    }

  }
  return(T2matrix)
}

