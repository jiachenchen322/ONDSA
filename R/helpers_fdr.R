#' FDR Control for Differential Structures (Step 1)
#'
#' Internal function to control the FDR for identifying differential omics network structures across multiple precision matrices
#'
#' @param T1matrix A matrix of test statistics for differential structures
#' @param p The number of variables (nodes) in each network
#' @param alpha The FDR control level
#'
#' @keywords internal
#' @noRd
FDR_control_t1<-function(T1matrix,p,alpha){
  II = diag(p)
  w.upper = which(upper.tri(II))
  z=seq(-floor(2*sqrt(log(p))+1),floor(2*sqrt(log(p))+1),by=0.01)
  s=length(z)
  rr=1
  k=1
  while(rr>alpha & k<s){
    P0=(2*pnorm(1)-1)
    Q0=dnorm(1)*sqrt(2)
    P0hat=2*sum(abs(T1matrix[w.upper])<=1)/(p*(p-1))
    A=(P0-P0hat)/Q0
    At=1+abs(A)*abs(z[k])*dnorm(z[k])/(sqrt(2)*(1-pnorm(z[k])))
    rr=At*((p^2-p)/2)*(1-pnorm(z[k]))/max(sum(T1matrix[w.upper]>=z[k]),1)
    k=k+1
  }
  k=k-1
  selected.list=which(T1matrix[w.upper]>=z[k])
  selected.listT2=which(T1matrix[w.upper]<z[k])
  return(list(selected.list,selected.listT2))
}

#' FDR Control for Similar Structures (Step 2)
#'
#' Internal function to control the FDR for identifying similar omics network structures across multiple precision matrices
#'
#' @param T2matrix A matrix of test statistics for similar structures
#' @param p The number of variables (nodes) in each network
#' @param selected.listT2 A list of non-differential edges selected from the first FDR control step
#' @param alpha The FDR control level
#'
#' @keywords internal
#' @noRd
FDR_control_t2<-function(T2matrix,p,selected.listT2,alpha){
  w.upper = which(upper.tri(T2matrix))
  D=length(selected.listT2)
  z=seq(-floor(2*sqrt(log(p))+1),floor(2*sqrt(log(p))+1),by=0.01)
  s=length(z)
  rr=1
  k=1
  while(rr>alpha & k<s){
    P0=(2*pnorm(1)-1)
    Q0=dnorm(1)*sqrt(2)
    P0hat=sum(abs(T2matrix[w.upper[selected.listT2]])<=1)/D
    A=(P0-P0hat)/Q0
    At=1+abs(A)*abs(z[k])*dnorm(z[k])/(sqrt(2)*(1-pnorm(z[k])))
    rr=At*(D)*(1-pnorm(z[k]))/max(sum(T2matrix[w.upper[selected.listT2]]>=z[k]),1)
    k=k+1
  }
  k=k-1
  selected.list2=which(T2matrix[w.upper[selected.listT2]]>=z[k])
  return(selected.list2)
}
