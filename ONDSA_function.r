standardize = function(x)
{
  return((x-mean(x))/sd(x))
}

CochranQ<-function(arraymean,arraysd2){
  weight=1/arraysd2
  d=sum(arraymean*weight)/sum(weight)
  Q=sum(weight*(arraymean-d)^2)
  return(c(Q=Q,d=d))
}

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


###FDR control procedure Step 1: differential structures
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


###FDR control procedure Step 2: similar structures
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


#Omega is a list storing group-specific precision matrices
#p is number of variables, all three subgroups should have the same p
#K is number of subgroups
#N is an array storing sample size for each subgroup
#alpha: a pre-specified FDR control level for both steps
#varname is an array of names for variables
ONDSA<-function(Omega,p,K,N,alpha,varname){
  
  Diff.Chisq=Diff.Chisq.omega(K,p,Omega,N)
  T1matrix=Diff.Chisq[[2]]
  selected.list=FDR_control_t1(T1matrix,p,alpha)[[1]]
  selected.listT2=FDR_control_t1(T1matrix,p,alpha)[[2]]
  
  upperindex=which(upper.tri(Omega[[1]]))
  upperindex[selected.list]
  
  show=matrix(0,p,p)
  show[upperindex[selected.list]]=1
  colnames(show)=varname
  rownames(show)=varname
  mat_nonzero <- which(show != 0, arr.ind = T, useNames = TRUE)
  mat_nonzero=as.data.frame(mat_nonzero)
  mat_nonzero$colname=varname[mat_nonzero[,2]]
  differential_structures=data.frame(node1index=mat_nonzero$row,node2index=mat_nonzero$col,node1name=varname[mat_nonzero[,1]],node2name=mat_nonzero$colname)
  
  
  T2matrix=Diff.Z.omega(K,p,selected.listT2,Omega,N)
  commonedge=FDR_control_t2(T2matrix,p,selected.listT2,alpha)
  show=matrix(0,p,p)
  show[upperindex[selected.listT2[commonedge]]]=1
  colnames(show)=varname
  rownames(show)=varname
  mat_nonzero <- which(show != 0, arr.ind = T, useNames = TRUE)
  mat_nonzero=as.data.frame(mat_nonzero)
  mat_nonzero$colname=varname[mat_nonzero[,2]]
  similar_structures=data.frame(node1index=mat_nonzero$row,node2index=mat_nonzero$col,node1name=varname[mat_nonzero[,1]],node2name=mat_nonzero$colname)
  res=list(differential_structures=differential_structures,similar_structures=similar_structures)
  return(res)
}

