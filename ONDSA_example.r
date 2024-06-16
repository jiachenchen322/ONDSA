library(huge)
library(BDgraph)
library(mvtnorm)
library(FastGGM)
source("ONDSA_function.r")
##load simulated data of three subgroups
load("simdata.rda")
##this example have three subgroups, with 500 variables, and sample size to be 300/300/300
#standardize data
alpha=0.05
p=500
K=3
Data_1Std = standardize(simdata1)
Data_2Std = standardize(simdata2)
Data_3Std = standardize(simdata3)
n=dim(simdata1)[1]+dim(simdata1)[2]+dim(simdata1)[3]
############################################################
########FastGGM for group-specific precision matrices#######
############################################################
fastggm_1=FastGGM_Parallel(Data_1Std,lambda=sqrt(2*log(p/sqrt(n))/n))
fastggm_2=FastGGM_Parallel(Data_2Std,lambda=sqrt(2*log(p/sqrt(n))/n))
fastggm_3=FastGGM_Parallel(Data_3Std,lambda=sqrt(2*log(p/sqrt(n))/n))

Omega=NULL
Omega[[1]]=fastggm_1$precision
Omega[[2]]=fastggm_2$precision
Omega[[3]]=fastggm_3$precision
N=c(dim(simdata1)[1],dim(simdata2)[1],dim(simdata3)[1])
varname=colnames(simdata1)

#Omega is a list storing group-specific precision matrices
#p is number of variables, all three subgroups should have the same p
#K is number of subgroups
#N is an array storing sample size for each subgroup
#alpha: a pre-specified FDR control level for both steps
#varname is an array of names for variables
res=ONDSA(Omega,p,K,N,alpha,varname)
##res will return differential_structures and similar_structures, respectively
