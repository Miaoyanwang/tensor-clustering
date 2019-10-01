source('functions_synthesis_all_miaoyan.R')


seed=24
cons="no"
whole_shape=c(20,20,20)
core_shappe=c(3,3,3)
p=c(0,0,4)
dist="binary"
dup=10
signal=10



## one-side covariates 
data = gene_data_all(seed, whole_shape, core_shape,p,dist, dup, signal)


result=update_binary_all(tsr,NULL,NULL,X_covar3,core_shape,Nsim,cons,lambda,alpha,solver,dist)
