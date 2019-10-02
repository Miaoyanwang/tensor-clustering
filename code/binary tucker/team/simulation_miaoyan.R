.libPaths( c( .libPaths(), "/workspace/miaoyan/x86_64-pc-linux-gnu-library/3.4") )
source('functions_synthesis_all_miaoyan.R')

##BIC=sele_rank(tsr,X_covar1 , X_covar2 ,X_covar3 ,rank = 2:4, Nsim,cons = 'non',dist)

data=load("HCP.RData")
tsr=tensor
table(attr[,5])
#####
##22-25 26-30 31+
##35    58    43 

X=attr[,4:5]
levels(X[,2])=c("22-25","26-30","31+","31+") ## three groups

contrasts(X[,1]) <- contr.sum 
contrasts(X[,2]) <- contr.sum 

X_covar3=model.matrix(~as.factor(X[,1])+as.factor(X[,2])) ## baseline female and age 22-25

### intercept, F 1, M,-1; 22-25 [1,0]; 26-30 [0 1]
X_covar1=X_covar2=NULL

##BIC=sele_rank(tsr,NULL, NULL ,X_covar3 ,rank1 = 4:11,rank2 = 4:11,rank3 = 2:4, Nsim,cons = 'non',dist)

core_shape=c(10,10,4)
Nsim=10
cons="non"
lambda=0.1
alpha=10
solver="GC"
dist="binary"

result=update_binary_all(tsr,NULL,NULL,X_covar3,core_shape,Nsim,cons,lambda,alpha,solver,dist)

save(result,file="output.RData")
