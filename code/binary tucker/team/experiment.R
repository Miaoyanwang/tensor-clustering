source('functions_synthesis_all_miaoyan.R')
library(ggplot2)

seed=24

dist="binary"
dup=2
signal=10
d=c(20,30,40,50,60)


## one-side covariates  Figure 1
data = gene_data_all(seed, whole_shape, core_shape,p,dist, dup, signal)
X_covar1=NULL
X_covar2=NULL
Nsim=20
cons="no"
lambda=0.1
alpha=sigma
solver="GC"
core_shape=c(3,3,3)


p=c(0,0,4)
error_matrix=matrix(nrow=5,ncol=dup)
KL_matrix=matrix(nrow=5,ncol=dup)

for(i in 1:5){
whole_shape=rep(d[i],3)
data = gene_data_all(seed, whole_shape, core_shape,p,dist, dup, signal)
X_covar3=data$X_covar3
error=error_KL=NULL
for(j in 1:dup){
tsr=data$tsr[[j]]
result=update_binary_all(tsr,X_covar1,X_covar2,X_covar3,core_shape,Nsim,cons,lambda,alpha,solver,dist)
error=c(error,mean((result$U-data$U)^2))
error_KL=c(error_KL,KL(result$U,data$U,dist,sigma_est=result$sigma))
}
error_matrix[i,]=error
KL_matrix[i,]=error_KL
}

ggplot(re, aes(x = d, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 1.5)  +
geom_point( size = 3) +
theme(axis.text=element_text(size=12),
axis.title=element_text(size=14,face="bold"))

