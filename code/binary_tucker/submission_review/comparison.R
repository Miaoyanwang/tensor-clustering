rm(list = ls())
source('methods.r')
library(tensorregress)
seed=25
d=c(40,50,60,70,80,90,100,120)
dist="normal"

core_shape=rbind(c(3,3,3),c(6,6,6))#,c(6,8,8))

#############

#############
dup=1
final=array(0,dim=c(2,8,3))
for(s in 1:2){
  for(i in 1:8){
    whole_shape = c(d[i],20,20)
    data=sim_data(seed, whole_shape = whole_shape, core_shape=core_shape[s,],p=c(12,0,0),dist=dist, dup=dup, signal=4)
    #rank = expand.grid(rank1,rank2,rank3)
    #rank=rank[which(apply(rank,1,function(x){max(x)<= sqrt(prod(x))})==1),]
    #rank_range=rank[which(apply(rank,1,function(x){sum(x<=p)==3})),]
    table=NULL
    
    ## GTR
    res = tensor_regress(data$tsr[[1]],X_covar1 = data$X_covar1,core_shape=core_shape[s,],Nsim=10, cons = 'non', dist = dist)
    err_res = sum((res$C_ts - data$C_ts)^2)
    
    ## HOLRR
    holrr = HOLRR(X = data$X_covar1, data$tsr[[1]], core_shape = core_shape[s,])
    err_holrr = sum((holrr$C_ts - data$C_ts)^2)
    
    ## TPG
    tpg = TPG(X = data$X_covar1, data$tsr[[1]], core_shape[s,1], eta = 0.5, iter = 800, tol = 0.05, seed = 1)
    err_tpg = sum((tpg - data$C_ts)^2)

    
    final[s,i,1] = err_res; final[s,i,2] = err_holrr 
    final[s,i,3] = err_tpg 
  }
}

rank333 = data.frame(d = rep(c(40,50,60,70,80,90,100,120),3), algorithm = c(rep('gtr',length(d)),rep('holrr',length(d)),rep('tpg',length(d))), 
                     MSE = c(final[1,,1],final[1,,2], final[1,,3]))
rank333 

rank666 = data.frame(d = rep(c(40,50,60,70,80,90,100,120),3), algorithm = c(rep('gtr',length(d)),rep('holrr',length(d)),rep('tpg',length(d))), 
                     MSE = c(final[2,,1],final[2,,2], final[2,,3]))
rank666 



library(ggplot2)

ggplot(rank333, aes(x = d, y = MSE)) + geom_line(aes(color = algorithm),size = 1.5)  +
  geom_point( size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  #xlab('Sample Size') + 
  labs(title = 'rank:3 3 3',size = 5) + 
  theme(plot.title = element_text(hjust = 0.5,size = 28))

ggplot(rank666, aes(x = d, y = MSE)) + geom_line(aes(color = algorithm),size = 1.5)  +
  geom_point( size = 3) +
  labs(title = 'rank:6 6 6',size = 5) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5,size = 28))


