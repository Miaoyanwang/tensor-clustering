rm(list = ls())
library(tensorregress)
seed=24
d=c(40,60,80)
dist="normal"

core_shape=rbind(c(3,3,3),c(4,4,6))#,c(6,8,8))

#############

#############
dup=1
final=array(0,dim=c(2,3,2))
for(s in 1:2){
  for(i in 1:3){
    whole_shape = c(d[i],20,20)
    data=sim_data(seed, whole_shape = whole_shape, core_shape=core_shape[s,],p=c(12,0,0),dist=dist, dup=dup, signal=4)
    #rank = expand.grid(rank1,rank2,rank3)
    #rank=rank[which(apply(rank,1,function(x){max(x)<= sqrt(prod(x))})==1),]
    #rank_range=rank[which(apply(rank,1,function(x){sum(x<=p)==3})),]
    table=NULL
    
    res = tensor_regress(data$tsr[[1]],X_covar1 = data$X_covar1,core_shape=core_shape[s,],Nsim=10, cons = 'non', dist = dist)
    err_res = sum((res$C_ts - data$C_ts)^2)
    
    holrr = HOLRR(X = data$X_covar1, data$tsr[[1]], core_shape = core_shape[s,])
    err_holrr = sum((holrr$C_ts - data$C_ts)^2)
    
    table=rbind(table,res$rank)
    
    
    final[s,i,1] = err_res; final[s,i,2] = err_holrr 
  }
}

rank333 = data.frame(d = rep(c(40,60,80),2), alg = c(rep('glr',3),rep('holrr',3)), 
                     MSE = c(final[1,,1],final[1,,2]))
rank333 

rank446 = data.frame(d = rep(c(40,60,80),2), alg = c(rep('glr',3),rep('holrr',3)), 
                     MSE = c(final[2,,1],final[2,,2]))
rank446 

 

library(ggplot2)

ggplot(rank333, aes(x = d, y = MSE)) + geom_line(aes(color = as.factor(alg)),size = 1.5)  +
  geom_point( size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



