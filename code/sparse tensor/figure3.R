rm(list=ls())
require(tensorsparse)
require(ggplot2)

########################################
###   Figure 3    ######################
########################################

mse = function(bires, data){
  npq = dim(bires$judgeX)
  n = npq[1]; p = npq[2]; q = npq[3]
  return(sum((bires$judgeX-data$truthX)^2)/n/p/q)
} 


mse.re1 = c()

p1 = c(); q1 = c()
range.k = c(4,4,8,8); range.r = c(4,8,8,8); range.l = c(4,8,8,12)
for (n in seq(20,70,by=5)){
  for (i in 1:4){
    set.seed(n)
    k = range.k[i]; r = range.r[i]; l = range.l[i]
    p = floor(n*log(k)/log(r)); q = floor(n*log(k)/log(l))
    p1 = c(p1, p); q1 = c(q1, q)
    data = get.data(n,p,q,k,r,l)
    bires = classify2(data$x, k, r, l)
    mse.re1 = c(mse.re1, mse(bires, data))
  }
}

figure1 = data.frame(npq = rep(seq(20,70,by=5),each=4), krl = as.factor(rep(c("d1=4,d2=4,d3=4","d1=4,d2=8,d2=8","d1=8,d2=8,d3=8","d1=8,d2=8,d3=12"), times=11)), 
                     sqrtmse = sqrt(mse.re1), rescalednpq = sqrt(rep(seq(20,70,by=5),each=4)/log(rep(c(4,8,8,12),times=11))))

pdf("non-scale.pdf")
ggplot(data=figure1, aes(x=npq,y=sqrtmse))+geom_line(aes(color=krl))+geom_point(aes(shape=krl))+
  labs(x='sample size n1', y='sqrt(mse)')
dev.off()

pdf("rescale.pdf")
ggplot(data=figure1, aes(x=rescalednpq,y=sqrtmse))+geom_line(aes(color=krl))+geom_point(aes(shape=krl))+
  labs(x='rescaled sample size N', y='sqrt(mse)')
dev.off()

