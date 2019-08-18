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

mse4 = function(bires, data){
  npqs = dim(bires$judgeX)
  n = npqs[1]; p = npqs[2]; q = npqs[3]; s = npqs[4]
  return(sum((bires$judgeX-data$truthX)^2)/n/p/q/s)
}

mse.re1 = c()

p1 = c(); q1 = c()
range.k = c(4,4,4,8,8,14,19); range.r = c(4,4,8,8,8,14,19); range.l = c(4,8,8,8,12,14,19)
for (n in seq(20,70,by=5)){
      cat("Round", n, "\n")
  for (i in 1:7){
    set.seed(i)
    k = range.k[i]; r = range.r[i]; l = range.l[i]
    p = floor(n*log(k)/log(r)); q = floor(n*log(k)/log(l))
    p1 = c(p1, p); q1 = c(q1, q)
    mse_in = c()
    for (j in 1:50){
      data = get.data(n,p,q,k,r,l)
      bires = classify2(data$x, k, r, l)
      mse_in = c(mse_in, mse(bires, data))
    }
    mse.re1 = c(mse.re1, mean(mse_in))
  }
}

p2 = c(); q2 = c(); s2 = c()
mse.re2 = c()
for (n in seq(20,55,by=5)){
  cat("n=",n,":\n")
  set.seed(n)
  k = 4; r = 4; l = 4; m = 4
  p = floor(n*log(k)/log(r)); q = floor(n*log(k)/log(l)); s = floor(n*log(k)/log(m))
  p2 = c(p2, p); q2 = c(q2, q); s2 = c(s2, s)
  mse_in = c()
  for (j in 1:50){
    cat("Round", j, "\n")
    data = get.data4(n,p,q,s,k,r,l,m)
    bires = classify4(data$x, k, r, l, m)
    mse_in = c(mse_in, mse(bires, data))
  }
  mse.re2 = c(mse.re2, mean(mse_in))
}


p3 = c(); q3 = c(); s3 = c()
range.k = c(4,4,4,4,8); range.r = c(4,4,4,8,8); range.l = c(4,4,8,8,8);range.m=c(4,8,8,8,8)
mse.re3 = c()
for (n in seq(15,25,by=2)){
  cat("n=",n,":\n")
  for(i in 1:5){
  set.seed(i)
  k = range.k[i]; r = range.r[i]; l = range.l[i];m=range.m[i]
  p = floor(n*log(k)/log(r)); q = floor(n*log(k)/log(l)); s = floor(n*log(k)/log(m))
  p3 = c(p3, p); q3 = c(q3, q); s3 = c(s3, s)
  mse_in = c()
  for (j in 1:50){
    cat("Round", j, "\n")
    data = get.data4(n,p,q,s,k,r,l,m)
    bires = classify4(data$x, k, r, l, m)
    mse_in = c(mse_in, mse(bires, data))
  }
  mse.re3 = c(mse.re3, mean(mse_in))
} 
}

figure1 = data.frame(npq = rep(seq(20,70,by=5),each=7), d1d2d3 = as.factor(rep(c('(4,4,4)','(4,4,8)','(4,8,8)',"(8,8,8)","(8,8,12)","(14,14,14)","(19,19,19)"), times=11)), 
                     sqrtmse = sqrt(mse.re1), rescalednpq = sqrt(p1*q1/log(rep(c(4,4,4,8,8,14,19),times=11))))

order4 = data.frame(npqs = seq(20,55,by=5), sqrtmse = sqrt(mse.re2), 
                    rescalednpqs = sqrt(p2*q2*s2/log(rep(4,8))))

rs.order4 = data.frame(npqs = seq(15,25,by=2), sqrtmse = sqrt(mse.re3), 
                       rescalednpqs = sqrt(p3*q3*s3/log(rep(4,6))))

ggplot(data=order4, aes(x=npqs,y=sqrtmse))+geom_line()

pdf("non-scale.pdf",width=4,height=4)
ggplot(data=figure1, aes(x=npq,y=sqrtmse))+geom_line(aes(color=d1d2d3))+
scale_shape_manual(values=seq(0,16))+
geom_point(aes(shape=d1d2d3))+
labs(x=expression(Dimension~"in"~the~first~mode~d[1]), y="Root mean squared error (RMSE)", 
color="Number of clusters", shape="Number of clusters")+
geom_point(data=order4,aes(x=npqs,y=sqrtmse,shape="(4,4,4,4)"))+
geom_line(data=order4,aes(x=npqs,y=sqrtmse,color="(4,4,4,4)"))
dev.off()

pdf("rescale.pdf",width=4,height=4)
ggplot(data=figure1, aes(x=rescalednpq,y=sqrtmse))+geom_line(aes(color=d1d2d3))+
  scale_shape_manual(values=seq(0,15))+
  geom_point(aes(shape=d1d2d3))+
  labs(x='Rescaled sample size N', y="Root mean squared error (RMSE)", 
       color="Number of clusters", shape="Number of clusters")+
  geom_point(data=rs.order4,aes(x=rescalednpqs,y=sqrtmse,shape="(4,4,4,4)"))+
  geom_line(data=rs.order4,aes(x=rescalednpqs,y=sqrtmse,color="(4,4,4,4)"))
dev.off()




rs.order4 = data.frame(npqs = seq(15,25,by=2), d1d2d3 = as.factor(rep(c('(4,4,4,4)','(4,4,4,8)','(4,4,8,8)','(4,8,8,8)','(8,8,8,8)'),times=6)),sqrtmse = sqrt(mse.re3), 
rescalednpq = sqrt(p3*q3*s3/log(rep(c(4,4,4,4,8),6))))





pdf("rescale_order4.pdf",width=5,height=3)
ggplot(data=rs.order4, aes(x=rescalednpq,y=sqrtmse))+geom_line(aes(color=d1d2d3))+
scale_shape_manual(values=seq(0,15))+
geom_point(aes(shape=d1d2d3))+
labs(x='Rescaled sample size N', y="Root mean squared error (RMSE)", 
color="Number of clusters", shape="Number of clusters")
dev.off()


