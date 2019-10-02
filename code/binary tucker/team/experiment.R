source('functions_synthesis_all_miaoyan.R')
library(ggplot2)

seed=24

## one-side covariates ## Figure 1
p_pre=c(5,10,15,20)
d_pre=c(20,30,40,50,60)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(0,0,p_pre)

table=conv_rate(seed,signal=10,Nsim=1,cons="no",lambda = 1,alpha=10,solver ="GC",core_shape=c(3,3,3),dist="binary",dup=3,d_range,p_range)

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
re=cbind(rep(d_pre,length(p_pre)),rep(p_pre,rep(length(d_pre),length(p_pre))),c(MSE_matrix),c(sd_matrix))
re=data.frame(re)
colnames(re)=c("d","p","MSE","sd")


p=ggplot(re, aes(x = d, y = MSE)) + geom_line(aes(color = as.factor(p)),size = 1)  +
geom_point(size = 1) +
theme(axis.text=element_text(size=12),
axis.title=element_text(size=14,face="bold"))

p=p+geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.5,position=position_dodge(0.05))


## one-side covariates ## Figure 1
d_pre=c(20,25,30,35,40)
p_pre=round(0.4*d_pre)
c_pre=c(2,4,6)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(0,0,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

table=conv_rate(seed,signal=10,Nsim=1,cons="non",lambda = 1,alpha=10,solver ="CG",c_range,dist="poisson",dup=2,d_range,p_range,match_dp=TRUE)

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(pl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

figure=ggplot(res, aes(x = (d^2*p)/(2*d+p), y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)  +
geom_point(size = 1) +
theme(axis.text=element_text(size=12),
axis.title=element_text(size=14,face="bold"))
figure

figure=figure+geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))
figure



#### all covariates ## Figure 2
d_pre=c(20,30,40,50,60)
p_pre=round(0.2*d_pre)
c_pre=c(3,5,10)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(p_pre,p_pre,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

table=conv_rate(seed,signal=1,Nsim=20,cons="non",lambda = 1,alpha=10,solver ="CG",c_range,dist="normal",dup=30,d_range,p_range,match_dp=TRUE)

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(dl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

figure=ggplot(res, aes(x =p^3/(3*p), y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)  +
geom_point(size = 1) +
theme(axis.text=element_text(size=12),
axis.title=element_text(size=14,face="bold"))
figure=figure+geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))
figure



data=gene_data_all(seed, whole_shape = c(30,30,30), core_shape =c(4,4,4),p=c(0,0,6),"poisson", dup=1, signal=1)

res=update_binary_all (data$tsr[[1]],X_covar1 = data$X_covar1, X_covar2 = data$X_covar2,X_covar3 = data$X_covar3, core_shape=c(4,4,4), Nsim=10, cons="no", lambda = 0.1, alpha = 20, solver = "CG","poisson")
    
plot(data$U,res$U)    
    
