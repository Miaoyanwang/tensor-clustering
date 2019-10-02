source('functions_synthesis_all_miaoyan.R')
library(ggplot2)

### no covariates
seed=24
d_pre=c(20,40,60,80,100)
p_pre=round(0.5*d_pre)
c_pre=c(3,5,10)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(p_pre,p_pre,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

table=conv_rate(seed,signal=1,Nsim=2,cons="non",lambda = 1,alpha=10,solver ="CG",c_range,dist="normal",dup=1,d_range,p_range,match_dp=TRUE)

save(table,file="table.RData")

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
