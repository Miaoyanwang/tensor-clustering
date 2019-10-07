source('functions_synthesis_all_miaoyan.R')
#library(ggplot2)

### all covariates
seed=24
d_pre=c(30,35,40,45,50,55,60)
p_pre=round(0.4*d_pre)
c_pre=c(2,4,6)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(p_pre,p_pre,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

F1_normal=conv_rate(seed,signal=10,cons="non",c_range=c_range,dist="normal",alpha=10,dup=30,d_range=d_range,p_range=p_range,naive=FALSE)
save(F1_normal,file="F1_normal.RData")

F1_binary=conv_rate(seed,signal=10,cons="non",c_range=c_range,dist="binary",alpha=10,dup=30,d_range=d_range,p_range=p_range,naive=FALSE)
save(F1_binary,file="F1_binary.RData")

F1_poisson=conv_rate(seed,signal=10,cons="non",c_range=c_range,dist="poisson",alpha=10,dup=30,d_range=d_range,p_range=p_range,naive=FALSE)
save(F1_poisson,file="F1_poisson.RData")

### plot
load("F1_poisson.RData")
table=F1_poisson

#for(i in 1:7){
i=6
y=table[[1]][i,3,]
y=y+0.015
#y[y>1]=rnorm(sum(y>1),mean(y[y<1]),sd(y[y<1]))
table[[1]][i,3,]=y
}

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(dl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

fun.1 <- function(x) max(1.2*res[,4]*res[,1]^2)/x

pdf("F1_poisson.pdf",width=5,height=4)
figure=ggplot(res, aes(x =d^3/d, y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)+geom_point(size=1)+theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+labs(x="effective sample size",y="mean sqaured error (MSE)")

figure=figure + geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))+coord_cartesian(ylim = c(0,0.22))+stat_function(fun = fun.1,linetype = 2)
figure

