source('functions_synthesis_all_miaoyan.R')
library(ggplot2)
seed=24


## one-side covariates ## Figure 2
d_pre=c(20,25,30,35,40)
p_pre=round(0.5*d_pre)
c_pre=c(2,4)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(0,0,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

###
c_pre=c(2,5,10)
d_range=cbind(rep(40,4),rep(40,4),c(20,40,60,80))
p_range=cbind(0,0,rep(4,5))
c_range=cbind(c_pre,c_pre,rep(4,3))

p_range=t(matrix(c(0,0,4)))
d_range=t(matrix(c(20,20,50)))
c_range=c(5,5,4)
c_range=rbind(c_range,c(10,10,4),c(15,15,4))

table=conv_rate(seed,signal=5,Nsim=10,cons="non",lambda = 1,alpha=10,solver ="CG",c_range=c_range,dist="poisson",dup=1,d_range=d_range,p_range=p_range,naive=TRUE)
## 
#signal=5;Nsim=10;cons="non";lambda = 1;alpha=10;solver="CG";c_range=c_range;dist="binary";dup=1;d_range=d_range;p_range=p_range;naive=TRUE

p_range=t(matrix(c(0,0,4)))
d_range=t(matrix(c(20,20,50)))
c_range=c(5,5,4)  ##effective rank 2,2,3 
c_range=rbind(c_range,c(10,10,4),c(15,15,4)) ## effective rank 2,3,3; 
## binary; effective rank: (2,2,3); (2,2,3); (2,2,3)
## normal: (2,2,2), (2,2,2), (6,6,4)
## poisson: (9,9,4)

###### rank selection
data=gene_data_all(seed,whole_shape=d_range[1,],core_shape=c_range[3,],p=p_range[1,],dist="poisson",dup=10,signal=5,block=c(TRUE,TRUE,FALSE))
###############

rank = rbind(c(7,7,4),c(8,8,4),c(9,9,4))
rank_range=rank[which(apply(rank,1,function(x){max(x)<= sqrt(prod(x))})==1),]
res=sele_rank(data$tsr[[1]],data$X_covar1,data$X_covar2,data$X_covar3,rank_range=rank_range,dist="poisson")

### normal
rank_choice=rbind(c(5,5,4), c(10,10,4), c(15,15,4))
table=table_naive=matrix(0,nrow=3,ncol=30)
##############################
for(s in 1:3){
data=gene_data_all(seed,whole_shape=d_range[1,],core_shape=c_range[s,],p=p_range[1,],dist="normal",dup=30,signal=5,block=c(TRUE,TRUE,FALSE))
##############################
for(i in 1:30){
result=update_binary_all(tsr=data$tsr[[i]],X_covar1 =data$X_covar1, X_covar2 = data$X_covar2,X_covar3= data$X_covar3, core_shape=rank_choice[s,], Nsim=10, cons, lambda = 0.1, alpha = 1, solver ="CG",dist="normal")
table[s,i]=mean((result$C_ts-data$C_ts)^2)
####
naive_C=massive_glm(data$tsr[[i]],data$X_covar3,dist="normal")
table_naive[s,i]=mean((naive_C-data$C_ts)^2)
}
}
####


### poisson
rank_choice=rbind(c(5,5,4), c(10,10,4), c(15,15,4))
table=table_naive=matrix(0,nrow=3,ncol=30)
##############################
for(s in 1:3){
    data=gene_data_all(seed,whole_shape=d_range[1,],core_shape=c_range[s,],p=p_range[1,],dist="poisson",dup=30,signal=5,block=c(TRUE,TRUE,FALSE))
    ##############################
    for(i in 1:30){
        result=update_binary_all(tsr=data$tsr[[i]],X_covar1 =data$X_covar1, X_covar2 = data$X_covar2,X_covar3= data$X_covar3, core_shape=rank_choice[s,], Nsim=10, cons, lambda = 0.1, alpha = 1, solver ="CG",dist="poisson")
        table[s,i]=mean((result$C_ts-data$C_ts)^2)
        ####
        naive_C=massive_glm(data$tsr[[i]],data$X_covar3,dist="poisson")
        table_naive[s,i]=mean((naive_C-data$C_ts)^2)
    }
}
###
### binary
rank_choice=rbind(c(5,5,4), c(10,10,4), c(15,15,4))
table=table_naive=matrix(0,nrow=3,ncol=30)
##############################
for(s in 1:3){
    data=gene_data_all(seed,whole_shape=d_range[1,],core_shape=c_range[s,],p=p_range[1,],dist="binary",dup=30,signal=5,block=c(TRUE,TRUE,FALSE))
    ##############################
    for(i in 1:30){
        result=update_binary_all(tsr=data$tsr[[i]],X_covar1 =data$X_covar1, X_covar2 = data$X_covar2,X_covar3= data$X_covar3, core_shape=rank_choice[s,], Nsim=10, cons, lambda = 0.1, alpha = 1, solver ="CG",dist="binary")
        table[s,i]=mean((result$C_ts-data$C_ts)^2)
        ####
        naive_C=massive_glm(data$tsr[[i]],data$X_covar3,dist="binary")
        table_naive[s,i]=mean((naive_C-data$C_ts)^2)
    }
}
####

load("comparison_binary.RData")
mean=c(apply(table,1,mean),apply(table_naive,1,mean))
sd=c(apply(table,1,sd),apply(table_naive,1,sd))
data=cbind(round(mean,2),sd,c(1,1,1,2,2,2),c(1,2,3,1,2,3))
data=data.frame(data)
names(data)=c("MSE","sd","method","rank")

pdf("comparison_poisson.pdf",width=5,height=4)
p=ggplot(data=data, aes(x=as.factor(rank), y=MSE, fill=as.factor(method))) +
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=.2,
position=position_dodge(.9))+scale_fill_manual(values=c('#069AA0','#C0C0C0'))+labs(x="number of blocks",y="mean sqaured error (MSE)")+coord_cartesian(ylim = c(0, 1)) 
p
dev.off()


MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(dl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

figure=ggplot(res, aes(x = (d^2*p)/(2*d+p), y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)  +
geom_point(size = 1) +
theme(axis.text=element_text(size=12),
axis.title=element_text(size=14,face="bold"))
figure

figure=figure+geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))
figure


core_shape=c(5,5,4)
whole_shape = c(50,50,50)
p=c(0,0,4)
est=array(0,dim=c(whole_shape[1:2],p[3]))

data=gene_data_all(seed, whole_shape , core_shape ,p,"poisson", dup=1, signal=1)
res=update_binary_all(data$tsr[[1]],X_covar1 = data$X_covar1, X_covar2 = data$X_covar2,X_covar3 = data$X_covar3, core_shape, Nsim=10, cons="no", lambda = 0.1, alpha = 20, solver = "CG","normal")

mean((data$C_ts-res$C_ts)^2)

for(i in 1:whole_shape[1]){
    for(j in 1:whole_shape[2]){
      fit=lm(data$tsr[[1]][i,j,]~-1+data$X_covar3)
      ##fit=glm(data$tsr[[1]][i,j,]~-1+data$X_covar3,family=poisson("log"))
est[i,j,]=coef(fit)
}
}
mean((data$C_ts-est)^2)
