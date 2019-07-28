install.packages("tensorsparse_1.0.tar.gz",repos=NULL,type="source")
library("tensorsparse")
set.seed(1)
n=60;p=60;q=60;k=5;r=5;l=5;
error=c(4,8,10,12,16)

new=cp=tucker=matrix(0,nrow=2,ncol=40)
new_e=cp_e=tucker_e=matrix(0,nrow=2,ncol=40)
new_t=cp_t=tucker_t=matrix(0,nrow=2,ncol=40)

for(i in 1:2){
    for(iter in 1:50){
        data = get.data(n,p,q,k,r,l,error[i])
## block model
##krl = choosekrl_bic(data$x,2:6,2:6,2:6)$estimated_krl
our_result=classify2(data$x,k,r,l)
ev = sparse.evaluate(our_result,data,CER=TRUE)
new[i,iter]=(ev$cerC+ev$cerD+ev$cerD)/3
new_e[i,iter]=ev$error
new_t[i,iter]=ev$cerTotal
## cp model
cp_result = cp_kmeans(data$x,k,r,l,multiplicative=k)
cp_ev = sparse.evaluate(cp_result,data,CER=TRUE,show=TRUE)
cp[i,iter]=(cp_ev$cerC+cp_ev$cerD+cp_ev$cerD)/3
cp_e[i,iter]=cp_ev$error
cp_t[i,iter]=cp_ev$cerTotal
## tucker model
tucker_result=tucker_means(data$x,k,r,l)
t_ev = sparse.evaluate(tucker_result,data,CER=TRUE)
tucker[i,iter]=(t_ev$cerC+t_ev$cerD+t_ev$cerD)/3
tucker_e[i,iter]=t_ev$error
tucker_t[i,iter]=t_ev$cerTotal
}
}


MSE=c(apply(new_e,1,mean),apply(tucker_e,1,mean),apply(cp_e,1,mean))
sd=c(apply(new_e,1,sd),apply(tucker_e,1,sd),apply(cp_e,1,sd))
err=rep(error,3)
method=c(rep("TvBM",7),rep("Tucker+k-means",7),rep("CP+k-means",7))
data=data.frame(MSE=MSE,err=err,sd=sd,method=method)

p=ggplot(data=data,aes(x=err,y=MSE))+geom_line(aes(color=method))+geom_point(aes(shape=method))+ labs(x='noise level', y='Root Mean Sqaured Error (RMSE)')

pdf("error_noise_estimation.pdf",width=5,height=4)
p
dev.off()


g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(3,2,1)]
p=p+geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.5,position=position_dodge(0.05),color=col[rep(rep(1:3),7)])


##We then generate using CP models
d=50
error=3
factor=shapes.two.moon(numObjects=d/2,shape1a=0,shape2b=1,shape1rFrom=1,shape1rTo=1,shape2rFrom=1, shape2rTo=1, outputCsv="", outputCsv2="",outputColNames=TRUE, outputRowNames=TRUE)

factor=cbind(c(rnorm(25,-3,0.1),rnorm(25,3,.1)),c(rnorm(25,-3,.1),rnorm(25,3,.1)))
tensor=outer(outer(factor[,1],factor[,1]),factor[,1])+outer(outer(factor[,2],factor[,2]),factor[,2])
#tensor=outer(outer(factor$data[,1],factor$data[,1]),factor$data[,1])+outer(outer(factor$data[,2],factor$data[,2]),factor$data[,2])
obs=tensor+array(rnorm(d^3,0,error),dim=c(d,d,d))

mus=array(0,dim=c(2,2,2))
truthCs=truthDs=truthEs=c(rep(1,(d/2)),rep(2,(d/2)))
for(i in 1:2){
for(j in 1:2){
    for(k in 1:2){
mus[i,j,k]=mean(tensor[truthCs==i,truthDs==j,truthEs==k])
}
}
}
data=list(x=obs,truthX=tensor,truthCs=truthCs,truthDs=truthDs,truthEs=truthEs,binaryX=obs[tensor==0],mus=mus)

###
our_result=classify2(data$x,k,r,l)
ev = sparse.evaluate(our_result,data,CER=TRUE)
new[i,iter]=(ev$cerC+ev$cerD+ev$cerD)/3
new_e[i,iter]=ev$error
new_t[i,iter]=ev$cerTotal
## cp model
cp_result = cp_kmeans(data$x,k,r,l,multiplicative=k)
cp_ev = sparse.evaluate(cp_result,data,CER=TRUE,show=TRUE)
cp[i,iter]=(cp_ev$cerC+cp_ev$cerD+cp_ev$cerD)/3
cp_e[i,iter]=cp_ev$error
cp_t[i,iter]=cp_ev$cerTotal
## tucker model
tucker_result=tucker_means(data$x,k,r,l)
t_ev = sparse.evaluate(tucker_result,data,CER=TRUE)
tucker[i,iter]=(t_ev$cerC+t_ev$cerD+t_ev$cerD)/3
tucker_e[i,iter]=t_ev$error
tucker_t[i,iter]=t_ev$cerTotal
