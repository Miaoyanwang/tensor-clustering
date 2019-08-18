install.packages("tensorsparse_1.0.tar.gz",repos=NULL,type="source")
library("tensorsparse")
set.seed(1)
n=40;p=40;q=40;k=5;r=5;l=5;
error=c(4,6,8,10,12,14,16)

new_zero=new_total=cp_zero=cp_total=tucker_zero=tucker_total=matrix(0,nrow=7,ncol=50)
new_e=cp_e=tucker_e=matrix(0,nrow=7,ncol=50)
#new_t=cp_t=tucker_t=matrix(0,nrow=2,ncol=50)



for(i in 1:7){
    for(iter in 1:50){
        data = get.data(n,p,q,k,r,l,sparse.percent=0.8,error[i])
## block model
##krl = choosekrl_bic(data$x,2:6,2:6,2:6)$estimated_krl
        
        range=sqrt(prod(n*p*q)/(k*r*l))*seq(0,2,by=0.1)*error[i]
        lambda_bar = chooseLambda(data$x,k,r,l,lambda=range,method="L0")
        
        our_result=classify2(data$x,k,r,l,lambda_bar$lambda)
        ev = sparse.evaluate(our_result,data,CER=TRUE)
        
        new_zero[i,iter]=ev$correctzerorate
        new_total[i,iter]=ev$totalincorrectrate
        
        new_e[i,iter]=ev$error
#   new_t[i,iter]=ev$cerTotal
## cp model
        cp_result = cp_kmeans(data$x,k,r,l,multiplicative=k)
        cp_ev = sparse.evaluate(cp_result,data,CER=FALSE,show=TRUE)
        cp_zero[i,iter]=cp_ev$correctzerorate
        cp_total[i,iter]=cp_ev$totalincorrectrate
        
        cp_e[i,iter]=cp_ev$error
#cp_t[i,iter]=cp_ev$cerTotal
## tucker model
        tucker_result=tucker_means(data$x,k,r,l)
        t_ev = sparse.evaluate(tucker_result,data,CER=FALSE)
        tucker_zero[i,iter]=t_ev$correctzerorate
        tucker_total[i,iter]=t_ev$totalincorrectrate
        
        tucker_e[i,iter]=t_ev$error
#  tucker_t[i,iter]=t_ev$cerTotal
    }
}

#######
data=load("Figure2_Sparse.RData")
MSE=c(apply(new_e,1,mean),apply(tucker_e,1,mean),apply(cp_e,1,mean))
sd=c(apply(new_e,1,sd),apply(tucker_e,1,sd),apply(cp_e,1,sd))
err=rep(error,3)

method=c(rep("TvBM",7),rep("Tucker+k-means",7),rep("CP+k-means",7))
data=data.frame(MSE=MSE,err=err,sd=sd,method=method)

p=ggplot(data=data,aes(x=err,y=MSE))+geom_line(aes(color=method))+geom_point(aes(shape=method))+ labs(x='noise', y='Root Mean Sqaured Error (RMSE)')

g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(3,2,1)]
p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.5,
                  position=position_dodge(0.05),color=col[rep(rep(1:3),7)])
pdf("clustering_404040_sparse.pdf",width=6,height=4)
p
dev.off()

