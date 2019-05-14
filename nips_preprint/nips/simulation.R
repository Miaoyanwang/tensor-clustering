install.packages("tensorsparse_1.0.tar.gz",repos=NULL,type="source")
library("tensorsparse")
set.seed(1)
n=40;p=40;q=40;k=3;r=4;l=5;

error=c(8,10,12,14,16)


data = get.data(n,p,q,k,r,l,error[i])
krl = choosekrl_bic(data$x,2:6,2:6,2:6)$estimated_krl
our_result=classify2(data$x,krl[1],krl[2],krl[3])
ev = sparse.evaluate(our_result,data,CER=TRUE)
new=c(new,(ev$cerC+ev$cerD+ev$cerD)/3)

cp_result = cp_kmeans(data$x,k,r,l,max.s=10)
cp_ev = sparse.evaluate(cp_result,data,CER=TRUE,show=TRUE)
cp=c(cp,(cp_ev$cerC+cp_ev$cerD+cp_ev$cerD)/3)

tucker_result=tucker_means(data$x,k,r,l)

