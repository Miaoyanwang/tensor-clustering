#####################################################
# Table 3 - Evaluation of biclustering methods
# The mean matrix is not sparse 
# Simulation for Various clustering methods on multiple scenarios
#####################################################
library(sparseBC)
library(clues)


Do<-function(n,p,k,r,lambda,lambda2,lambda3,iterations){
# Initialize some variables
kcerrow<-c(rep(NA,iterations))
bcerrow<-c(rep(NA,iterations))
Lbcerrow<-c(rep(NA,iterations))
Lbcerrow2<-c(rep(NA,iterations))
Lbcerrow3<-c(rep(NA,iterations))
kcercol<-c(rep(NA,iterations))
bcercol<-c(rep(NA,iterations))
Lbcercol<-c(rep(NA,iterations))
Lbcercol2<-c(rep(NA,iterations))
Lbcercol3<-c(rep(NA,iterations))

SparseLb<-c(rep(NA,iterations))
SparseLb2<-c(rep(NA,iterations))
SparseLb3<-c(rep(NA,iterations))

for(a in 1:iterations){
cat("time",a,fill=TRUE)
set.seed(a)
  mus<-runif(k*r,-2,2)
  mus<-matrix(c(mus),nrow=k,ncol=r,byrow=FALSE)
  truthCs<-sample(1:k,n,rep=TRUE)
  truthDs<-sample(1:r,p,rep=TRUE)
  x<-matrix(rnorm(n*p,mean=0,sd=4),nrow=n,ncol=p)
# This is used to calculate the sparsity error rate
  truthmatrix<-matrix(NA,nrow=n,ncol=p)
  for(i in 1:max(truthCs)){
    for(j in 1:max(truthDs)){
      x[truthCs==i,truthDs==j]<-x[truthCs==i,truthDs==j]+mus[i,j]
      truthmatrix[truthCs==i,truthDs==j]<-mus[i,j]
    } 
  }	
  x<-x-mean(x)
  kmeanres<-kmeans(x,k,nstart=20)
  kmeanres2<-kmeans(t(x),r,nstart=20)
  bires<-sparseBC(x,k,r,0)
  Lbires<-sparseBC(x,k,r,lambda)
  Lbires2<-sparseBC(x,k,r,lambda2)
  Lbires3<-sparseBC(x,k,r,lambda3)

  kcerrow[a]<-1-adjustedRand(truthCs,kmeanres$cluster,randMethod=c("Rand"))
  bcerrow[a]<-1-adjustedRand(truthCs,bires$Cs,randMethod=c("Rand"))
  Lbcerrow[a]<-1-adjustedRand(truthCs,Lbires$Cs,randMethod=c("Rand"))
  Lbcerrow2[a]<-1-adjustedRand(truthCs,Lbires2$Cs,randMethod=c("Rand"))
  Lbcerrow3[a]<-1-adjustedRand(truthCs,Lbires3$Cs,randMethod=c("Rand"))

  kcercol[a]<-1-adjustedRand(truthDs,kmeanres2$cluster,randMethod=c("Rand"))
  bcercol[a]<-1-adjustedRand(truthDs,bires$Ds,randMethod=c("Rand"))
  Lbcercol[a]<-1-adjustedRand(truthDs,Lbires$Ds,randMethod=c("Rand"))
  Lbcercol2[a]<-1-adjustedRand(truthDs,Lbires2$Ds,randMethod=c("Rand"))
  Lbcercol3[a]<-1-adjustedRand(truthDs,Lbires3$Ds,randMethod=c("Rand"))
 
# Record the sparsity of estimated mus

  SparseLb[a]<-sum(abs(Lbires$mus)<1e-6)/(n*p)

  SparseLb2[a]<-sum(abs(Lbires2$mus)<1e-6)/(n*p)
  
  SparseLb3[a]<-sum(abs(Lbires3$mus)<1e-6)/(n*p)

  
 
}
return(list(kmeanscer=c(mean(kcerrow),sd(kcerrow)/sqrt(iterations)),Biclustercer=c(mean(bcerrow),sd(bcerrow)/sqrt(iterations)),Lbicluster=c(mean(Lbcerrow),sd(Lbcerrow)/sqrt(iterations)),Lbicluster2=c(mean(Lbcerrow2),sd(Lbcerrow2)/sqrt(iterations)),Lbicluster3=c(mean(Lbcerrow3),sd(Lbcerrow3)/sqrt(iterations)),kmeanscer=c(mean(kcercol),sd(kcercol)/sqrt(iterations)),Biclustercer=c(mean(bcercol),sd(bcercol)/sqrt(iterations)),Lbicluster=c(mean(Lbcercol),sd(Lbcercol)/sqrt(iterations)),Lbicluster2=c(mean(Lbcercol2),sd(Lbcercol2)/sqrt(iterations)),Lbicluster3=c(mean(Lbcercol3),sd(Lbcercol3)/sqrt(iterations)),Sparsityrate1=c(mean(SparseLb),sd(SparseLb)/sqrt(iterations)),Sparsityrate2=c(mean(SparseLb2),sd(SparseLb2)/sqrt(iterations)),Sparsityrate3=c(mean(SparseLb3),sd(SparseLb3)/sqrt(iterations))))
}




