#####################################################
# Table 3 - Evaluation of biclustering methods
# The mean matrix is not sparse 
# Simulation for biclustering method with lambda automatically chosen as in Section 5.2
#####################################################
library(sparseBC)
library(clues)

Do<-function(n,p,k,r,iteration,lambda,standarddeviation=4){
# Initialize some variables	
cerrow<-c(rep(NA,iteration))	
cercol<-c(rep(NA,iteration))
totalzero<-c(rep(NA,iteration))
correctzero<-c(rep(NA,iteration))
correctone<-c(rep(NA,iteration))
totalincorrect<-c(rep(NA,iteration))
selectedlambda<-c(rep(NA,iteration))



for(iter in 1:iteration){
cat("Iteration",iter,fill=TRUE)
set.seed(iter)
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
  binaryX<-(truthmatrix!=0)*1
	
  x<-x-mean(x)
  
	  selectlambda<-sparseBC.BIC(x,k,r,lambda)$lambda
	  selectedlambda[iter]<-selectlambda
      biclustering<-sparseBC(x,k,r,selectlambda)
      resmatrix<-(abs(biclustering$mus)>1e-6)
      cerrow[iter]<-1-adjustedRand(truthCs,biclustering$Cs,randMethod=c("Rand"))
      cercol[iter]<-1-adjustedRand(truthDs,biclustering$Ds,randMethod=c("Rand"))
  
    
totalzero[iter]<-sum(resmatrix==0)/(n*p)
correctzero[iter]<-sum(resmatrix[which(binaryX==0)]==0)/sum(binaryX==0)
correctone[iter]<-sum(resmatrix[which(binaryX==1)]==1)/(n*p-sum(binaryX==0))
totalincorrect[iter]<-sum(abs(resmatrix-binaryX))/(n*p)

}
return(list(cerrow=c(mean(cerrow),sd(cerrow)/sqrt(iteration)),cercol=c(mean(cercol),sd(cercol)/sqrt(iteration)),totalzero=c(mean(totalzero),sd(totalzero)/sqrt(iteration)),correctzero=c(mean(correctzero),sd(correctzero)/sqrt(iteration)),correctone=c(mean(correctone),sd(correctone)/sqrt(iteration)),totalincorrect=c(mean(totalincorrect),sd(totalincorrect)/sqrt(iteration)),selectedlambda=selectedlambda))	
}