#####################################################
# Table 7 - Evaluation of biclustering methods under the assumption of matrix variate normal
# The mean matrix is sparse
# Simulation for MVN biclustering with unknown sigma and delta and with lambda automatically chosen
#####################################################
library(sparseBC)
library(clues)
################################
# Generate a positive definite covariance matrix
################################
Covariance<-function(p,sparsity,offdiagonal)
{
  sparse<-rbinom(p*p,1,1-sparsity)
  Sigma<-(matrix(data=sparse,nrow=p,ncol=p))
  Sigma[lower.tri(Sigma,diag=FALSE)]<-0 
  Sigma<-Sigma+t(Sigma)
  Sigma[which(Sigma==1)]<-offdiagonal
  Sigma<-(Sigma+t(Sigma))/2
  diag(Sigma)<-1
  ee<-min(eigen(Sigma,only.values=T)$values)
  diag(Sigma)<-diag(Sigma)+ifelse(ee < 0, -ee + 0.01, 0.01)
  return(Sigma)
}

################################
# Data Matrix from Matrix-variate distribution
################################

CreateData<-function(n,p,k,r,sd){		
  mus<-runif(k*r,-2,2)
#  sam<-sample(c(0,1),k*r,replace=TRUE)
  sam<-rbinom(k*r,1,0.5)
  for(i in 1:(k*r)){
  	if(sam[i]==0){
      mus[i]<-0  
  	}
  }
  mus<-matrix(c(mus),nrow=k,ncol=r,byrow=FALSE)
  truthCs<-rep(1:k,each=(n/k))
  truthDs<-rep(1:r,each=(p/r))

  truthmatrix<-matrix(NA,nrow=n,ncol=p)
  Sigma<-matrix(0,nrow=n,ncol=n)
  Delta<-matrix(0,nrow=p,ncol=p)

  for(i in 1:max(truthCs)){
	for(j in 1:max(truthDs)){
      truthmatrix[truthCs==i,truthDs==j]<-mus[i,j]
	} 
  }
  
  for(i in 1:k){
  	Cov<-Covariance(n/k,0.2,runif(100,-0.5,0.5))
  	Sigma[truthCs==i,truthCs==i]<-Cov
  	}
  	
  for(j in 1:r){
  	Cov<-Covariance(p/r,0.2,runif(100,-0.5,0.5))
  	Delta[truthDs==j,truthDs==j]<-Cov
  	}
  	  
  x<-matrix(rnorm(n*p,mean=0,sd=sd),nrow=n,ncol=p)
  eigenSigma<-eigen(Sigma)
  sqrtSigma<-eigenSigma$vectors%*%sqrt(diag(eigenSigma$values))%*%t(eigenSigma$vectors)
  eigenDelta<-eigen(Delta)
  sqrtDelta<-eigenDelta$vectors%*%sqrt(diag(eigenDelta$values))%*%t(eigenDelta$vectors)

  y<-sqrtSigma%*%x%*%sqrtDelta+truthmatrix
  
  return(list(truthmus=truthmatrix, truthCs=truthCs, truthDs=truthDs,matrix=y,Delta=Delta,Sigma=Sigma))		
}


Do<-function(n,p,k,r,iteration,lambda,standarddeviation=2,max.iter){
# Initialize some variables	
cerrow<-NULL
cercol<-NULL
totalzero<-NULL
correctzero<-NULL
correctone<-NULL
totalincorrect<-NULL
selectedlambda<-NULL

set.seed(iteration)
  Matrix<-CreateData(n,p,k,r,standarddeviation)	
  x<-Matrix$matrix
  mus<-Matrix$truthmus
  Sigma<-Matrix$Sigma
  Delta<-Matrix$Delta
  truthCs<-Matrix$truthCs
  truthDs<-Matrix$truthDs 	  
  binaryX<-(mus!=0)*1
  
	  selectlambda<-matrixBC.BIC(x,k,r,lambda,Sigma.init=Sigma,Delta.init=Delta)$lambda
      biclustering<-matrixBC(x,k,r,selectlambda,alpha=0.05,beta=0.05,max.iter=max.iter)
      resmatrix<-(abs(biclustering$mus)>1e-6)
      cerrow<-1-adjustedRand(truthCs,biclustering$Cs,randMethod=c("Rand"))
      cercol<-1-adjustedRand(truthDs,biclustering$Ds,randMethod=c("Rand"))
  
    
totalzero<-sum(resmatrix==0)/(n*p)
correctzero<-sum(resmatrix[which(binaryX==0)]==0)/sum(binaryX==0)
correctone<-sum(resmatrix[which(binaryX==1)]==1)/(n*p-sum(binaryX==0))
totalincorrect<-sum(abs(resmatrix-binaryX))/(n*p)


return(list(cerrow=cerrow,cercol=cercol,totalzero=totalzero,correctzero=correctzero,correctone=correctone,totalincorrect=totalincorrect,selectlambda=selectlambda))	
}


