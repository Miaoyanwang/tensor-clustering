#####################################################
# Table 7 - Evaluation of biclustering methods under the assumption of matrix variate normal
# The mean matrix is sparse
# Simulation for Various clustering methods on multiple scenarios
#####################################################

library(sparseBC)
library(s4vd)
library(clues)
library(biclust)

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



Convert<-function(X,member){
membermatrix<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  for(i in 1:length(member))	
  {
    membermatrix[member[[i]][[1]],member[[i]][[2]]]<-1
  } 	
  return(membermatrix)	
}

Do<-function(n,p,K,R,iteration,lambda,method,rank,standarddeviation,alpha,beta,max.iter){
  
cerrow<-c(rep(NA,iteration))	
cercol<-c(rep(NA,iteration))
totalzero<-c(rep(NA,iteration))
correctzero<-c(rep(NA,iteration))
correctone<-c(rep(NA,iteration))
totalincorrect<-c(rep(NA,iteration))

  for(iter in 1:iteration){
    set.seed(iter)
    cat("iteration ",iter,fill=TRUE)
  
    Matrix<-CreateData(n,p,K,R,standarddeviation)
    x<-Matrix$matrix
    Sigma<-Matrix$Sigma
    Delta<-Matrix$Delta
    mus<-Matrix$truthmus
    Cs<-Matrix$truthCs
    Ds<-Matrix$truthDs 	
    binaryX<-(mus!=0)*1
    kmeansCs<-kmeans(x,K,nstart=20)$cluster
  	kmeansDs<-kmeans(t(x),R,nstart=20)$cluster
  	if(method=="kmeans"){
  	  #kmeansCs<-kmeans(x,K,nstart=20)$cluster
  	  #kmeansDs<-kmeans(t(x),R,nstart=20)$cluster
  	  cerrow[iter]<-1-adjustedRand(kmeansCs,Cs,randMethod=c("Rand"))
  	  cercol[iter]<-1-adjustedRand(kmeansDs,Ds,randMethod=c("Rand"))	
  	  resmatrix<-binaryX*0
  	  }
  	  
	if(method=="biclustering"){
	  biclustering<-sparseBC(x,K,R,lambda,nstart=20,Cs.init=kmeansCs,Ds.init=kmeansDs)
	  resmatrix<-(abs(biclustering$mus)>1e-6)
	  cerrow[iter]<-1-adjustedRand(biclustering$Cs,Cs,randMethod=c("Rand"))
      cercol[iter]<-1-adjustedRand(biclustering$Ds,Ds,randMethod=c("Rand"))
	  }

	if(method=="matrix"){
	  matrixbiclust<-matrixBC(x,K,R,lambda,alpha,beta,nstart=20,Cs.init=kmeansCs,Ds.init=kmeansDs,max.iter=max.iter)
	  resmatrix<-(abs(matrixbiclust$mus)>1e-6)
	  cerrow[iter]<-1-adjustedRand(matrixbiclust$Cs,Cs,randMethod=c("Rand"))
      cercol[iter]<-1-adjustedRand(matrixbiclust$Ds,Ds,randMethod=c("Rand"))	  	  
	  }

	if(method=="matrixknown"){
	  matrixbiclust<-matrixBC(x,K,R,lambda,alpha,beta,nstart=20,Cs.init=kmeansCs,Ds.init=kmeansDs,max.iter=max.iter,Sigma.init=Sigma,Delta.init=Delta)
	  resmatrix<-(abs(matrixbiclust$mus)>1e-6)
	  cerrow[iter]<-1-adjustedRand(matrixbiclust$Cs,Cs,randMethod=c("Rand"))
      cercol[iter]<-1-adjustedRand(matrixbiclust$Ds,Ds,randMethod=c("Rand"))	  	  
	  }
	  
  if(method=="plaid1"){
    plaid<-biclust(x,BCPlaid(),fit.model=y~m+a+b,verbose=FALSE,row.release=0.5,col.release=0.5,background=TRUE)
    plaid.member<-biclusternumber(plaid)
    resmatrix<-Convert(x,plaid.member)
    }

  if(method=="plaid2"){
    plaid<-biclust(x,BCPlaid(),fit.model=y~m,verbose=FALSE,row.release=0.5,col.release=0.5,background=TRUE)
    plaid.member<-biclusternumber(plaid)
    resmatrix<-Convert(x,plaid.member)
    }

  if(method=="ssvd"){
    ssvd<-biclust(x,BCssvd(),K=rank,gamu=2,gamv=2)
    ssvd.member<-biclusternumber(ssvd)
    resmatrix<-Convert(x,ssvd.member)
    }
    
  if(method=="LAS"){
	write.table(x,file=paste("LAS",iter,".txt",sep=""),sep="\t",col.names=FALSE,row.names=FALSE)
	write.table(binaryX,file=paste("binaryX",iter,".txt",sep=""),sep="\t",col.names=FALSE,row.names=FALSE)
 
    resmatrix<-binaryX*0
  	
  	}
  
totalzero[iter]<-sum(resmatrix==0)/(n*p)
correctzero[iter]<-sum(resmatrix[which(binaryX==0)]==0)/sum(binaryX==0)
correctone[iter]<-sum(resmatrix[which(binaryX==1)]==1)/(n*p-sum(binaryX==0))
totalincorrect[iter]<-sum(abs(resmatrix-binaryX))/(n*p)




}
return(list(cerrow=c(mean(cerrow),sd(cerrow)/sqrt(iteration)),cercol=c(mean(cercol),sd(cercol)/sqrt(iteration)),totalzero=c(mean(totalzero),sd(totalzero)/sqrt(iteration)),correctzero=c(mean(correctzero),sd(correctzero)/sqrt(iteration)),correctone=c(mean(correctone),sd(correctone)/sqrt(iteration)),totalincorrect=c(mean(totalincorrect),sd(totalincorrect)/sqrt(iteration))))	
}