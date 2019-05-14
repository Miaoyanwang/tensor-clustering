#####################################################
# Table 4 - Evaluation of biclustering methods
# The mean matrix is sparse in this case
# Simulation for various clustering methods on multiple scenarios
#####################################################
library(sparseBC)
library(clues)
library("s4vd")
library("biclust")

Convert<-function(X,member){
membermatrix<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  for(i in 1:length(member))	
  {
    membermatrix[member[[i]][[1]],member[[i]][[2]]]<-1
  } 	
  return(membermatrix)	
}

Do<-function(n,p,k,r,iteration,lambda,method,rank,standarddeviation=4){
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
  mus<-runif(k*r,1.5,2.5)*sample(c(-1,1),k*r,replace=TRUE)*rbinom(k*r,1,prob=0.5)
  mus<-matrix(c(mus),nrow=k,ncol=r,byrow=FALSE)
  truthCs<-sample(1:k,n,rep=TRUE)
  truthDs<-sample(1:r,p,rep=TRUE)
  x<-matrix(rnorm(n*p,mean=0,sd=standarddeviation),nrow=n,ncol=p)  
	
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
	
  if(method=="kmeans"){
  	  kmeanres<-kmeans(x,k,nstart=20)
  	  kmeanres2<-kmeans(t(x),r,nstart=20)
      cerrow[iter]<-1-adjustedRand(truthCs,kmeanres$cluster,randMethod=c("Rand"))
      cercol[iter]<-1-adjustedRand(truthDs,kmeanres2$cluster,randMethod=c("Rand"))
	  resmatrix<-binaryX*0  	  
  }	

  if(method=="biclustering"){
      biclustering<-sparseBC(x,k,r,lambda)
      resmatrix<-(abs(biclustering$mus)>1e-6)
      cerrow[iter]<-1-adjustedRand(truthCs,biclustering$Cs,randMethod=c("Rand"))
      cercol[iter]<-1-adjustedRand(truthDs,biclustering$Ds,randMethod=c("Rand"))
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
return(list(cerrow=c(mean(cerrow),sd(cerrow)/sqrt(iteration)),cercol=c(mean(cercol),sd(cercol)/sqrt(iteration)),totalzero=c(mean(totalzero),sd(totalzero)/sqrt(iteration)),correctzero=c(mean(correctzero),sd(correctzero)/sqrt(iteration)),correctone=c(mean(correctone),sd(correctone)/sqrt(iteration)),totalincorrect=c(mean(totalincorrect),sd(totalincorrect)/sqrt(iteration)),selectedlambda=selectedlambda))	
}