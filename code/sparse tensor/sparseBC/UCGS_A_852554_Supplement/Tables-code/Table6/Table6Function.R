#####################################################
# Table 6 - Evaluation of biclustering methods under the assumption of ssvd
# The mean matrix is sparse
# This is a rank 2 matrix with overlapping biclusters
# Simulation for Various clustering methods on multiple scenarios
#####################################################
library(sparseBC)
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

Do<-function(iteration,method,lambda){
	
u1<-c(10,9,8,7,6,5,4,3,rep(2,17),rep(0,75))
v1<-c(10,-10,8,-8,5,-5,rep(3,5),rep(-3,5),rep(0,34))
u1<-u1/sqrt(sum(u1^2))
v1<-v1/sqrt(sum(v1^2))
u2<-c(rep(0,13),10,9,8,7,6,5,4,3,rep(2,17),rep(0,62))
v2<-c(rep(0,9),10,-9,8,-7,6,-5,rep(4,5),rep(-3,5),rep(0,25))
u2<-u2/sqrt(sum(u2^2))
v2<-v2/sqrt(sum(v2^2))
d<-50
mus<-(d*u1%*%t(v1))+(d*u2%*%t(v2))
binaryX<-(mus!=0)*1

totalzero<-NULL
correctzero<-NULL
correctone<-NULL
totalincorrect<-NULL

for(i in 1:iteration){
  cat("iteration",i,fill=TRUE)
  set.seed(i)
  X<-mus+matrix(rnorm(100*50),100,50)
  X<-X-mean(X)

  if(method=="biclustering"){
    KR<-sparseBC.choosekr(X,1:10,1:10,0,0.1)
    k<-KR$estimated_kr[1]
    r<-KR$estimated_kr[2]
    
      lambda2<-lambda
      biclustering<-sparseBC(X,k,r,lambda2)
      resmatrix<-(abs(biclustering$mus)>1e-15)

}

  if(method=="plaid1"){
    plaid<-biclust(X,BCPlaid(),fit.model=y~m+a+b,verbose=FALSE,row.release=0.5,col.release=0.5,background=FALSE)
    plaid.member<-biclusternumber(plaid)
    resmatrix<-Convert(X,plaid.member)
    }


  if(method=="plaid2"){
    plaid<-biclust(X,BCPlaid(),fit.model=y~m,verbose=FALSE,row.release=0.5,col.release=0.5,background=FALSE)
    plaid.member<-biclusternumber(plaid)
    resmatrix<-Convert(X,plaid.member)
    }

  if(method=="ssvd"){
    ssvd<-biclust(X,BCssvd(),K=2,gamu=2,gamv=2)
    ssvd.member<-biclusternumber(ssvd)
    resmatrix<-Convert(X,ssvd.member)
    }
  if(method=="LAS"){
    write.table(X,file=paste("LAS",i,".txt",sep=""),sep="\t",col.names=FALSE,row.names=FALSE) 
    resmatrix<-matrix(0,100,50)
  	}

  totalzero<-c(totalzero,sum(resmatrix==0))
  correctzero<-c(correctzero, sum(resmatrix[which(binaryX==0)]==0))
  correctone<-c(correctone,sum(resmatrix[which(binaryX==1)]==1))
  totalincorrect<-c(totalincorrect,sum(abs(resmatrix-binaryX)))

  }
return(list(totalzero=mean(totalzero)/5000,sdtotalzero=sd(totalzero/5000)/sqrt(100),correctzero=mean(correctzero)/4284,sdcorrectzero=sd(correctzero/4284)/sqrt(100),correctone=mean(correctone)/716,sdcorrectone=sd(correctone/716)/sqrt(100),totalincorrect=mean(totalincorrect)/5000,sdtotalincorrect=sd(totalincorrect/5000)/sqrt(100)))
}
