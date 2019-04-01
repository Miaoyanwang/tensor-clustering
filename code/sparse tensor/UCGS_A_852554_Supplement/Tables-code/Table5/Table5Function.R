library(sparseBC)
library("s4vd")
library("biclust")



# Function to calculate the results
Convert<-function(X,member){
membermatrix<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  for(i in 1:length(member))	
  {
    membermatrix[member[[i]][[1]],member[[i]][[2]]]<-1
  } 	
  return(membermatrix)	
}
###########################################
# New Simulation as in SSVD paper Mihee Lee
# Function for the papers simulation
###########################################

Do<-function(iteration,method,lambda){
	
u<-c(10,9,8,7,6,5,4,3,rep(2,17),rep(0,75))
v<-c(10,-10,8,-8,5,-5,rep(3,5),rep(-3,5),rep(0,34))
u<-u/sqrt(sum(u^2))
v<-v/sqrt(sum(v^2))
d<-50
mus<-(d*u%*%t(v))
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
    KR<-sparseBC.choosekr(X,1:6,1:6,0,0.1)
    k<-KR$estimated_kr[1]
    r<-KR$estimated_kr[2]
    
      lambda2<-lambda
      biclustering<-sparseBC(X,k,r,lambda2)
      resmatrix<-(abs(biclustering$mus)>1e-10)
    
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
    ssvd<-biclust(X,BCssvd(),K=1,gamu=2,gamv=2)
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
return(list(totalzero=mean(totalzero)/5000,sdtotalzero=sd(totalzero/5000)/sqrt(100),correctzero=mean(correctzero)/4600,sdcorrectzero=sd(correctzero/4600)/sqrt(100),correctone=mean(correctone)/400,sdcorrectone=sd(correctone/400)/sqrt(100),totalincorrect=mean(totalincorrect)/5000,sdtotalincorrect=sd(totalincorrect/5000)/sqrt(100)))
}
