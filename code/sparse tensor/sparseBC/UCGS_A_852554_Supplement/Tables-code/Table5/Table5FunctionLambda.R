#####################################################
# Table 5 - Evaluation of biclustering methods under the assumption of ssvd
# The mean matrix is sparse
# This is a rank 1 matrix 
# Simulation for Various clustering methods on multiple scenarios
# Biclustering with lambda automatically chosen with the BIC criterion
#####################################################

library(sparseBC)

Do<-function(iteration,lambda){
	
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
selectedlambda<-NULL

for(i in 1:iteration){
  cat("iteration",i,fill=TRUE)
  set.seed(i)
  X<-mus+matrix(rnorm(100*50),100,50)
  X<-X-mean(X)


    KR<-sparseBC.choosekr(X,1:6,1:6,0,0.1)
    k<-KR$estimated_kr[1]
    r<-KR$estimated_kr[2]

      lambda2<-sparseBC.BIC(X,k,r,lambda)$lambda
      biclustering<-sparseBC(X,k,r,lambda2)
      resmatrix<-(abs(biclustering$mus)>1e-6)



  selectedlambda<-c(selectedlambda,lambda2)
  totalzero<-c(totalzero,sum(resmatrix==0))
  correctzero<-c(correctzero, sum(resmatrix[which(binaryX==0)]==0))
  correctone<-c(correctone,sum(resmatrix[which(binaryX==1)]==1))
  totalincorrect<-c(totalincorrect,sum(abs(resmatrix-binaryX)))

  }
return(list(totalzero=mean(totalzero)/5000,sdtotalzero=sd(totalzero/5000)/sqrt(100),correctzero=mean(correctzero)/4600,sdcorrectzero=sd(correctzero/4600)/sqrt(100),correctone=mean(correctone)/400,sdcorrectone=sd(correctone/400)/sqrt(100),totalincorrect=mean(totalincorrect)/5000,sdtotalincorrect=sd(totalincorrect/5000)/sqrt(100),lambda=mean(selectedlambda)))
}
