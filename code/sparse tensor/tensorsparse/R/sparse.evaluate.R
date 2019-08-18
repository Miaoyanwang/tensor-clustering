#' Evaluate the accuracy of clustering result
#' 
#' Given the input tensor, perform tensor clustering and evaluate the accuracy of the clustering result. 
#' @param bires the return value of label2()
#' @param data the return value of get.data()
#' @param CER logic value. If true, it would return CER also
#' @param show logic value. If true, it would print the result
#' @export

sparse.evaluate = function(bires, data, CER=TRUE, show=TRUE){
  npq=dim(data$x);n=npq[1];p=npq[2];q=npq[3]
  restensor<-(abs(bires$judgeX)>1e-10)*1
  totalzero<-sum(restensor==0)/(n*p*q)
  correctzero<-sum(restensor[which(data$binaryX==0)]==0)/sum(data$binaryX==0)
  correctone<-sum(restensor[which(data$binaryX==1)]==1)/(n*p*q-sum(data$binaryX==0))
  totalincorrect<-sum(abs(restensor-data$binaryX))/(n*p*q)
  error = sum((bires$judgeX-data$truthX)^2)/prod(dim(data$x))
  result = list(sparsityrate=totalzero, correctzerorate=correctzero,correctonerate=correctone,totalincorrectrate=totalincorrect,error=error)
  if (CER == TRUE){
    cerC<-1-adjustedRand(data$truthCs,bires$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(data$truthDs,bires$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(data$truthEs,bires$Es,randMethod=c("Rand"))
     blocktruth=cluster2block(data$mus,data$truthCs,data$truthDs,data$truthEs)
     blockest=cluster2block(bires$mus,bires$Cs,bires$Ds,bires$Es)
    
     cerTotal<-1-adjustedRand(as.numeric(as.factor(blocktruth)),as.numeric(as.factor(blockest)),randMethod=c("Rand"))
     #if (show == TRUE) cat("CerC is", cerC, ", cerD is", cerD, ", cerE is", cerE, ", cerTotal is", cerTotal, ".\n")
     if (show == TRUE) cat("CerC is", cerC, ", cerD is", cerD, ", cerE is", cerE, ".\n")
    result = list(sparsityrate=totalzero, correctzerorate=correctzero,correctonerate=correctone,totalincorrectrate=totalincorrect,cerC=cerC,cerD=cerD,cerE=cerE,cerTotal=cerTotal,error=error)
  }
  if(show == TRUE) cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect, ", MSE is", error,".\n")
  return(result)
}
