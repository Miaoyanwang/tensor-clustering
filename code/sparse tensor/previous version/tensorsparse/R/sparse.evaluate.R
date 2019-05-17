#' evaluate the accuracy of the model when inputting a sparse tensor
#' 
#' evaluate the accuracy of the model when inputting a sparse tensor
#' @param bires the return list of label2()
#' @param data the return list of get.data()

sparse.evaluate = function(bires, data){
  npq=dim(data$x);n=npq[1];p=npq[2];q=npq[3]
  restensor<-(abs(bires$judgeX)>1e-10)*1
  totalzero<-sum(restensor==0)/(n*p*q)
  correctzero<-sum(restensor[which(data$binaryX==0)]==0)/sum(data$binaryX==0)
  correctone<-sum(restensor[which(data$binaryX==1)]==1)/(n*p*q-sum(data$binaryX==0))
  totalincorrect<-sum(abs(restensor-data$binaryX))/(n*p*q)
  cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect,".\n")
}