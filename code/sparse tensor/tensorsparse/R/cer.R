#' Calculate the clustering error rate (CER)
#' 
#' Calculate clustering error rate (CER) by comparing the clustering result and the true underlying mean signal tensor.
#' @param bires List. The return value of label2() or classify2() which consists of judgeX, Cs, Ds, Es, objs, mus.
#' @param data List. The return value of get.data().
#' 
#' @return A list which consists of the CER of all three modes.
#' 
#' @export
cer = function(bires,data){
  cerC<-1-adjustedRand(data$truthCs,bires$Cs,randMethod=c("Rand"))
  cerD<-1-adjustedRand(data$truthDs,bires$Ds,randMethod=c("Rand"))
  cerE<-1-adjustedRand(data$truthEs,bires$Es,randMethod=c("Rand"))
  cat("The CER(clustering error rate) is ",cerC,",",cerD,",",cerE,".\n")
  return(list(cerC = cerC, cerD = cerD, cerE = cerE))
}
