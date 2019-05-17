#' calculate clustering error rate
#' 
#' calculate clustering error rate
#' @param bires the return list of label2()
#' @param data the return list of get.data()
#' 
#' @return the cer in all modes.
#' 
#' @export
cer = function(bires,data){
  cerC<-1-adjustedRand(data$truthCs,bires$Cs,randMethod=c("Rand"))
  cerD<-1-adjustedRand(data$truthDs,bires$Ds,randMethod=c("Rand"))
  cerE<-1-adjustedRand(data$truthEs,bires$Es,randMethod=c("Rand"))
  cat("The CER(clustering error rate) is ",cerC,",",cerD,",",cerE,".\n")
  return(list(cerC = cerC, cerD = cerD, cerE = cerE))
}
