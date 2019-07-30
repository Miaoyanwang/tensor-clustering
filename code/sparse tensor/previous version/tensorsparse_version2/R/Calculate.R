#' calculate the correct rate of the result of sparse.choosekrl()
#' 
#' calculate the correct rate of the result of sparse.choosekrl()
#' @param true true k, r, l
#' @param results result returned by sparse.choosekrl()
#' 
#' @export
Calculate<-function(true,results){
  real<-matrix(true,ncol=3)
  percent<-0
  for(i in 1:length(results)){
    if(nrow(results[[i]])>1){
      for(a in 1:nrow(results[[i]])){
        #  	  cat("iteration",a,fill=TRUE)
        if(sum(results[[i]][a,]==real)==3){
          percent<-percent+(1/nrow(results[[i]]))
        }
      }
    }
    else if(nrow(results[[i]])<2){
      if(sum(results[[i]]==real)==3){
        percent<-percent+1
      }	
    }
  }
  return(percent/length(results))
}
