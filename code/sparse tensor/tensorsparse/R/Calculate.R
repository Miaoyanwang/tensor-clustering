#' Calculate the correct rate of the result of selecting \eqn{d_1}, \eqn{d_2}, \eqn{d_3} 
#' 
#' Calculate the correct rate of the result of sparse.choosekrl().
#' @param true Vector. The true number of clusters in each mode (eg. c(3,3,3)).
#' @param results List. A list consists of 3*1 matrix and each matrix is the estimated c(d_1,d_2,d_3). (The return value of sparse.choosekrl().)
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
