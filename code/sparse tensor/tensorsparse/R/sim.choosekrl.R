#' simulation: choose the best k,r,l to do clustering when the ranges are given
#' 
#' simulation: choose the best k,r,l to do clustering when the ranges are given
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k true clusters number of mode 1
#' @param r true clusters number of mode 2
#' @param l true clusters number of mode 3
#' @param error noise when producing the data
#' @param sim.times simulation times
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @param mode two options: "bic", "crossvalidation".
#' @param seed whether set seed to each simulation.
#' @return A list of the krl result in different iteration.
#' 


sim.choosekrl <- function(n,p,q,k,r,l,error=1,sim.times=5,method="L0",mode="bic",seed=TRUE){
  classification<-list()
  if(mode == "bic"){
    for(a in 1:sim.times){
      cat("Starting", a, fill=TRUE)
      if(seed == TRUE) set.seed(a)
      x = get.data(n,p,q,k,r,l,error=error)$x
      classification[[a]]<-choosekrl_bic(x,k=2:5,r=2:5,l=2:5,method=method)$estimated_krl#estimate clusters
    }
  }
  if(mode == "crossvalidation"){
    for(a in 1:sim.times){
      cat("Starting", a, fill=TRUE)
      if(seed == TRUE) set.seed(a)
      x = get.data(n,p,q,k,r,l,error=error)$x
      classification[[a]]<-sparse.choosekrl(x,k=2:5,r=2:5,l=2:5,method=method)$estimated_krl#estimate clusters
    }
  }
  if(mode!="bic" & mode!="crossvalidation") stop("No such kind of mode:", mode, ".\n")
  return(classification)
}