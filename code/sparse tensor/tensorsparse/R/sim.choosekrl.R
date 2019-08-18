#' Perform simulation: selecting the best \eqn{d_1}, \eqn{d_2} and \eqn{d_3}
#' 
#' Simulation: selecting the best \eqn{d_1}, \eqn{d_2} and \eqn{d_3}.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param error positive numeric value. The noise when producing the data
#' @param sim.times positive integer. Simulation times
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
#' @param mode two options: "bic": selecting cluster numbers by bic; "crossvalidation":selecting cluster numbers by cross validation.
#' @param seed logic value. Whether set seed to each simulation.
#' @return A list consists of the estimating result in different iteration.
#' 
#' @export

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
