#' Perform tensor clustering via TBM method
#' 
#' Perform tensor clustering with TBM method in a more stable way. Repeat the classify2() many times and select the one with lowest MSE.
#' @param x a three-dimensional array
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param lambda a positive numeric value. The coefficient of the regularized term.
#' @param max.iter a positive integer. The Maximum times of iteration.
#' @param threshold a positive small numeric value which determines whether the algorithm converges or not.
#' @param sim.times the times of calling classify2() with different seeds. 
#' @param trace logic value. If true, it would print the iteration situation.
#' @param Cs.init vector or NULL. Initial clsuter result of mode 1.
#' @param Ds.init vector or NULL. Initial clsuter result of mode 2.
#' @param Es.init vector or NULL. Initial clsuter result of mode 3.
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
#' @return a list   
#' \code{judgeX} estimated underlying mean signal.  
#' 
#'                \code{Cs} clustering result of mode 1.  
#'                
#'                \code{Ds} clustering result of mode 2.  
#'                
#'                \code{Es} clustering result of mode 3.  
#'                
#'                \code{mus} estimated underlying mean signal of each cluster.  
#'                
#' @export
label2 = function(x,k,r,l,lambda=0,max.iter=1000,threshold = 1e-10,sim.times=1,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,method="L0"){
  #x=test;lambda=1e-3;max.iter=200;threshold = 5e-3;sim.times=10
  if (sim.times == 1) return(classify2(x,k,r,l,lambda=lambda,max.iter = max.iter,threshold = threshold,Cs.init = Cs.init,Ds.init = Ds.init,Es.init = Es.init,method=method))
  if (.Platform$OS.type == "windows") {
    result = lapply(rep(list(x),sim.times), classify2, k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,method=method)
    objs = unlist(lapply(result, function(result){result$objs}))
  } else {
    result = mclapply(rep(list(x),sim.times), classify2, k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,nstart = sample(1:1000,1),method=method,mc.cores = n.cores)
    objs = unlist(lapply(result, function(result){result$objs}))
  }
  result = result[[which(objs == min(objs))[1]]]
  return(result)
}
