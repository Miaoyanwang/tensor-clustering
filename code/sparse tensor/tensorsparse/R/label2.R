#' perform tensor clustering
#' 
#' perform tensor clustering in a more stable way
#' @param x a three-dimensional array
#' @param k the number of clusters in mode 1
#' @param r the number of clusters in mode 2
#' @param l the number of clusters in mode 3
#' @param lambda the sparsity parameter in the sparse clustering
#' @param max.iter maximum times of iteration
#' @param threshold ...
#' @param sim.times the times of calling classify2() with different seeds.
#' @param trace ...
#' @param Cs.init ...
#' @param Ds.init ...
#' @param Es.init ...
#' @param nstart ...
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicates Lasso penalty, and "L0" indicates sparse subset penalty
#' @param center ...
#' @return a list \code{judgeX}
#'                \code{Cs} clustering result in mode 1 
#'                \code{Ds} clustering result in mode 2
#'                \code{Es} clustering result in mode 3
#'                \code{mus}
#' 
#' 
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
