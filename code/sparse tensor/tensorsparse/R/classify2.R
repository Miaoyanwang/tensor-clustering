#' Perform tensor clustering via TBM method
#' 
#' This function performs sparse clustering on a three-dimensional tensor via TBM method.
#' @param x a three-dimensional array
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param lambda a positive numeric value. The coefficient of the regularized term.
#' @param max.iter a positive integer. The Maximum times of iteration.
#' @param threshold a positive small numeric value which determines whether the algorithm converges or not.
#' @param trace logic value. If true, it would print the iteration situation.
#' @param Cs.init vector or NULL. Initial clsuter result of mode 1.
#' @param Ds.init vector or NULL. Initial clsuter result of mode 2.
#' @param Es.init vector or NULL. Initial clsuter result of mode 3.
#' @param nstart positive interger. The same as the "nstart" in kmeans().
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
#' @param center logic value that indicates whether run "x = x-mean(x)" before performing clustering.
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
#' 
#' @export
classify2 = function(x,k,r,l,lambda=0,max.iter=1000,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,nstart=20,method="L0",center=FALSE){
  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
  if(center == TRUE) {
    mustemp <- mean(x)
    x <- x-mustemp
  }
  if(is.null(Cs.init)){
    if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor.unfold(x,1),k,nstart = nstart)$cluster}
  } else {
    Cs = Cs.init
  }
  if(is.null(Ds.init)){
    if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor.unfold(x,2),r,nstart = nstart)$cluster}
  } else {
    Ds = Ds.init
  }
  if(is.null(Es.init)){
    if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor.unfold(x,3),l,nstart = nstart)$cluster}
  } else {
    Es = Es.init
  }
  objs <- 1e+15
  improvement <- 1e+10
  i <- 1
  mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
  while((improvement > threshold) & (i <= max.iter)){
    Cs = UpdateClusters.tensor(tensor.unfold(x),tensor.unfold(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Cs <- ReNumber(Cs)
    k = length(unique(Cs))
    #mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds = UpdateClusters.tensor(tensor.unfold(x,2),tensor.unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    #mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Es = UpdateClusters.tensor(tensor.unfold(x,3),tensor.unfold(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Es <- ReNumber(Es)
    l = length(unique(Es))
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda, method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    improvement <- abs(objs[length(objs)] - objs[length(objs) - 
                                                   6])/abs(objs[length(objs) - 6])
    i <- i + 1
    if(trace) cat("step",i,",improvement=",improvement,".\n")
    if (is.na(improvement)) break
  }
  if (i > max.iter) {
    warning("The algorithm has not converged by the specified maximum number of iteration.\n")
  }
  if(center==TRUE){
    mu.array <- mu.array + mustemp
  }
  mu.array[abs(mu.array)<=1e-6] = 0
  return(list("judgeX"=mu.array[Cs,Ds,Es, drop=FALSE],"Cs"=Cs,"Ds"=Ds,"Es"=Es,"objs"=objs[length(objs)], "mus"=mu.array))
}
