#' perform tensor clustering
#' 
#' perform tensor clustering
#' @param x a three-dimensional array
#' @param k the number of clusters in mode 1
#' @param r the number of clusters in mode 2
#' @param l the number of clusters in mode 3
#' @param lambda the sparsity parameter in the sparse clustering
#' @param max.iter maximum times of iteration
#' @param threshold ...
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
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds = UpdateClusters.tensor(tensor.unfold(x,2),tensor.unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
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
