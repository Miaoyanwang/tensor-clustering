#' Perform tensor clustering via TBM method
#' 
#' This function performs sparse clustering on a three-dimensional tensor via TBM method.
#' @param x a four-dimensional array
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param m \eqn{d_4}: the clusters number of mode 4
#' @param lambda a positive numeric value. The coefficient of the regularized term.
#' @param max.iter a positive integer. The Maximum times of iteration.
#' @param threshold a positive small numeric value which determines whether the algorithm converges or not.
#' @param trace logic value. If true, it would print the iteration situation.
#' @param Cs.init vector or NULL. Initial clsuter result of mode 1.
#' @param Ds.init vector or NULL. Initial clsuter result of mode 2.
#' @param Es.init vector or NULL. Initial clsuter result of mode 3.
#' @param Fs.init vector or NULL. Initial clsuter result of mode 4.
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
#'                \code{Fs} clustering result of mode 4. 
#'                
#'                \code{mus} estimated underlying mean signal of each cluster.  
#' 
#' 
#' @export


classify4 = function(x,k,r,l,m,lambda=0,max.iter=1000,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,Fs.init=NULL,nstart=20,method="L0",center=FALSE){
  #x=try$x;lambda=0;max.iter=1000;threshold=1e-15;trace=FALSE;Cs.init=NULL;Ds.init=NULL;Es.init=NULL;nstart=20;method="L0";center=FALSE
  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]; s = dim(x)[4]
  if(center == TRUE) {
    mustemp <- mean(x)
    x <- x-mustemp
  }
  if(is.null(Cs.init)){
    if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor.unfold4(x,1),k,nstart = nstart)$cluster}
  } else {
    Cs = Cs.init
  }
  if(is.null(Ds.init)){
    if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor.unfold4(x,2),r,nstart = nstart)$cluster}
  } else {
    Ds = Ds.init
  }
  if(is.null(Es.init)){
    if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor.unfold4(x,3),l,nstart = nstart)$cluster}
  } else {
    Es = Es.init
  }
  if(is.null(Es.init)){
    if(m==1) Es = rep(1,s) else {Fs  = kmeans(tensor.unfold4(x,4),m,nstart = nstart)$cluster}
  } else {
    Fs = Fs.init
  }
  objs <- 1e+15
  improvement <- 1e+10
  i <- 1
  mu.array = UpdateMus.tensor4(x,Cs,Ds,Es,Fs,lambda,method=method)
  while((improvement > threshold) & (i <= max.iter)){
    #Cs = UpdateClusters.tensor(tensor.unfold4(x),tensor.unfold4(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
    Cs = UpdateClusters.tensor(tensor.unfold4(x),tensor.unfold4(mu.array),Cs,
                               rep(Ds,times=q*s)+r*(rep(rep(Es,each=p),times=s)-1)+r*l*(rep(Fs,each=p*q)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Cs <- ReNumber(Cs)
    k = length(unique(Cs))
    
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    #Ds = UpdateClusters.tensor(tensor.unfold4(x,2),tensor.unfold4(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    Ds = UpdateClusters.tensor(tensor.unfold4(x,2),tensor.unfold4(mu.array,2),Ds,
                               rep(Cs,times=q*s)+k*(rep(rep(Es,each=n),times=s)-1)+k*l*(rep(Fs,each=n*q)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    #Es = UpdateClusters.tensor(tensor.unfold4(x,3),tensor.unfold4(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    Es = UpdateClusters.tensor(tensor.unfold4(x,3),tensor.unfold4(mu.array,3),Es,
                               rep(Cs,times=p*s)+k*(rep(rep(Ds,each=n),times=s)-1)+k*r*(rep(Fs,each=n*p)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Es <- ReNumber(Es)
    l = length(unique(Es))
    
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    #Fs = UpdateClusters.tensor(tensor.unfold4(x,3),tensor.unfold4(mu.array,3),Fs,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    Fs = UpdateClusters.tensor(tensor.unfold4(x,4),tensor.unfold4(mu.array,4),Fs,
                               rep(Cs,times=p*q)+k*(rep(rep(Ds,each=n),times=q)-1)+k*r*(rep(Es,each=n*p)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Fs <- ReNumber(Fs)
    m = length(unique(Fs))
    
    mu.array = UpdateMus.tensor4(x,Cs,Ds,Es,Fs,lambda, method=method)
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    improvement <- abs(objs[length(objs)] - objs[length(objs) - 
                                                   7])/abs(objs[length(objs) - 7])
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
  return(list("judgeX"=mu.array[Cs,Ds,Es,Fs, drop=FALSE],"Cs"=Cs,"Ds"=Ds,"Es"=Es,"Fs"=Fs,"objs"=objs[length(objs)], "mus"=mu.array))
}