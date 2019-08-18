#' Perform tuning parameter (\eqn{d_1}, \eqn{d_2}, \eqn{d_3}) selection for sparse tensor clustering via BIC criterion
#' 
#' Select the best \eqn{d_1}, \eqn{d_2}, \eqn{d_3}  to perform clustering. A range of values of \eqn{d_1}, \eqn{d_2}, \eqn{d_3} is usually considered - value that results in the lowest BIC is selected.
#' @param x a three-dimensional array
#' @param k the range of \eqn{d_1}: a vector, the possible clusters numbers of mode 1
#' @param r the range of \eqn{d_2}: a vector, the possible clusters numbers of mode 2
#' @param l the range of \eqn{d_3}: a vector, the possible clusters numbers of mode 3
#' @param lambda a numeric value. The coefficient of the regularization term.
#' @param sim.times the simulation times when perform clustering of classify2() in label2().
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty, "L1" indicating Lasso penalty.
#' @return a list   
#' 
#' \code{estimated_krl} a 1*3 matrix which is the estimated c(d_1,d_2,d_3).   
#' 
#'                \code{BIC} a vector which contains the BIC values of all the combination of given range.  
#' 
#' @export
choosekrl_bic = function (x,k,r,l,lambda=0,sim.times=1,method="L0"){
  #k = 2:5;r=2:5;l=2:5;lambda=0;sim.times=1;method="L0"
    ## x = x - mean(x) ## commented out
  if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0 || sum(diff(l) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
  n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
  krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
                 rep(1:length(r),times=length(k)*length(l)),
                 rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
               nrow=3)
  krl_list = as.list(as.data.frame(krl))
  if (.Platform$OS.type == "windows") {
    bires = apply(krl,MARGIN=2,label_for_krl,k,r,l,sim.times=sim.times,lambda=lambda,xmiss=x,method=method,crossvalidation=FALSE)
    CBIC = lapply(bires,tensor.calculateBIC, x=x,method=method)
    BIC = unlist(CBIC)
  } else {
    bires = mclapply(krl_list, label_for_krl,k,r,l,sim.times=sim.times,lambda=lambda,xmiss=x,method=method,crossvalidation=FALSE,mc.cores=n.cores)
    CBIC = mclapply(bires,tensor.calculateBIC, x=x,method=method,mc.cores = n.cores)
    BIC = unlist(CBIC)
  }
  names(BIC) = apply(krl,MARGIN=2,FUN=function(x)paste(k[x[1]],r[x[2]],l[x[3]]))
  best = krl_list[[which(BIC == min(BIC))[1]]]
  return(list(estimated_krl = t(as.matrix(c(k[best[1]],r[best[2]],l[best[3]]))), BIC = BIC))
}
