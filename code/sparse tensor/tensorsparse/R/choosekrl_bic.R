#' choose the best k,r,l to do clustering when the ranges are given
#' 
#' choose the best k,r,l to do clustering when the ranges are given by minimizing the BIC.
#' @param x a three-dimensional array
#' @param k a vector, the possible clusters number of mode 1
#' @param r a vector, the possible clusters number of mode 2
#' @param l a vector, the possible clusters number of mode 3
#' @param lambda a numeric value
#' @param sim.times the same as label2()
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @return a list \code{estimated_krl} ...\\
#'                \code{BIC} ...
#' 
choosekrl_bic = function (x,k,r,l,lambda=0,sim.times=1,method="L0"){
  #k = 2:5;r=2:5;l=2:5;lambda=0;sim.times=1;method="L0"
    ## x = x - mean(x) ## commented out by Miaoyan
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
