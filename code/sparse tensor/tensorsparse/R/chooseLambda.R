#' Perform tuning parameter (lambda) selection for sparse tensor clustering via BIC criterion
#' 
#' We assume that \eqn{d_1}, \eqn{d_2}, \eqn{d_3} are known. A range of values of lambda is usually considered - value that results in the lowest BIC is selected.
#' @param x a three-dimensional array
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param lambda a vector of possible lambda, eg: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200).
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
#' @return a list   
#'                \code{lambda} the lambda with lowest BIC;  
#' 
#'                \code{BIC} the corresponding BIC for each lambda in the given range;  
#'                
#'                \code{nonzeromus} the number clusters with non-zero mean.
#'   
#' @export
chooseLambda = function (x, k, r, l, lambda=NULL,method="L0") {
  ##x = x - mean(x) commented out
  if (is.null(lambda)){
    n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
    if (method == "L0") lambda = sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.1)
    if (method == "L1") lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.1)
    if (is.null(lambda)) stop("No such kind of method:", method, ".\n")
  } 
  if (.Platform$OS.type == "windows") {
    bires = lapply(lambda,FUN=classify2,x=x,k=k,r=r,l=l,method=method)
    CBIC = lapply(bires,tensor.calculateBIC, x=x,method=method)
    BIC = unlist(CBIC)
    nonzero = unlist(lapply(bires, FUN=function(bires){return(sum(bires$mus!=0))}))
  } else {
    bires = mclapply(lambda,FUN=classify2,x=x,k=k,r=r,l=l,method=method,mc.cores = n.cores)
    CBIC = mclapply(bires,tensor.calculateBIC, x=x,method=method,mc.cores = n.cores)
    BIC = unlist(CBIC)
    nonzero = unlist(mclapply(bires, FUN=function(bires){return(sum(bires$mus!=0))},mc.cores = n.cores))
  }
  return(list(lambda = lambda[which(BIC == min(BIC))[1]], BIC = BIC, 
              nonzeromus = nonzero))
}
