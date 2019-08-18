#' Perform simulation: selecting the lambda
#' 
#' 
#' Simulation: selecting the lambda.  Select the lambda with lowest BIC while given a certain range.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param sparse the sparse percent of data
#' @param iteration iteration times
#' @param lambda a vector of possible lambda, for example: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)
#' @param standarddeviation the standard deviation when producing data
#' @param center if True, then x = x- mean(x)
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
#' @return \code{selectedlambda}: a vector of selecting lambdas with lowest BIC in each iteration.
#' 
#' @export
sim.chooseLambda = function(n,p,q,k,r,l,sparse,iteration,lambda,standarddeviation=4,center = FALSE,method="L0"){
  selectedLambda = 1:iteration
  for(iter in 1:iteration){
    cat("Iteration",iter,fill=TRUE)
    set.seed(iter)
    smp = get.data(n,p,q,k,r,l,error=standarddeviation,sort=FALSE,sparse.percent = sparse)$x
    if(center == TRUE)smp = smp - mean(smp)
    selectedLambda[iter] = chooseLambda(smp,k,r,l,lambda=lambda,method=method)$lambda
  }
  return(selectedLambda)
}
