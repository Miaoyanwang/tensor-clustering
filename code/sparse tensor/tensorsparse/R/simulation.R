#' Perform simulation: calculating the mean CER over iterations
#' 
#' Perform simulation: calculating the mean CER over iterations.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param error the standard deviation when producing data
#' @param lambda ...
#' @param iteration iteration times
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @export
simulation  = function(n,p,q,k,r,l,error,lambda,iteration=1,method="L0"){
  cer = c()
  for (i in 1:iteration){
    set.seed(i)
    data = get.data(n,p,q,k,r,l,error=error)
    test = data$x
    truthCs = data$truthCs
    truthDs = data$truthDs
    truthEs = data$truthEs
    sim = classify2(test,k,r,l,method=method)
    cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
    cer = c(cer,cerC,cerD,cerE)
  }
  meancer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,mean)
  sdcer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,sd)
  return(list("meancer"=meancer,"sdcer"=sdcer))
}
