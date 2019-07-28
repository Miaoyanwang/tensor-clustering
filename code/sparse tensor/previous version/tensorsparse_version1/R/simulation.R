#' simulation: choose the best lambda when the range is given
#' 
#' choose the best lambda for a certain range.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k the clusters number of mode 1
#' @param r the clusters number of mode 2
#' @param l the clusters number of mode 3
#' @param error the standard deviation when producing data
#' @param lambda a vector of possible lambda, for example: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)
#' @param iteration iteration times
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
simulation  = function(n,p,q,k,r,l,error,lambda,iteration=1,method="L0"){
  cer = c()
  for (i in 1:iteration){
    data = get.data(n,p,q,k,r,l,error=error)
    test = data$x
    truthCs = data$truthCs
    truthDs = data$truthDs
    truthEs = data$truthEs
    sim = label2(test,k,r,l,threshold=5e-2,lambda=0,trace=FALSE,method=method)
    cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
    cer = c(cer,cerC,cerD,cerE)
  }
  cer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,mean)
  print(cer)
}