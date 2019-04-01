#' generate a random tensor
#' 
#' generate a random tensor
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k the clusters number of mode 1 
#' @param r the clusters number of mode 2
#' @param l the clusters number of mode 3 
#' @param error noise in rnorm
#' @param sort if TRUE, the data belongs to the same clusters would assemble together
#' @param sparse.percent the proportion of 0 mean
#' @param center if center sum(x) = 0
#' @param mumin ...
#' @param mumax ...
#' @return a list \code{x}
#'                \code{truthX}
#'                \code{truthCs} true distribution in mode 1
#'                \code{truthDs} true distribution in mode 2
#'                \code{truthEs} true distribution in mode 3
#'                \code{mus} 
#'                \code{binaryX}
#' 
get.data = function(n,p,q,k,r,l,error=3,sort=TRUE,sparse.percent=0,center=FALSE,mumin = -3, mumax = 3){
  mus = runif(k*r*l,mumin,mumax)#take the mean of k*r*l biclusters/cubes
  if(sparse.percent!=0) mus[1:floor(k*r*l*sparse.percent)] = 0
  mus = array(mus,c(k,r,l))
  if (k!=1) truthCs = ReNumber(sample(1:k,n,replace=TRUE)) else truthCs = rep(1,n)
  if (r!=1) truthDs = ReNumber(sample(1:r,p,replace=TRUE)) else truthDs = rep(1,p)
  if (l!=1) truthEs = ReNumber(sample(1:l,q,replace=TRUE)) else truthEs = rep(1,q)
  
  ##### added by Miaoyan
  if(sort==TRUE){
      truthCs=sort(truthCs)
      truthDs=sort(truthDs)
      truthEs=sort(truthEs)
  }
  ######
  
  x = array(rnorm(n*p*q,mean=0,sd=error),dim = c(n,p,q))
  truthX = array(rep(0,n*p*q),c(n,p,q))
  for(i in 1:max(truthCs)){
    for(j in 1:max(truthDs)){
      for(m in 1:max(truthEs)){
        x[truthCs==i, truthDs==j, truthEs==m] = x[truthCs==i, truthDs==j, truthEs==m] + mus[i,j,m]
        truthX[truthCs==i, truthDs==j, truthEs==m] =  mus[i,j,m]
      }
    }
  }
  if (center == TRUE) x = x - mean(x)
  
  #### removed by Miaoyan
  #if(sort==TRUE){
  # truthX = reorderClusters(truthX,truthCs,truthDs,truthEs)$x
  # neworder = reorderClusters(x,truthCs,truthDs,truthEs)
  # binaryX = (truthX!=0)*1
  # result = list("x"=neworder$x,"truthX"=truthX,"truthCs"=neworder$Cs,"truthDs"=neworder$Ds,"truthEs"=neworder$Es,"mus"=mus,"binaryX"=binaryX)
  #} else {
    binaryX = (truthX!=0)*1
    result = list("x"=x,"truthX"=truthX,"truthCs"=truthCs,"truthDs"=truthDs,"truthEs"=truthEs,"mus"=mus,"binaryX"=binaryX)
   #}
  
  return(result)
}
