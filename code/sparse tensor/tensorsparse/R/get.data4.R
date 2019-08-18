#' Generate a random 4th-order tensor
#' 
#' Generate a random 4th-order tensor.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param s the dimension of mode 4
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param m \eqn{d_4}: the clusters number of mode 4
#' @param error positive numeric value: the noise in rnorm()
#' @param sort if TRUE, the data belongs to the same clusters would assemble together after sorting
#' @param sparse.percent the proportion of 0 mean in normal case; the proportion of 0 of each vector in multiplicative case
#' @param center logic value that indicates whether run "x = x-mean(x)" before performing clustering.
#' @param seed default is NULL, otherwise would set seed to corresponding point
#' @param mumin numeric value. The lower bound of mu when sampling mu
#' @param mumax numeric value. The lower bound of mu when sampling mu
#' @return a list \code{x} the tensor   
#' 
#'                \code{truthX} the tensor before adding the noise   
#'                
#'                \code{truthCs} true distribution in mode 1   
#'                
#'                \code{truthDs} true distribution in mode 2   
#'                
#'                \code{truthEs} true distribution in mode 3   
#'                
#'                \code{mus} the mean signal of all clusters   
#'                
#'                \code{binaryX} the 0-1 tensor (0:the mean signal = 0; 1:the mean signal != 0)   
#'                
#' @export
#' 
#' 
get.data4 = function(n,p,q,s,k=NULL,r=NULL,l=NULL,m=NULL,error=3,sort=TRUE,sparse.percent=0,center=FALSE,seed=NULL,mumin = -3, mumax = 3){
  if(!is.null(seed)) set.seed(seed)
  #print(m)
  mus = runif(k*r*l*m,mumin,mumax)#take the mean of k*r*l biclusters/cubes
  if(sparse.percent!=0) mus[sample(k*r*l*m,floor(k*r*l*m*sparse.percent),replace=F)]= 0
  mus = array(mus,c(k,r,l,m))
  
  if(is.null(mus)) stop("multiplicative must be a positive integer!") 
  if (k!=1) truthCs = ReNumber(sample(1:k,n,replace=TRUE)) else truthCs = rep(1,n)
  if (r!=1) truthDs = ReNumber(sample(1:r,p,replace=TRUE)) else truthDs = rep(1,p)
  if (l!=1) truthEs = ReNumber(sample(1:l,q,replace=TRUE)) else truthEs = rep(1,q)
  if (m!=1) truthFs = ReNumber(sample(1:m,s,replace=TRUE)) else truthFs = rep(1,s)
  
  ##### added
  if(sort==TRUE){
    truthCs=sort(truthCs)
    truthDs=sort(truthDs)
    truthEs=sort(truthEs)
    truthFs=sort(truthFs)
  }
  ######
  
  x = array(rnorm(n*p*q*s,mean=0,sd=error),dim = c(n,p,q,s))
  truthX = array(rep(0,n*p*q*s),c(n,p,q,s))
  for(i1 in 1:max(truthCs)){
    for(i2 in 1:max(truthDs)){
      for(i3 in 1:max(truthEs)){
        for(i4 in 1:max(truthFs)){
          x[truthCs==i1, truthDs==i2, truthEs==i3, truthFs==i4] = x[truthCs==i1, truthDs==i2, truthEs==i3, truthFs==i4] + mus[i1,i2,i3,i4]
          truthX[truthCs==i1, truthDs==i2, truthEs==i3, truthFs==i4] =  mus[i1,i2,i3,i4]
        }
      }
    }
  }
  if (center == TRUE) x = x - mean(x)
  binaryX = (truthX!=0)*1
  result = list("x"=x,"truthX"=truthX,"truthCs"=truthCs,"truthDs"=truthDs,"truthEs"=truthEs,"truthFs"=truthFs,"mus"=mus,"binaryX"=binaryX)
  #}
  return(result)
}
