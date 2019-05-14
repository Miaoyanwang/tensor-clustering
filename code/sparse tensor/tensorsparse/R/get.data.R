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
#' @param sparse.percent the proportion of 0 mean in normal case; the proportion of 0 of each vector in multiplicative case
#' @param multiplicative if !=0, then produce overlapping multiplicative data, the number refers to the number of overlapping
#' @param center if center sum(x) = 0
#' @param seed default is NULL, otherwise would set seedd to corresponding point
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
get.data = function(n,p,q,k=NULL,r=NULL,l=NULL,error=3,sort=TRUE,sparse.percent=0,multiplicative=0,center=FALSE,seed=NULL,mumin = -3, mumax = 3){
  if(!is.null(seed)) set.seed(seed)
  if(multiplicative == 0){
    mus = runif(k*r*l,mumin,mumax)#take the mean of k*r*l biclusters/cubes
    if(sparse.percent!=0) mus[1:floor(k*r*l*sparse.percent)] = 0
    mus = array(mus,c(k,r,l))
  
  if(is.null(mus)) stop("multiplicative must be a positive integer!") 
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
  } else {
    n = as.integer(n)
    m1 = list();m2 = list();m3 = list()
    mus = array(rep(0,n*p*q),c(n,p,q))
    l1 = sample(1:k,n,replace = TRUE); l2 = sample(1:r,p,replace = TRUE); l3 = sample(1:l,q,replace = TRUE)
    s1 = sample(1:k,floor(k*sparse.percent)); s2 = sample(1:r,floor(r*sparse.percent)); s3 = sample(1:l,floor(l*sparse.percent))
    for (i in 1:multiplicative){

        r1 = runif(k,mumin,mumax);r2 = runif(r,mumin,mumax);r3 = runif(l,mumin,mumax)
        r1[s1] = 0; r2[s2] = 0; r3[s3] = 0
        m1 = r1[l1] ; m2 = r2[l2] ; m3 = r3[l3]
        m1t = new("Tensor",3L,c(n,1L,1L),data=m1)
        m2t = matrix(m2,ncol=1)
        m3t = matrix(m3,ncol=1)

        mus = mus + ttm(ttm(m1t,m2t,2),m3t,3)@data
    }
    truthX = mus
    noise = array(rnorm(n*p*q,0,error),dim = c(n,p,q))
    x = truthX + noise
    if (center == TRUE) x = x - mean(x)
    binaryX = (truthX!=0)*1
    return(list(x = x, truthX = truthX, binaryX = binaryX, truthCs = l1, truthDs = l2, truthEs = l3))
  }
}
