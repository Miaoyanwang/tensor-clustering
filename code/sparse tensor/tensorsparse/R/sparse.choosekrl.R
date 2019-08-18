#' Perform tuning parameter (\eqn{d_1}, \eqn{d_2} and \eqn{d_3}) selection for sparse tensor clustering via cross validation
#' 
#' Select the best \eqn{d_1}, \eqn{d_2} and \eqn{d_3} to perform clustering. A range of values of d[1], d[2], d[3] is usually considered.
#' @param x a three-dimensional array
#' @param k the range of \eqn{d_1}: a vector, the possible clusters numbers of mode 1
#' @param r the range of \eqn{d_2}: a vector, the possible clusters numbers of mode 2
#' @param l the range of \eqn{d_3}: a vector, the possible clusters numbers of mode 3
#' @param lambda a numeric value. The coefficient of the regularization term.
#' @param percent a numeric value between 0 and 1
#' @param trace  trace logic value. If true, it would print the iteration situation.
#' @param nstart positive interger. The same as the "nstart" in kmeans().
#' @param sim.times the same as label2(): the times of calling classify2() with different seeds.
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
#' @return a list   
#' \code{estimated_krl} a 1*3 matrix which is the estimated c(d_1,d_2,d_3).  
#' 
#'                \code{results.se} the standard error of each possible combination of the cluster numbers.  
#'                
#'                \code{results.mean} the mean residual of each possible combination of the cluster numbers.  
#'                
#' 
#' @export
sparse.choosekrl = function (x,k,r,l,lambda=0,percent=0.2,trace=FALSE,nstart=20,sim.times=1,method="L0") {
  #x=test;l=range.l;lambda=0;percent=0.2;trace=TRUE
  #k=2:4;r=2:4;l=2:4
  if ((1%%percent) != 0) 
    stop("1 must be divisible by the specified percentage")
  if (percent <= 0) 
    stop("percentage cannot be less than or equal to 0")
  if (percent >= 1) 
    stop("percentage cannot be larger or equal to 1")
  
  #Returns suitably lagged and iterated differences.
  if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0 || sum(diff(l) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
  n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
  miss <- sample(1:(n*p*q), n*p*q, replace = FALSE)
  numberoftimes <- 1/percent
  allresults <- array(NA, dim = c(numberoftimes, length(k), 
                                  length(r), length(l)))
  Cs.init <- matrix(NA, nrow = dim(x)[1], ncol = length(k))
  #put the kmeans results into columns
  for (i in 1:length(k)) {
    Cs.init[, i] <- kmeans(tensor.unfold(x), k[i], nstart = nstart)$cluster
  }
  
  Ds.init <- matrix(NA, nrow = dim(x)[2], ncol = length(r))
  for (j in 1:length(r)) {
    Ds.init[, j] <- kmeans(tensor.unfold(x,2), r[j], nstart = nstart)$cluster
  }
  
  Es.init <- matrix(NA, nrow = dim(x)[3], ncol = length(l))
  for (j in 1:length(l)) {
    Es.init[, j] <- kmeans(tensor.unfold(x,3), l[j], nstart = nstart)$cluster
  }
  
  for (i in 1:numberoftimes) {#i = 1 
    if (trace == TRUE) 
      cat("Iteration", i, fill = TRUE)
    xmiss <- x
    missing <- miss[1:round(n*p*q*percent)]
    xmiss[missing] <- NA
    xmiss[missing] <- mean(xmiss, na.rm = TRUE)
    
    
    krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
                   rep(1:length(r),times=length(k)*length(l)),
                   rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
                 nrow=3)
    
    
    if (.Platform$OS.type == "windows") {
      res = apply(krl,MARGIN=2,label_for_krl,k,r,l,Cs.init,Ds.init,Es.init,sim.times=sim.times,lambda=lambda,xmiss=xmiss,method=method)
    } else {
      krl_list = as.list(as.data.frame(krl))
      res = mclapply(krl_list, label_for_krl,k,r,l,Cs.init,Ds.init,Es.init,sim.times=sim.times,lambda=lambda,xmiss=xmiss,method=method,mc.cores=n.cores)
    }
    
    
    for (a in 1:dim(krl)[2]){
      #print(krl[,a])
      #print(sum((x[missing]-res[[a]][[1]][missing])^2))
      allresults[i,krl[,a][1],krl[,a][2],krl[,a][3]] = sum((x[missing]-res[[a]][[1]][missing])^2)
    }
    miss <- miss[-1:-(dim(x)[1] * dim(x)[2] * dim(x)[3]/numberoftimes)]
  }
  results.se <- apply(allresults, MARGIN=c(2, 3, 4), sd)/sqrt(numberoftimes)
  results.mean <- apply(allresults, c(2, 3, 4), mean)
  #comparing every mean with the mean with higher k,r,l
  IndicatorArray <- 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
                                                              1), 1:(length(l) - 1)] <= results.mean[2:length(k), 2:length(r), 2:length(l)] + results.se[2:length(k), 2:length(r), 2:length(l)])
  if (max(IndicatorArray) == 0) {
    warning("The krl has reached the upper boundary. Please enlarger the range.")
    return(list(estimated_krl = c(max(k),max(r),max(l))))}
  
  #outer(outer(matrix(1:4,2),matrix(5:8,2)),matrix(1:4,2))
  ModeIndex <- outer(outer(k[-length(k)], r[-length(r)],  "*"), l[-length(l)], "*")
  smallestIndicatorTrue <- min(ModeIndex[IndicatorArray == TRUE])
  out <- which(IndicatorArray == TRUE & ModeIndex == 
                 smallestIndicatorTrue, arr.ind = TRUE)
  out <- array(c(k[out[,1]], r[out[,2]], l[out[,3]]), dim=dim(out))
  tempmode1 <- NULL
  tempmode2 <- NULL
  tempmode3 <- NULL
  for (i in 1:length(k)) {
    tempmode1 <- c(tempmode1, paste("K = ", k[i], sep = ""))
  }
  for (i in 1:length(r)) {
    tempmode2 <- c(tempmode2, paste("R = ", r[i], sep = ""))
  }
  for (i in 1:length(l)) {
    tempmode3 <- c(tempmode3, paste("L = ", l[i], sep = ""))
  }
  
  
  dimnames(results.se)[[1]] = tempmode1
  dimnames(results.se)[[2]] = tempmode2
  dimnames(results.se)[[3]] = tempmode3
  dimnames(results.mean)[[1]] = tempmode1
  dimnames(results.mean)[[2]] = tempmode2
  dimnames(results.mean)[[3]] = tempmode3
  return(list(estimated_krl = out, results.se = results.se, 
              results.mean = results.mean))
}
