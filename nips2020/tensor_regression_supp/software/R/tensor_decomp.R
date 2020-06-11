###Tensor Decompositions

#'(Truncated-)Higher-order SVD
#'
#'Higher-order SVD of a K-Tensor. Write the K-Tensor as a (m-mode) product of a core Tensor (possibly smaller modes) and K orthogonal factor matrices. Truncations can be specified via \code{ranks} (making them smaller than the original modes of the K-Tensor will result in a truncation). For the mathematical details on HOSVD, consult Lathauwer et. al. (2000).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors.
#'@name hosvd
#'@rdname hosvd
#'@aliases hosvd
#'@param tnsr Tensor with K modes
#'@param ranks a vector of desired modes in the output core tensor, default is \code{tnsr@@modes}
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{core tensor with modes speficied by \code{ranks}}
#'\item{\code{U}}{a list of orthogonal matrices, one for each mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)} - if there was no truncation, then this is on the order of mach_eps * fnorm. }
#'}
#'@seealso \code{\link{tucker}}
#'@references L. Lathauwer, B.Moor, J. Vandewalle, "A multilinear singular value decomposition". Journal of Matrix Analysis and Applications 2000, Vol. 21, No. 4, pp. 1253–1278.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(6,7,8))
#'hosvdD <-hosvd(tnsr)
#'hosvdD$fnorm_resid
#'hosvdD2 <-hosvd(tnsr,ranks=c(3,3,4))
#'hosvdD2$fnorm_resid
hosvd <- function(tnsr,ranks=NULL){
  stopifnot(is(tnsr,"Tensor"))
  if (sum(ranks<=0)!=0) stop("ranks must be positive")
  if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  num_modes <- tnsr@num_modes
  #no truncation if ranks not provided
  if(is.null(ranks)){
    ranks <- tnsr@modes
  }else{
    if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
  }
  #progress bar
  pb <- txtProgressBar(min=0,max=num_modes,style=3)
  #loops through and performs SVD on mode-m matricization of tnsr
  U_list <- vector("list",num_modes)
  for(m in 1:num_modes){
    #temp_mat <- rs_unfold(tnsr,m=m)@data
    temp_mat = unfold(tnsr, row_idx = m, col_idx = (1:num_modes)[-m])@data  ## add by Zhuoyan, fixed bugs on hosvd
    U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
    setTxtProgressBar(pb,m)
  }
  close(pb)
  #computes the core tensor
  Z <- ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
  est <- ttl(Z,U_list,ms=1:num_modes)
  #resid <- fnorm(est-tnsr)
  resid = sqrt(sum((est@data-tnsr@data)^2))
  #put together the return list, and returns
  list(Z=Z,U=U_list,est=est,fnorm_resid=resid)	
}

#'Tucker Decomposition
#'
#The Tucker decomposition of a tensor. Approximates a K-Tensor using a n-mode product of a core tensor (with modes specified by \code{ranks}) with orthogonal factor matrices. If there is no truncation in one of the modes, then this is the same as the MPCA, \code{\link{mpca}}. If there is no truncation in all the modes (i.e. \code{ranks = tnsr@@modes}), then this is the same as the HOSVD, \code{\link{hosvd}}. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the Tucker decomposition, consult Kolda and Bader (2009).
#'The Tucker decomposition of a tensor. Approximates a K-Tensor using a n-mode product of a core tensor (with modes specified by \code{ranks}) with orthogonal factor matrices. If there is no truncation in all the modes (i.e. \code{ranks = tnsr@@modes}), then this is the same as the HOSVD, \code{\link{hosvd}}. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the Tucker decomposition, consult Kolda and Bader (2009).

#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure also known as Higher-Order Orthogonal Iteration (HOOI). Intialized using a (Truncated-)HOSVD. A progress bar is included to help monitor operations on large tensors.
#'@name tucker
#'@rdname tucker
#'@aliases tucker
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the modes of the output core Tensor
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{the core tensor, with modes specified by \code{ranks}}
#'\item{\code{U}}{a list of orthgonal factor matrices - one for each mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#@seealso \code{\link{hosvd}}, \code{\link{mpca}}
#'@seealso \code{\link{hosvd}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009, Vol. 51, No. 3 (September 2009), pp. 455-500. URL: https://www.jstor.org/stable/25662308
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(4,4,4,4))
#'tuckerD <- tucker(tnsr,ranks=c(2,2,2,2))
#'tuckerD$conv 
#'tuckerD$norm_percent
#'plot(tuckerD$all_resids)
tucker <- function(tnsr,ranks=NULL,max_iter=25,tol=1e-5){
  stopifnot(is(tnsr,"Tensor"))
  if(is.null(ranks)) stop("ranks must be specified")
  if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks<=0)!=0) stop("ranks must be positive")
  if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  #initialization via truncated hosvd
  num_modes <- tnsr@num_modes
  U_list <- vector("list",num_modes)
  for(m in 1:num_modes){
    temp_mat <- rs_unfold(tnsr,m=m)@data
    U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
  }
  tnsr_norm <- sqrt(sum(tnsr@data*tnsr@data))
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z,U_list){
    est <- ttl(Z,U_list,ms=1:num_modes)
    curr_resid <- sqrt(sum(((tnsr - est)@data)^2))
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{return(FALSE)}
  }
  #progress bar
  pb <- txtProgressBar(min=0,max=max_iter,style=3)
  #main loop (until convergence or max_iter)
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)	
    modes <- tnsr@modes
    modes_seq <- 1:num_modes
    for(m in modes_seq){
      #core Z minus mode m
      X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
      #truncated SVD of X
      #U_list[[m]] <- (svd(rs_unfold(X,m=m)@data,nu=ranks[m],nv=prod(modes[-m]))$u)[,1:ranks[m]]
      U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
    }
    #compute core tensor Z
    Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)
    
    #checks convergence
    if(CHECK_CONV(Z, U_list)){
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)	
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  close(pb)
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
  est <- ttl(Z,U_list,ms=1:num_modes)
  invisible(list(Z=Z, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}


