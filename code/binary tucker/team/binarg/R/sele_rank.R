#' Select rank of the coefficient tensor given binary tensor observation
#'
#' This is the function for the selecting the rank of coefficient tensor i.e. size of the core tensor
#'             \eqn{B = G \times_1 W_1 \times_2 W_2 \times_3 C}
#'             by BIC criterion.
#' @param tsr    the data vector
#' @param X_covar1    the covariate on first mode, unsupervised version please see "details"
#' @param X_covar2    the covariate on second mode, unsupervised version please see "details"
#' @param rank        a vector containing the several rank candidates
#' @param Nsim        simulation times in each update
#' @param linear      whether use linear regression to get initialization, only used when covariate is identity matrix
#' @param cons        the constrain method, "non" for without constrain, "vanilla" for global scale down once at each iteration,
#'
#'                    "penalty" for add log-barrier penalty to object function.
#' @return     a list containing following components:
#'
#'                    \code{W1} coefficient matrix on mode-1
#'
#'                    \code{W2} coefficient matrix on mode-2
#'
#'                    \code{C} factor matrix on mode-3
#'
#'                    \code{G} core tensor
#'
#'                    \code{lglk} the lig-likelihood
#'
#' @details    For selecting rank, recommend use non-constrain version updating.
#'
#'             When you select the rank of binary in unsupervised verison, please leave two covariate to blank,
#'             or you can set them as identity matrix.
#'
#'
#' @export




sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, rank = c(3,5), Nsim,linear = FALSE, cons = 'non'){
    rank = expand.grid(rank,rank,rank)
    colnames(rank) = NULL
    rank = as.matrix(rank)
    whole_shape = dim(tsr)
    rank = lapply(1:dim(rank)[1], function(x) rank[x,]) ## turn rank to a list
    if((is.null(X_covar1)) & (is.null(X_covar2)) | (all.equal(X_covar1,diag(dim(tsr)[1])) & all.equal(X_covar2,diag(dim(tsr)[2])))){
        upp = lapply(rank, FUN= update_binary_un,tsr = tsr, Nsim = Nsim, cons = cons)
        
    }
    else {
        upp = lapply(rank, FUN= update_binary,tsr = tsr,X_covar1 = X_covar1, 
        X_covar2 = X_covar2, Nsim = Nsim, linear = linear, cons = cons)
    }
    log_lik = unlist(lapply(seq(length(upp)), function(x) max(upp[[x]]$lglk)))
    if(is.null(X_covar1)) X_covar1 = diag(d1)
    if(is.null(X_covar2)) X_covar2 = diag(d2)
    
    BIC = unlist(lapply(seq(length(rank)), 
    function(x) (prod(rank[[x]]) + sum((c(dim(X_covar1)[2],dim(X_covar2)[2],
    whole_shape[3])-1)*
    rank[[x]])) * log(prod(whole_shape))))
    BIC = -2*log_lik + BIC
    # BIC = rep(0,length(rank))
    # for(i in 1:length(rank)){
    #     BIC[i] =  (prod(rank[[i]]) + sum((c(dim(X_covar1)[2],dim(X_covar2)[2],whole_shape[3])-1)*
    #                                  rank[[i]]))*log(prod(whole_shape))
    #   }
    res=cbind(data.frame(t(sapply(rank,c))),BIC)
    return(list(rank = rank[[which(BIC == min(BIC))]],res=res))
}




