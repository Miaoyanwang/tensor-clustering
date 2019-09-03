#' Solve supervised binary tucker decomposition with covariate on two modes
#'
#' This is main updating function for the update of factor matrces and  core tensor
#'             \eqn{B = G \times_1 W_1 \times_2 W_2 \times_3 C}
#' @param tsr    the data vector
#' @param X_covar1    the covariate on first mode, if you do not have covariate on mode 1 please leave it
#' @param X_covar2    the covariate on second mode, if you do not have covariate on mode 2 please leave it
#' @param core_shape  the tucker rank of the given tensor
#' @param Nsim        simulation times
#' @param linear      whether use linear regression to get initialization, only used when covariate is identity matrix
#' @param cons        the constrain method, "non" for without constrain, "vanilla" for global scale down once at each iteration,
#'
#'                    "penalty" for add log-barrier penalty to object function.
#' @param lambda      penalty coefficient when use "penalty" constrain.
#' @param alpha       max norm consrain on ground truth
#' @param solver      solver for solving object function when using "penalty" constrain, see "details"
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
#' @details    When you choose constrain to be "penalty". It actually add log-barrier as regularizer to
#'            general object function (negative log-likelihood). It use solver in function "optim" to
#'            solve the objective function. You can choose any solver contained in function "optim".
#'
#' @export





update_binary = function(tsr, X_covar1 = NULL, X_covar2 = NULL, core_shape, Nsim, linear = TRUE,
                         cons = 'vanilla', lambda = 1, alpha = 1, solver = NULL){

  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  p_1 = dim(X_covar1)[2] ; p_2 = dim(X_covar2)[2]
  if(is.null(X_covar1)) X_covar1 = diag(d1)
  if(is.null(X_covar2)) X_covar2 = diag(d2)

  ## get initialization
  C_ts = glm_two_mat(tsr@data, X_covar1, t(X_covar2), ini = TRUE,linear=linear) ## add the linear model option for initilization

  C_ts = as.tensor(C_ts)
  tckr = tucker(C_ts, ranks = core_shape)
  W1 = tckr$U[[1]] ; W2 = tckr$U[[2]] ; C = tckr$U[[3]]
  G = tckr$Z
  A = X_covar1%*%W1
  B = X_covar2%*%W2

  lglk = rep(0,4*Nsim)

  for(n in 1:Nsim){

    ###### update W1
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data

    mod_re = glm_two(Y = Y_1, X1 = X_covar1, X2 = G_BC1, start = W1)
    W1 = mod_re[[1]]
    lglk[4*n-3] = mod_re[[2]]

    ## orthogonal W1*
    C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))
    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    C = tuk$U[[3]]
    print("W1 Done------------------")

    ##### calculate A
    A = X_covar1%*%W1

    ##### calculate B
    B = X_covar2%*%W2

    ##### update W2
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data

    mod_re = glm_two(Y_2, X_covar2, G_AC2, start = W2)
    W2 = mod_re[[1]]
    lglk[4*n-2] = mod_re[[2]]

    ## orthogonal W2*
    C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))
    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    C = tuk$U[[3]]
    print("W2 Done------------------")

    ##### calculate A
    A = X_covar1%*%W1

    ##### calculate B
    B = X_covar2%*%W2

    ###### update C
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data

    re = glm_mat(t(Y_3),start = t(C),t(G_AB3))

    C = t(re[[1]])
    lglk[4*n - 1] = re[[2]]



    ###  then we apply out constrain
    ######---- differnent version of contrains
    U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data

    if(cons == 'non'){C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))}
    else if(max(abs(U)) <= alpha){C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))}
    else if(cons == 'vanilla'){
      U = U/max(abs(U))*alpha
      C_ts = glm_two_mat(U, X_covar1, t(X_covar2), ini = TRUE,linear=linear, lm = TRUE) ## add the linear model option for initilization
      C_ts = as.tensor(C_ts)
      print("Violate constrain ------------------")
    }
    else{
      U = U/max(abs(U))*(alpha*0.01)
      C_ts = glm_two_mat(U, X_covar1, t(X_covar2), ini = TRUE,linear=linear, lm = TRUE) ## add the linear model option for initilization
      C_ts = as.tensor(C_ts)
      print("Violate constrain ------------------")
    }


    ## orthogonal C*

    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    C = tuk$U[[3]]
    print("C Done------------------")

    ###### added by miaoyan ######
    ##### calculate A
    A = X_covar1%*%W1

    ##### calculate B
    B = X_covar2%*%W2
    ############################

    ##### update G
    M_long = kronecker_list(list(C,B,A)) ## form M_long

    if(cons == 'penalty'){
      mod_re = optim(par = as.vector(G@data),loss,loss_gr,y = as.vector(tsr@data),
                     X = M_long, lambda = lambda, alpha = alpha, method = solver)
      coe = mod_re$par
      G = as.tensor(array(data = coe,dim = core_shape))
      lglk[4*n] = -mod_re$value
    }
    else {
      mod_re = glm_modify(as.vector(tsr@data), M_long, as.vector(G@data))
      coe = mod_re[[1]]
      G = as.tensor(array(data = coe,dim = core_shape))
      lglk[4*n] = mod_re[[2]]
    }


    print("G Done------------------")
    print(n)

    print(paste(n,"-th  iteration ---- when dimension is ",d1,"-- rank is ",r1," -----------------"))
    #if(lglk[4*n] - lglk[4*n-1] <= 0.00005) break
    if(abs(lglk[4*n-1] - lglk[4*n-2]) <= 0.0005) break



  }
  return(list(W1 = W1,W2 = W2,C = C,G = G,lglk = lglk))
}
