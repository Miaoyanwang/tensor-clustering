#' Generalized Tensor Regression
#'
#' Tensor-response regression given covariates on multiple modes. Main function in the package. The function takes a response tensor, multiple covariate matrices, and a desired Tucker rank as input. The output is a constrained
#' MLE for the coefficient tensor.
#'
#' @param tsr    response tensor with 3 modes
#' @param X_covar1    covariate on first mode
#' @param X_covar2    covariate on second mode
#' @param X_covar3    covariate on third mode
#' @param core_shape  the Tucker rank of the regression coefficients
#' @param Nsim        max number of iterations if update does not convergence
#' @param cons        the constraint method, "non" for without constraint, "vanilla" for global scale down at each iteration,
#'
#'                    "penalty" for adding log-barrier penalty to object function
#' @param lambda      penalty coefficient for "penalty" constraint
#' @param alpha       max norm constraint on linear predictor
#' @param solver      solver for solving object function when using "penalty" constraint, see "details"
#' @param dist        distribution of the response tensor, see "details"
#' @return     a list containing the following:
#'
#'                  \code{W} {a list of orthogonal coefficient matrices - one for each mode, with the number of columns given by \code{core_shape}}
#'
#'                  \code{G}  {core tensor, with the size specified by \code{core_shape}}
#'
#'                  \code{C_ts}  {coefficient tensor, Tucker product of \code{G},\code{A},\code{B},\code{C}}
#'
#'                  \code{U} {linear predictor,i.e. Tucker product of \code{C_ts},\code{X_covar1},\code{X_covar2},\code{X_covar3}}
#'
#'                  \code{lglk} {loglikelihood at convergence}
#'
#'                  \code{sigma} {estimated error variance (for Gaussian tensor) or dispersion parameter (for Bernoulli and Poisson tensors)}
#'
#'                  \code{violate} {a vector listing whether each iteration violates the max norm constraint on the linear predictor, \code{1} indicates violation}
#'
#'
#'
#' @details   Constraint \code{penalty} adds log-barrier regularizer to
#'            general object function (negative log-likelihood). The main function uses solver in function "optim" to
#'            solve the objective function. The "solver" passes to the argument "method" in function "optim".
#'
#'            \code{dist} specifies three distributions of response tensor: binary, poisson and normal distribution.
#'
#' @export

tensor_regress = function(tsr,X_covar1 = NULL, X_covar2 = NULL,X_covar3 = NULL, core_shape, Nsim=20, cons = c("non","vanilla","penalty"), lambda = 0.1, alpha = 1, solver ="CG",dist = c("binary", "poisson","normal")){


  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]

  ###  check whether unsupervised on each mode
  un_m1 = FALSE ; un_m2 = FALSE ; un_m3 = FALSE
  if(is.null(X_covar1)|(identical(X_covar1,diag(d1)))) {X_covar1 = diag(d1) ; un_m1 = TRUE}
  if(is.null(X_covar2)|(identical(X_covar2,diag(d2)))) {X_covar2 = diag(d2) ; un_m2 = TRUE}
  if(is.null(X_covar3)|(identical(X_covar3,diag(d3)))) {X_covar3 = diag(d3) ; un_m3 = TRUE}
  p1 = dim(X_covar1)[2] ; p2 = dim(X_covar2)[2] ; p3 = dim(X_covar3)[2]


  if(dist=="binary"){
    tsr.transform=as.tensor(2*tsr@data-1)
  }else if(dist=="poisson"){
    tsr.transform=as.tensor(log(tsr@data+0.1))###?? new initilization
  }else if (dist=="normal"){
    tsr.transform=tsr
  }

  C_ts=ttl(tsr.transform,list(ginv(X_covar1),ginv(X_covar2),ginv(X_covar3)),ms=c(1,2,3))

  tckr = tucker(C_ts, ranks = core_shape)
  W1 = tckr$U[[1]] ; W2 = tckr$U[[2]] ; W3 = tckr$U[[3]] ## tucker factors
  G = tckr$Z
  A = X_covar1%*%W1
  B = X_covar2%*%W2
  C = X_covar3%*%W3

  core=update_core(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist)
  G=core$G
  lglk=core$lglk
  violate=core$violate

  for(n in 1:Nsim){
    ## parameter from previous step

    W10 = W1 ; W20 = W2 ; W30 = W3 ; G0=G; A0=A;B0=B;C0=C;lglk0=tail(lglk,1);

    ###### update W1
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data

    if(un_m1) {re = glm_mat(t(Y_1),t(G_BC1),dist=dist) ## no covariate
    } else {re = glm_two(Y = Y_1, X1 = X_covar1, X2 = G_BC1, dist=dist)}

    if(dim(re[[1]])[1]==1) W1=t(re[[1]]) else W1 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])

    ## orthogonal W1*
    qr_res=qr(W1)
    W1=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),1)
    #print("W1 Done------------------")

    ##### calculate A
    A = X_covar1%*%W1;
    ##### update W2
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data

    if(un_m2) {re = glm_mat(t(Y_2),t(G_AC2),dist=dist)
    } else {re = glm_two(Y_2, X_covar2, G_AC2, dist=dist)}

    if(dim(re[[1]])[1]==1) W2=t(re[[1]]) else W2 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])

    ## orthogonal W2*
    qr_res=qr(W2)
    W2=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),2)
    #print("W2 Done------------------")

    ##### calculate B
    B = X_covar2%*%W2;


    ###### update W3
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data

    if(un_m3) {re = glm_mat(t(Y_3),t(G_AB3),dist=dist)
    } else {re = glm_two(Y_3, X_covar3, G_AB3,dist=dist)}

    if(dim(re[[1]])[1]==1) W3=t(re[[1]]) else W3 = as.matrix(re[[1]])

    lglk = c(lglk,re[[2]])

    ## orthogonal W3*
    qr_res=qr(W3)
    W3=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),3)
    #print("W3 Done------------------")

    ##### calculate C
    C = X_covar3%*%W3;

    #########-----------------------------------------------
    ###  obtain core tensor under constraint
    core=update_core(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist)
    G=core$G
    lglk=c(lglk,core$lglk)
    violate=c(violate,core$violate)

    #print("G Done------------------")

    print(paste(n,"-th  iteration -- when dimension is",d1,d2,d3,"- rank is ",r1,r2,r3," -----------------"))
    #print(paste(n,"-th  iteration"))

    if ((tail(lglk,1)-lglk0)/abs(lglk0)<= 0.0001 & tail(lglk,1)>= lglk0 ){
      print(paste(n,"-th iteration: convergence"))
      break
    } else if (tail(lglk,1)-lglk0 < 0) {
      W1 = W10 ; W2 = W20 ; W3 = W30; G=G0; lglk=lglk[-c((length(lglk)-3):length(lglk))];
      A=A0;B=B0;C=C0;
      break
    }

  }

  U=ttl(G,list(A,B,C),ms = c(1,2,3))@data

  sigma_est=mean((tsr@data-U_to_mean(U,dist))^2)

  return(list(W = list(W1 = W1,W2 = W2,W3 = W3),G = G,U=U, C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data,lglk = lglk, sigma=sigma_est,violate = violate))
}
