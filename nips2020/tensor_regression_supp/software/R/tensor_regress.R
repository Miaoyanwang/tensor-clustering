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
#' @param cons        the constraint method, "non" for without constraint, "vanilla" for global scale down at each iteration, "penalty" for adding log-barrier penalty to object function
#' @param lambda      penalty coefficient for "penalty" constraint
#' @param alpha       max norm constraint on linear predictor
#' @param solver      solver for solving object function when using "penalty" constraint, see "details"
#' @param dist        distribution of the response tensor, see "details"
#' @return     a list containing the following:
#'
#'                  \code{W} {a list of orthogonal coefficient matrices - one for each mode, with the number of columns given by \code{core_shape}}
#'
#'                  \code{G}  {an array, core tensor with the size specified by \code{core_shape}}
#'
#'                  \code{C_ts}  {an array, coefficient tensor, Tucker product of \code{G},\code{A},\code{B},\code{C}}
#'
#'                  \code{U} {linear predictor,i.e. Tucker product of \code{C_ts},\code{X_covar1},\code{X_covar2},\code{X_covar3}}
#'
#'                  \code{lglk} {a vector containing loglikelihood at convergence}
#'
#'                  \code{sigma} {a scalar, estimated error variance (for Gaussian tensor) or dispersion parameter (for Bernoulli and Poisson tensors)}
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
#' @author Z. Xu \email{zxu444@wisc.edu}  and J. Hu \email{jhu267@wisc.edu}
#' @references Z. Xu, J. Hu, and M. Wang, "Generalized Tensor Regression with Covariates on Multiple Modes". 2019. <arXiv:1910.09499>. URL: https://arxiv.org/abs/1910.09499
#' @export
#' @examples 
#' seed = 34
#' dist = 'binary'
#' data=sim_data(seed, whole_shape = c(20,20,20), core_shape=c(3,3,3),
#' p=c(5,5,5),dist=dist, dup=5, signal=4)
#' re = tensor_regress(data$tsr[[1]],data$X_covar1,data$X_covar2,data$X_covar3,
#' core_shape=c(3,3,3),Nsim=10, cons = 'non', dist = dist)

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

    message(paste(n,"-th  iteration -- when dimension is",d1,d2,d3,"- rank is ",r1,r2,r3," -----------------"))
    #print(paste(n,"-th  iteration"))

    if ((tail(lglk,1)-lglk0)/abs(lglk0)<= 0.0001 & tail(lglk,1)>= lglk0 ){
      message(paste(n,"-th iteration: convergence"))
      break
    } else if (tail(lglk,1)-lglk0 < 0) {
      W1 = W10 ; W2 = W20 ; W3 = W30; G=G0; lglk=lglk[-c((length(lglk)-3):length(lglk))];
      A=A0;B=B0;C=C0;
      break
    }

  }

  U=ttl(G,list(A,B,C),ms = c(1,2,3))@data

  sigma_est=mean((tsr@data-U_to_mean(U,dist))^2)

  return(list(W = list(W1 = W1,W2 = W2,W3 = W3),G = G@data,U=U, C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data,lglk = lglk, sigma=sigma_est,violate = violate))
}




#' Simulation of tensor regression models
#'
#' Generate response tensors with multiple covariates under different simulation models, specifically for tensors with 3 modes
#' @param seed         a random seed for generating data
#' @param whole_shape  a vector containing dimension of the tensor
#' @param core_shape   a vector containing Tucker rank of the coefficient tensor
#' @param p            a vector containing numbers of covariates on each mode, see "details"
#' @param dist         distribution of response tensor, see "details"
#' @param dup          number of simulated tensors from the same linear predictor
#' @param signal       a scalar controlling the max norm of the linear predictor
#' @param block        a vector containing boolean variables, see "details"
#' @return     a list containing the following:
#'
#' \code{tsr} {a list of simulated tensors, with the number of replicates specified by \code{dup}}
#'
#' \code{X_covar1}  {a matrix, covariate on first mode}
#'
#' \code{X_covar2}  {a matrix, covariate on second mode}
#'
#' \code{X_covar3}  {a matrix, covariate on third mode}
#'
#' \code{W} {a list of orthogonal coefficient matrices - one for each mode, with the number of columns given by \code{core_shape}}
#'
#' \code{G}  {an array, core tensor with size specified by \code{core_shape}}
#'
#' \code{C_ts}  {an array, coefficient tensor, Tucker product of \code{G},\code{A},\code{B},\code{C}}
#'
#' \code{U} {an array, linear predictor,i.e. Tucker product of \code{C_ts},\code{X_covar1},\code{X_covar2},\code{X_covar3}}
#'
#' @details    By default non-positive entry in \code{p} indicates no covariate on the corresponding mode of the tensor.
#'
#'             \code{dist} specifies three distributions of response tensor: binary, poisson or normal distribution.
#'
#'             \code{block} specifies whether the coefficient factor matrix is a membership matrix, set to \code{TRUE} when utilizing the stochastic block model
#'
#'
#' @export
#' @examples 
#' seed = 34
#' dist = 'binary'
#' data=sim_data(seed, whole_shape = c(20,20,20), core_shape=c(3,3,3),
#' p=c(5,5,5),dist=dist, dup=5, signal=4)


###----  functions for simulation
#####---- This is the function used for generating data through different distribution
#         of core tensor in  semi-supervised setting
## p is the dimension of the covaraite. p = 0 represents the case without covaraites
sim_data = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),p=c(3,3,0),dist, dup, signal,block=rep(FALSE,3)){
  
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  p1 = p[1]; p2 = p[2]; p3 = p[3];
  
  #### warning for r should larger than 0
  if(r1<=0 | r2 <= 0|r3<= 0){
    warning("the rank of coefficient tensor should be larger than 0",immediate. = T)
    return()
  }
  
  if(p1 > d1|p2 > d2|p3 > d3){
    warning("the number of covariates at each mode should be no larger than the dimension of the tensor",immediate. = T)
  }
  
  #### warning for p should larger than r
  if((p1<r1 & p1>0)|(p2<r2 & p2>0)|(p3<r3 & p3>0)){
    warning("the rank of coefficient tensor should be no larger than the number of covariates",immediate. = T)
  }
  
  
  
  ####-------- generate data
  set.seed(seed)  # 24 # 37  #  347
  X_covar1 = X_covar2 = X_covar3 = NULL
  
  if(p1<=0){
    X_covar1=diag(1,d1)
    p1=d1
  }else{
    X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 1/sqrt(d1)),d1,p1)
  }
  
  if(p2<=0){
    X_covar2=diag(1,d2)
    p2=d2
  }else{
    X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd =1/sqrt(d2)),d2,p2)
  }
  
  
  if(p3<=0){
    X_covar3=diag(1,d3)
    p3=d3
  }else{
    X_covar3 = matrix(rnorm(d3*p3,mean = 0, sd = 1/sqrt(d3)),d3,p3)
  }
  
  #### if use block, p and r should larger than 1 of smaller than 0
  if(block[1] == T & (p1 ==1|r1 == 1)){
    warning("number of groups should be larger than 1",immediate. = T)
  }
  
  if(block[2] == T & (p2 ==1|r2 == 1)){
    warning("number of groups should be larger than 1",immediate. = T)
  }
  
  if(block[3] == T & (p3 ==1|r3 == 1)){
    warning("number of groups should be larger than 1",immediate. = T)
  }
  
  
  if (block[1]==TRUE){
    b1=sort(sample(1:r1,p1,replace=TRUE))
    if(length(unique(b1) )== 1){
      warning("mode-1 degenerates to 1",immediate. = T)
    }
    W1=model.matrix(~-1+as.factor(b1))
    r1=dim(W1)[2]
  }else W1 =as.matrix(randortho(p1)[,1:r1])
  
  A = X_covar1%*%W1 ## factor matrix
  
  if (block[2]==TRUE){
    b2=sort(sample(1:r2,p2,replace=TRUE))
    if(length(unique(b2)) == 1){
      warning("mode-2 degenerates to 1",immediate. = T)
    }
    W2=model.matrix(~-1+as.factor(b2))
    r2=dim(W2)[2]
  }else W2 = as.matrix(randortho(p2)[,1:r2])
  B = X_covar2%*%W2 ## factor matrix
  
  if (block[3]==TRUE){
    b3=sort(sample(1:r3,p3,replace=TRUE))
    if(length(unique(b3)) == 1){
      warning("mode-3 degenerates to 1",immediate. = T)
    }
    W3=model.matrix(~-1+as.factor(b3))
    r3=dim(W3)[2]
  }else W3 = as.matrix(randortho(p3)[,1:r3])
  C= X_covar3%*%W3 ## factor matrix
  
  
  ### G: core tensor
  G = as.tensor(array(runif(r1*r2*r3,min=-1,max=1),dim = c(r1,r2,r3)))
  
  
  ### U: linear predictor
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  G=G/max(abs(U))*signal ## rescale subject to entrywise constraint
  U=U/max(abs(U))*signal
  
  C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data ## coefficient
  
  ### tsr:binary tensor
  if(dist=="binary"){
    tsr = lapply(seq(dup), function(x) array(rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)))),dim = c(d1,d2,d3)))}
  else if (dist=="normal"){
    tsr = lapply(seq(dup), function(x) array(rnorm(d1*d2*d3,U),dim = c(d1,d2,d3)))}#sd = 1
  else if (dist=="poisson"){
    tsr = lapply(seq(dup), function(x) array(rpois(d1*d2*d3,exp(U)),dim = c(d1,d2,d3)))}
  
  
  return(list(tsr = tsr,X_covar1 = X_covar1, X_covar2 = X_covar2,X_covar3 = X_covar3,
              W = list(W1 = W1,W2 = W2,W3 = W3), G=G@data, U=U,C_ts=C_ts))
}




#' Rank selection
#'
#' Estimate the Tucker rank of coefficient tensor based on BIC criterion. The choice of BIC
#'  aims to balance between the goodness-of-fit for the data and the degree of freedom in the population model.
#' @param tsr    response tensor with 3 modes
#' @param X_covar1    covariate on first mode
#' @param X_covar2    covariate on second mode
#' @param X_covar3    covariate on third mode
#' @param rank_range  a matrix containing rank candidates on each row
#' @param Nsim        max number of iterations if update does not convergence
#' @param cons        the constraint method, "non" for without constraint, "vanilla" for global scale down at each iteration,
#'
#'                    "penalty" for adding log-barrier penalty to object function.
#' @param lambda      penalty coefficient for "penalty" constraint
#' @param alpha       max norm constraint on linear predictor
#' @param solver      solver for solving object function when using "penalty" constraint, see "details"
#' @param dist        distribution of response tensor, see "details"
#' @return     a list containing the following:
#'
#'                    \code{rank} a vector with selected rank with minimal BIC
#'
#'                    \code{result}  a matrix containing rank candidate and its loglikelihood and BIC on each row

#' @details    For rank selection, recommend using non-constraint version.
#'             
#'            Constraint \code{penalty} adds log-barrier regularizer to
#'            general object function (negative log-likelihood). The main function uses solver in function "optim" to
#'            solve the objective function. The "solver" passes to the argument "method" in function "optim".
#'            
#'             \code{dist} specifies three distributions of response tensor: binary, poisson and normal distributions.
#'
#'
#' @export
#' @examples 
#' seed=24
#' dist='binary'
#' data=sim_data(seed, whole_shape = c(20,20,20),
#' core_shape=c(3,3,3),p=c(5,5,5),dist=dist, dup=5, signal=4)
#' rank_range = rbind(c(3,3,3),c(3,3,2),c(3,2,2),c(2,2,2),c(3,2,3)) 
#' re = sele_rank(data$tsr[[1]],data$X_covar1,data$X_covar2,data$X_covar3,
#'  rank_range = rank_range,Nsim=10,cons = 'non',dist = dist)






sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,rank_range,Nsim=10,cons = 'non', lambda = 0.1, alpha = 1, solver ='CG',dist){
  whole_shape=dim(tsr)                      
  p=rep(0,3)
  if(is.null(X_covar1)) p[1]=whole_shape[1] else p[1]=dim(X_covar1)[2]
  if(is.null(X_covar2)) p[2]=whole_shape[2] else p[2]=dim(X_covar2)[2]
  if(is.null(X_covar3)) p[3]=whole_shape[3] else p[3]=dim(X_covar3)[2]
  
  
  rank_matrix=rank_range
  rank=as.matrix(rank_range)
  
  whole_shape = dim(tsr)
  rank = lapply(1:dim(rank)[1], function(x) rank[x,]) ## turn rank to a list
  upp = lapply(rank, FUN= tensor_regress,tsr = tsr,X_covar1 = X_covar1,X_covar2 = X_covar2,X_covar3 = X_covar3, Nsim = Nsim, cons = cons,lambda = lambda, alpha = alpha, solver = solver,dist=dist)
  
  lglk= unlist(lapply(seq(length(upp)), function(x) tail(upp[[x]]$lglk,1)))
  BIC = unlist(lapply(seq(length(rank)), function(x) (prod(rank[[x]]) + sum((p-rank[[x]])*rank[[x]])) * log(prod(whole_shape))))
  BIC = -2*lglk + BIC
  rank_matrix=cbind(rank_matrix,lglk,BIC)
  
  return(list(rank = rank[[which(BIC == min(BIC))]],result=rank_matrix))
}












