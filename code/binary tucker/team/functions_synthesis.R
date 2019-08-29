library(rTensor)
library(pracma)
library(gtools)
library(rgl)
library("RColorBrewer")

#######################
######----------   functions for update

## This part used for penalty (Conjugate gradient method)
loss = function(beta,y,X,lambda,alpha){
  U = X %*% beta
  p=plogis(X %*% beta) 
  p[p==1] = p[p==1] - 0.0001
  p[p==0] = p[p==0] + 0.0001
  # if U is too large/small exp(U)/(1+exp(U)) would give value 1/0
  # that leads to the penalty function can not be intialized
  L=-y*log(p)-(1-y)*log(1-p)
  L2 =sum(L) - lambda * sum(log(1 - (U / alpha)^2))  ## object function
  #L2 =sum(L) + lambda * sum(log((U / alpha)^2))  ## object function
  return(c(L2))
}

## This part follows Sec 2.1.2 to compute Gradient of object function
loss_gr <- function(beta,y,X,lambda,alpha){
  U = X %*% beta
  p=plogis(X %*% beta)
  L_g=t(X) %*% (p - y)
  penal_g = 2 * lambda * t(X) %*% (U/(alpha^2 - U^2))
  #penal_g = 2 * lambda * t(X) %*% (1/U)
  return(c(L_g) + c(penal_g))
}

# author: prof. Wang
glm_modify=function(y,x,start){
  
  ## initial coefficent
  ini_loglik=sum(log(inv.logit((2*y-1)*(x%*%start))))
  
  ## Option 1: glm fittig with default initilization
  fit1 = suppressWarnings(glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F)))
  
  ## Option 2: glm with user specified initilization
  fit2= suppressWarnings(glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F),start=start))
  
  ## report the result whichever gives the highest likelihood
  if(max(logLik(fit1),logLik(fit2))<ini_loglik) return (list(start, ini_loglik))
  else if(logLik(fit1)>logLik(fit2)) return(list(coef(fit1), logLik(fit1)))
  else return(list(coef(fit2), logLik(fit2)))
  
}

############################################################
glm_f = function(y,start,X)  
{
  re = glm_modify(y, X, start)
  coe = re[[1]]
  lglk = re[[2]]
  return(list(coe,lglk))
}

#    form U = X*Beta
glm_mat = function(Y,start,X){
  R = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  ma = mapply(glm_f, as.data.frame(Y),as.data.frame(start),MoreArgs = list(X))
  re = t(matrix(unlist(ma), nrow = p, byrow=T))
  beta = re[1:R,]
  lglk = sum(re[R + 1,])
  return(list(beta,lglk))
}

###########---------  GLM on two modes
##---  This section follows 1.3 of [Algorithm: Semi-Supervised Binary Tensor Factorization]
## add the linear model option for initilization
##---  function version
glm_two = function(Y, X1, X2, ini = FALSE, start, linear=FALSE, lm = FALSE # used for supervised scale down
){ ## Y_size = m * n
  #X2 = t(X2)
  # logit(E(Y)) = X1 %*% coe %*% X2
  m = dim(Y)[1] ; n = dim(Y)[2]
  q1 = dim(X1)[2] ; q2 = dim(X2)[1]
  
  N_long = kronecker_list(list(t(X2),X1))
  
  if(ini == TRUE){
    coe_start = rnorm(q1*q2)
  }
  else {coe_start = as.vector(start)}
  
  if(linear==TRUE){
    mod_re = lm(10*(2*as.vector(Y)-1)~-1+N_long)
    coe = matrix(mod_re[[1]], nrow = q1, ncol = q2)
    lglk= logLik(mod_re)
    return(list(coe = coe, lglk = lglk))
    
  }
  
  if(lm==TRUE){
    mod_re = lm(as.vector(Y)~-1+N_long)
    coe = matrix(mod_re[[1]], nrow = q1, ncol = q2)
    lglk= logLik(mod_re)
    return(list(coe = coe, lglk = lglk))
    
  }
  
  mod_re = glm_modify(as.vector(Y), N_long, coe_start)
  coe = mod_re[[1]]
  coe = matrix(coe, nrow = q1, ncol = q2)
  lglk= mod_re[[2]]
  return(list(coe = coe, lglk = lglk))
}

#####---- This function is a parallel version of GLM on two modes, used for initialization
## add the linear model option for initilization
glm_two_mat = function(Y, X1, X2, ini = TRUE, start = NULL,linear, lm = FALSE){
  Yl = lapply(seq(dim(Y)[3]), function(x) Y[ , , x])
  re = lapply(Yl, glm_two, X1, X2, ini = ini,linear=linear, lm = lm)
  
  coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
  coe = array(unlist(coe), dim = c(dim(X1)[2],dim(X2)[1],dim(Y)[3]))  ## form coe
  return(coe)
}



#################  update
# usage:
# where cons = 'non', 'vanilla', 'penalty'

# solver = any solver consisted in function optim
# recommend the method require gradient(since it has higher convergence rate)

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
      U = U/max(abs(U))*(alpha-0.01)
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


###
###----  This function shows how we select rank 
##   recommend to use non-constrain verison to select rank

sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, rank = c(3,5,6,8), Nsim,
                     linear = FALSE, cons = 'non'){
  BIC = rep(0,length(rank))
  whole_shape = dim(tsr)
  for(i in 1:length(rank)){
    if((is.null(X_covar1)) & (is.null(X_covar2)) | (all.equal(X_covar1,diag(dim(tsr)[1])) & all.equal(X_covar2,diag(dim(tsr)[2])))){
      upp = update_binary_un(tsr, rep(rank[i],3), Nsim, cons = cons)
      log_Lik = max(upp$lglk)
      BIC[i] = -2*log_Lik + (rank[i]^3 + sum((whole_shape-1)*rank[i]))*log(prod(whole_shape))
    }
    else {
      upp = update_binary(tsr, X_covar1, X_covar2, rep(rank[i],3), Nsim, linear, cons = cons)
      log_Lik = max(upp$lglk)
      BIC[i] = -2*log_Lik + (rank[i]^3 + sum((c(dim(X_covar1)[2],dim(X_covar2)[2],whole_shape[3])-1)*rank[i]))*log(prod(whole_shape))
    }
    
    
  }
  return(rank = which(BIC == min(BIC)))
}

#### This function is to select the lambda in the penalty constrain version
select_lambda = function(tsr,trueU,lambda, X_covar1, X_covar2, Nsim, linear = TRUE){
  #lambda is a sequence of lambda options 
  #tsr is an array
  #trueU is an array
  #have selected the rank
  
  d1 = dim(tsr)[1]; d2 = dim(tsr)[2]; d3 = dim(tsr)[3]
  r1 = dim(trueU)[1] ; r2 = dim(trueU)[2] ; r3 = dim(trueU)[3]
  
  len = length(lambda)
  rmse = c()
  
  for (i in 1:len) {
    upp2 = update_binary(tsr, X_covar1, X_covar2, core_shape = c(r1,r2,r3), Nsim, linear = linear,
                         cons = "penalty", lambda = lambda[i], alpha = 10*max(abs(trueU)), solver = "CG")
    U_est2 = ttl(upp2$G,list(upp2$A,upp2$B,upp2$C),ms = c(1,2,3))@data
    
    rmse[i] =  sum((U_est2 - trueU)^2)/(d1*d2*d3)
  }
  
  best_lambda = lambda[order(rmse)[1]]
  return(best_lambda)
}



###----  functions for simulation
#####---- This is the function used for generating data through different distribution
#         of core tensor in unsupervised setting
gene_data_un = function(whole_shape = c(20,20,20), core_shape = c(3,3,3),dis,
                     gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1){ 
  #dis can be "gaussian" and "uniform"
  
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  ####-------- generate data
  set.seed(24)  # 24 # 37  #  347
  A = randortho(d1)[,1:r1]     ## factor matrix
  B = randortho(d2)[,1:r2]     ## factor matrix
  C = randortho(d3)[,1:r3]     ## factor matrix
  
  ### G: core tensor
  if(dis == "gaussian"){
    G = as.tensor(array(data = rnorm(r1*r2*r3,mean = gs_mean,sd = gs_sd),dim = core_shape))
  }
  else if(dis == "uniform"){
    G = as.tensor(array(data = runif(r1*r2*r3,min = unf_a,max = unf_b),dim = core_shape))
  }
  
  ### U: ground truth
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  
  ### ts:binary tensor
  ts = list()
  for (i in 1:5) {
    binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
    ts[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data
  }
  
  return(list(U = U,ts = ts))
}

#####---- This is the function used for generating data through different distribution
#         of core tensor in supervised setting
#         dup: number of dupplicate generated from one ground truth
gene_data = function(whole_shape = c(20,20,20), core_shape = c(3,3,3),p1,p2,dis,
                     gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1, dup){
  #dis can be "gaussian" and "uniform"
  
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  ####-------- generate data
  set.seed(24)  # 24 # 37  #  347
  
  X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 10),d1,p1)
  X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd = 10),d2,p2)
  W1 = randortho(p1,type = c("orthonormal"))[,1:r1]
  W2 = randortho(p2,type = c("orthonormal"))[,1:r2]
  
  C = randortho(d3,type = c("orthonormal"))[,1:r3]
  
  ### G: core tensor
  if(dis == "gaussian"){
    G = as.tensor(array(data = rnorm(r1*r2*r3,mean = gs_mean,sd = gs_sd),dim = core_shape))
  }
  else if(dis == "uniform"){
    G = as.tensor(array(data = runif(r1*r2*r3,min = unf_a,max = unf_b),dim = core_shape))
  }
  
  ### U: ground truth
  U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data
  
  ### tsr:binary tensor
  tsr = list()
  for (i in 1:dup) {
    binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
    tsr[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data
  }
  
  return(list(X_covar1 = X_covar1, X_covar2 = X_covar2, U = U,tsr = tsr))
}



conv_rate = function(d,r, p1, p2, dis,gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1, 
                     dup, Nsim, linear = TRUE, cons = 'vanilla' ,lambda = 1,
                     alpha = 1, solver = NULL){
  #cons can be "non","vanilla","penalty"
  rate = rep(0,length(d))
  RMSE = rep(0,length(d))
  for (i in 1:length(d)) {
    data = gene_data(rep(d[i],3), rep(r[i],3), p1[i], p2[i], dis, gs_mean, gs_sd, unf_a, unf_b, dup)
    X_covar1 = data$X_covar1
    X_covar2 = data$X_covar2
    U = data$U
    tsr = data$tsr
    RMSEi = rep(0,dup)
    for (j in 1:dup) {
      upp = update_binary(tsr = tsr[[j]], X_covar1 = X_covar1, X_covar2 = X_covar2, 
                          core_shape =  rep(r[i],3), Nsim, linear, cons, lambda = 1, 
                          alpha = 10*max(abs(U)), solver = NULL)
      
      U_est = ttl(upp$G,list(X_covar1%*%upp$W1,X_covar2%*%upp$W2,upp$C),ms = c(1,2,3))@data
      RMSEi[j] = sqrt(sum((U_est - U)^2)/(d[i]^3))
      print(paste(j,"-th observation ---- when dimension is ",d[i],"-- rank is ",r[i]," ---------"))
    }
    RMSE[i] = mean(RMSEi)
    rate[i] = r[i]^2/d[i]^2
  }
  return(list(RMSE = RMSE, rate = rate))
}







