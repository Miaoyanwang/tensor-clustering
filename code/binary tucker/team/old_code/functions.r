library(rTensor)
library(pracma)
library(gtools)
library(rgl)
library("RColorBrewer")

#######################
######----------   functions for update
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
glm_two = function(Y, X1, X2, ini = FALSE, start, linear=FALSE){ ## Y_size = m * n
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
  
  mod_re = glm_modify(as.vector(Y), N_long, coe_start)
  coe = mod_re[[1]]
  coe = matrix(coe, nrow = q1, ncol = q2)
  lglk= mod_re[[2]]
  return(list(coe = coe, lglk = lglk))
}

#####---- This function is a parallel version of GLM on two modes, used for initialization
## add the linear model option for initilization
glm_two_mat = function(Y, X1, X2, ini = TRUE, start = NULL,linear){
  Yl = lapply(seq(dim(Y)[3]), function(x) Y[ , , x])
  re = lapply(Yl, glm_two, X1, X2, ini = ini,linear=linear)
  
  coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
  coe = array(unlist(coe), dim = c(dim(X1)[2],dim(X2)[1],dim(Y)[3]))  ## form coe
  return(coe)
}



#################  update

update_binary = function(tsr, X_covar1, X_covar2, core_shape, Nsim, linear = TRUE){
  
  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3] 
  p_1 = dim(X_covar1)[2] ; p_2 = dim(X_covar2)[2]
  
  
  ## get initialization
  C_ts = glm_two_mat(tsr@data, X_covar1, t(X_covar2), ini = TRUE,linear=linear) ## add the linear model option for initilization
  
  C_ts = as.tensor(C_ts)
  tckr = tucker(C_ts, ranks = core_shape)
  W1 = tckr$U[[1]] ; W2 = tckr$U[[2]] ; C = tckr$U[[3]]
  G = tckr$Z
  A = X_covar1%*%W1
  B = X_covar2%*%W2
  
  lglk = rep(0,4*Nsim)
  L_W_hand = rep(0,2*Nsim+1)
  L_W_logLik = rep(0,2*Nsim+1)
  
  for(n in 1:Nsim){
    
    ###---------------------------------------------- check w1
    U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    L_W_hand[2*n-1]=sum(log(inv.logit((2*tsr@data-1)*U)))
    
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
    
    ###---------------------------------------------- check w1
    U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data
    L_W_hand[2*n]=sum(log(inv.logit((2*tsr@data-1)*U)))
    
    L_W_logLik[2*n] = mod_re[[2]]
    
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
    ## orthogonal C*
    C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))
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
    
    
    mod_re = glm_modify(as.vector(tsr@data), M_long, as.vector(G@data))
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    lglk[4*n] = mod_re[[2]]
    print("G Done------------------")
    print(n)
    
    print(paste(n,"-th  iteration ---- when dimension is ",d1,"-- rank is ",r1," -----------------"))
    if(lglk[4*n] - lglk[4*n-1] <= 0.00005) break
    
    L_W_logLik[2*n + 1] = mod_re[[2]]  ### for checking W
    
  }
  return(list(W1 = W1,W2 = W2,C = C,G = G,lglk = lglk,L_W_logLik = L_W_logLik, L_W_hand = L_W_hand))
}




###----  functions for simulation

gene_data = function(whole_shape = c(20,20,20), core_shape = c(3,3,3),dis,
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


conv_rate = function(d,r,Nsim = 50,cons,lambda = 1){
  #cons can be "non","CG","vanilla"
  rate = rep(0,length(d))
  RMSE = rep(0,length(d))
  for (i in 1:length(d)) {
    data = gene_data(rep(d[i],3), rep(r[i],3))
    U = data$U
    ts = data$ts
    RMSEi = rep(0,5)
    for (j in 1:5) {
      if(cons == "non"){
        upp = update_binary_non(ts[[j]],rep(r[i],3),Nsim)
      }
      else if(cons == "CG"){
        alpha = 10*max(abs(U))
        upp = update_binary_cons(ts[[j]],rep(r[i],3),Nsim, lambda = lambda, alpha = alpha)
      }
      else if(cons == "vanilla"){
        alpha = 10*max(abs(U))
        upp = update_binary_vanilla(ts[[j]],rep(r[i],3),Nsim, alpha = alpha)
      }
      U_est = ttl(upp$G,list(upp$A,upp$B,upp$C),ms = c(1,2,3))@data
      RMSEi[j] = sqrt(sum((U_est - U)^2)/(d[i]^3))
      print(paste(j,"-th observation ---- when dimension is ",d[i],"-- rank is ",r[i]," ---------"))
    }
    RMSE[i] = mean(RMSEi)
    rate[i] = r[i]^2/d[i]^2
  }
  return(list(RMSE = RMSE, rate = rate))
}


###
###----  This function shows how we select rank in non-constrain setting
# tsr, X_covar1, X_covar2, core_shape, Nsim, linear = TRUE
sele_rank = function(ts, X_covar1, X_covar2, rank, linear){
  BIC = rep(0,length(rank))
  for(i in 1:length(rank)){
    upp = update_binary(ts, X_covar1, X_covar2, rep(rank[i],3),100, linear)
    log_Lik = max(upp$lglk)
    BIC[i] = -2*log_Lik + (rank[i]^3 + sum((whole_shape-1)*rank[i]))*log(prod(whole_shape))
  }
  return(rank(which(BIC = min(BIC))))
}

#### This function is to select the lambda in the CG constrain version
select_lambda = function(ts,trueU,lambda){
  #tsr is an array
  #trueU is an array
  #have selected the rank
  
  d1 = dim(ts)[1]; d2 = dim(ts)[2]; d3 = dim(ts)[3]
  r1 = dim(trueU)[1] ; r2 = dim(trueU)[2] ; r3 = dim(trueU)[3]
  
  len = length(lambda)
  rmse = c()
  
  for (i in 1:len) {
    upp2 = update_binary_cons(ts,c(r1,r2,r3),Nsim = 50, lambda = lambda[i], alpha = 10* max(abs(trueU)))
    U_est2 = ttl(upp2$G,list(upp2$A,upp2$B,upp2$C),ms = c(1,2,3))@data
    
    rmse[i] =  sum((U_est2 - trueU)^2)/(d1*d2*d3)
  }
  
  best_lambda = lambda[order(rmse)[1]]
  return(best_lambda)
}








