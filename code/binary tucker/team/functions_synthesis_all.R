library(rTensor)
library(pracma)
library(gtools)
#library(rgl)
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

loglike=function(data,linearpre){
  p=plogis(linearpre)
  p[p==1] = p[p==1] - 0.0001
  p[p==0] = p[p==0] + 0.0001
  L=data*log(p)+(1-data)*log(1-p)## log-likelihood
  return(sum(L)) 
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
glm_modify=function(y,x,start = NULL){
  
  if(is.null(start)){
    fit1 = suppressWarnings(glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F)))
    return(list(coef(fit1), logLik(fit1)))
  }
  
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

#    form U = X*Beta
glm_mat = function(Y,start,X){
  R = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  ma = mapply(glm_modify, y =  as.data.frame(Y),start = as.data.frame(start),MoreArgs = list(X))
  re = t(matrix(unlist(ma), nrow = p, byrow=T))
  beta = re[1:R,]
  lglk = sum(re[R + 1,])
  return(list(t(beta),lglk))
}

###########---------  GLM on two modes
##---  This section follows 1.3 of [Algorithm: Semi-Supervised Binary Tensor Factorization]
## add the linear model option for initilization
##---  function version
glm_two = function(Y, X1, X2, start = NULL, linear=FALSE, lm = FALSE # used for supervised scale down
){ ## Y_size = m * n
  #X2 = t(X2)
  # logit(E(Y)) = X1 %*% coe %*% X2
  m = dim(Y)[1] ; n = dim(Y)[2]
  q1 = dim(X1)[2] ; q2 = dim(X2)[1]
  
  N_long = kronecker_list(list(t(X2),X1))
  
  
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
  
  mod_re = glm_modify(as.vector(Y), N_long, as.vector(start))
  coe = mod_re[[1]]
  coe = matrix(coe, nrow = q1, ncol = q2)
  lglk= mod_re[[2]]
  return(list(coe = coe, lglk = lglk))
}

#####---- This function is a parallel version of GLM on two modes, used for initialization
## add the linear model option for initilization
glm_two_mat = function(Y, X1, X2, start = NULL,linear, lm = FALSE){
 Yl = lapply(seq(dim(Y)[3]), function(x) Y[ , , x])
 re = lapply(Yl, glm_two, X1, X2,linear=linear, lm = lm)
 
 coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
 coe = array(unlist(coe), dim = c(dim(X1)[2],dim(X2)[1],dim(Y)[3]))  ## form coe
 return(coe)
}

glm_two_ini = function(Y, X1, X2, X3, start = NULL,linear, lm = FALSE, un_mode = 3){
  if(un_mode == 1){
    Yl = lapply(seq(dim(Y)[1]), function(x) Y[ x, , ])
    re = lapply(Yl, glm_two, X2, t(X3),linear=linear, lm = lm)
    coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
    coe = aperm(array(unlist(Yl),dim = c(dim(X2)[2],dim(X3)[2],dim(Y)[1])), perm = c(3,1,2) )
    
  } else if(un_mode == 2) {
    Yl = lapply(seq(dim(Y)[2]), function(x) Y[ , x, ])
    re = lapply(Yl, glm_two, X1, t(X3),linear=linear, lm = lm)
    coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
    coe = aperm(array(unlist(Yl),dim = c(dim(X1)[2],dim(X3)[2],dim(Y)[2])), perm = c(1,3,2) )
    
  } else if(un_mode == 3) {
    Yl = lapply(seq(dim(Y)[3]), function(x) Y[ , , x])
    re = lapply(Yl, glm_two, X1, t(X2),linear=linear, lm = lm)
    coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
    coe = array(unlist(coe), dim = c(dim(X1)[2],dim(X2)[2],dim(Y)[3]))  ## form coe
    
  }
  
  # re = lapply(Yl, glm_two, X1, X2,linear=linear, lm = lm)
  # coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
  # coe = array(unlist(coe), dim = c(dim(X1)[2],dim(X2)[1],dim(Y)[3]))  ## form coe
  return(coe)
}


#################  update
###--------   unsupervised
update_binary_all = function(tsr,X_covar1 = NULL, X_covar2 = NULL,X_covar3 = NULL, 
                             core_shape, Nsim=20, linear = TRUE, cons, lambda = 0.1,
                             alpha = 1, solver = NULL){
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
  
  
  ## get initialization
  if(un_m1 & un_m2 & un_m3){
    C_ts = 10*(2*tsr - 1)
  }  else if(un_m1) {
    C_ts = glm_two_ini(tsr@data,X1 = NULL,X2 = X_covar2,X3 = X_covar3, un_mode = 1,linear=linear) ## add the linear model option for initilization
    C_ts = as.tensor(C_ts)
  } else if(un_m2) {
    C_ts = glm_two_ini(tsr@data,X1 = X_covar1,X2 = NULL,X3 = X_covar3, un_mode = 2,linear=linear) ## add the linear model option for initilization
    C_ts = as.tensor(C_ts)
  } else if(un_m3) {
    C_ts = glm_two_ini(tsr@data,X1 = X_covar1,X2 = X_covar2, X3 = NULL, un_mode = 3,linear=linear) ## add the linear model option for initilization
    C_ts = as.tensor(C_ts)
  } else {
    M_long = kronecker_list(list(X_covar3,X_covar2,X_covar1)) ## form M_long
    
    if(linear==TRUE){
      mod_re = lm(10*(2*as.vector(tsr@data)-1)~-1+M_long)
    } else {
      mod_re = glm_modify(as.vector(tsr@data), M_long)
    }
    C_ts = as.tensor(array(data = mod_re[[1]],dim = c(p1,p2,p3)))
  }
  
  tckr = tucker(C_ts, ranks = core_shape)
  W1 = tckr$U[[1]] ; W2 = tckr$U[[2]] ; W3 = tckr$U[[3]]
  G = tckr$Z
  A = X_covar1%*%W1
  B = X_covar2%*%W2
  C = X_covar3%*%W3
  
  violate = c() ## which iteration violate constrain
  lglk0 = -Inf
  lglk = c()
  
  for(n in 1:Nsim){
    ## parameter from previous step
    
    W10 = W1 ; W20 = W2 ; W30 = W3 ; G0=G
    if(n>=2) lglk0=tail(lglk, n=1)
    
    ###### update W1
    
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    
    if(un_m1) {re = glm_mat(t(Y_1),start = t(A),t(G_BC1))
    } else {re = glm_two(Y = Y_1, X1 = X_covar1, X2 = G_BC1, start = W1)}
    
    if(dim(W1)[2]==1) W1 = as.matrix(re[[1]]) else W1 = re[[1]]
    
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W1*
    C_ts = ttl(G,list(W1,W2,W3),ms = c(1,2,3))
    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    W3 = tuk$U[[3]]
    print("W1 Done------------------")
    
    ##### calculate A ; B ; C
    A = X_covar1%*%W1 ; B = X_covar2%*%W2 ; C = X_covar3%*%W3
    
    
    ##### update W2
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
    
    if(un_m2) {re = glm_mat(t(Y_2),start = t(B),t(G_AC2))
    } else {re = glm_two(Y_2, X_covar2, G_AC2, start = W2)}
    
    if(dim(W2)[2]==1) W2 = as.matrix(re[[1]]) else W2 = re[[1]]
    
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W2*
    C_ts = ttl(G,list(W1,W2,W3),ms = c(1,2,3))
    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    W3 = tuk$U[[3]]
    print("W2 Done------------------")
    
    ##### calculate A ; B ; C
    A = X_covar1%*%W1 ; B = X_covar2%*%W2 ; C = X_covar3%*%W3
    
    
    ###### update W3
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
    
    if(un_m3) {re = glm_mat(t(Y_3),start = t(C),t(G_AB3))
    } else {re = glm_two(Y_3, X_covar3, G_AB3, start = W3)}

    if(dim(W3)[2]==1) W3 = as.matrix(re[[1]]) else W3 = re[[1]]
    
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W3*
    C_ts = ttl(G,list(W1,W2,W3),ms = c(1,2,3))
    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    W3 = tuk$U[[3]]
    print("W3 Done------------------")
    
    ##### calculate A ; B ; C
    A = X_covar1%*%W1 ; B = X_covar2%*%W2 ; C = X_covar3%*%W3
    
    #########-----------------------------------------------
    ###  then we apply constrain
    ######---- differnent version of contrains
    
    U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    ## scale down the core tensor. 
    if(cons == 'non' | max(abs(U))<= alpha | cons == 'vanilla') {U=U
    } else if (cons == 'penalty'){
      G=G/max(abs(U))*(alpha-0.01)
      U = U/max(abs(U))*(alpha-0.01)
      print("Violate constrain ------------------")
      violate = c(violate,n)
    }
    
    ##### update G
    M_long = kronecker_list(list(C,B,A)) ## form M_long
    
    if(cons == 'penalty'){
      mod_re = optim(par = as.vector(G@data),loss,loss_gr,y = as.vector(tsr@data),
                     X = M_long, lambda = lambda, alpha = alpha, method = solver)
      coe = mod_re$par
      G = as.tensor(array(data = coe,dim = core_shape))
      lglk = c(lglk,-mod_re$value) 
    } else {
      mod_re = glm_modify(as.vector(tsr@data), M_long, as.vector(G@data))
      coe = mod_re[[1]]
      G = as.tensor(array(data = coe,dim = core_shape))
      U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
      
      if ((cons== 'vanilla')&(max(abs(U))>=alpha)){
        G=G/max(abs(U))*alpha
        U=U/max(abs(U))*alpha
        print("Violate constrain ------------------")
        violate = c(violate,n)
      }
      lglk = c(lglk,loglike(tsr@data,U))
    }
    
    print("G Done------------------")
    
    print(paste(n,"-th  iteration ---- when dimension is ",d1,"-- rank is ",r1," -----------------"))
    
    if (tail(lglk,1)-lglk0 <= 0.0005 & tail(lglk,1)>= lglk0 ) break
    else if (tail(lglk,1)-lglk0 < 0) {
      W1 = W10 ; W2 = W20 ; W3 = W30; G=G0; lglk=lglk; break
    } 
    
  }
  return(list(W1 = W1,W2 = W2,W3 = W3,G = G,lglk = lglk, violate = violate))
}


update_binary = function(tsr, X_covar1 = NULL, X_covar2 = NULL, core_shape, Nsim, linear = TRUE,
                         cons = 'vanilla', lambda = 1, alpha = 1, solver = NULL){
  
  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  if(is.null(X_covar1)) X_covar1 = diag(d1)
  if(is.null(X_covar2)) X_covar2 = diag(d2)
  p_1 = dim(X_covar1)[2] ; p_2 = dim(X_covar2)[2]
  
  
  ## get initialization
  C_ts = glm_two_mat(tsr@data, X_covar1, t(X_covar2),linear=linear) ## add the linear model option for initilization
  
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
    
    C = re[[1]]
    lglk[4*n - 1] = re[[2]]
    
    
    
    ###  then we apply out constrain
    ######---- differnent version of contrains
    U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data
    
    if(cons == 'non'){C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))}
    else if(max(abs(U)) <= alpha){C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))}
    else if(cons == 'vanilla'){
      U = U/max(abs(U))*alpha
      C_ts = glm_two_mat(U, X_covar1, t(X_covar2),linear=FALSE, lm = TRUE) ## add the linear model option for initilization
      C_ts = as.tensor(C_ts)
      print("Violate constrain ------------------")
    }
    else{
      U = U/max(abs(U))*(alpha-0.01)
      C_ts = glm_two_mat(U, X_covar1, t(X_covar2),linear=FALSE, lm = TRUE) ## add the linear model option for initilization
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
    
    #U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data
    #print(max(abs(U))/alpha)
    
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

# sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, rank = c(3,5,6,8), Nsim,
#                      linear = FALSE, cons = 'non'){
#   BIC = rep(0,length(rank))
#   whole_shape = dim(tsr)
#   for(i in 1:length(rank)){
#     if((is.null(X_covar1)) & (is.null(X_covar2)) | (all.equal(X_covar1,diag(dim(tsr)[1])) & all.equal(X_covar2,diag(dim(tsr)[2])))){
#       upp = update_binary_un(tsr, rep(rank[i],3), Nsim, cons = cons)
#       BIC[i] = -2*log_Lik + (rank[i]^3 + sum((whole_shape-1)*rank[i]))*log(prod(whole_shape))
#     }
#     else {
#       upp = update_binary(tsr, X_covar1, X_covar2, rep(rank[i],3), Nsim, linear, cons = cons)
#       BIC[i] = -2*log_Lik + (rank[i]^3 + sum((c(dim(X_covar1)[2],dim(X_covar2)[2],whole_shape[3])-1)*rank[i]))*log(prod(whole_shape))
#     }
#     log_Lik = max(upp$lglk)
#     
#     
#   }
#   return(rank = which(BIC == min(BIC)))
# }


sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,
                     rank = c(3,5), Nsim,linear = FALSE, cons = 'non'){
  rank = expand.grid(rank,rank,rank)
  colnames(rank) = NULL
  rank = as.matrix(rank)
  whole_shape = dim(tsr)
  rank = lapply(1:dim(rank)[1], function(x) rank[x,]) ## turn rank to a list
  upp = lapply(rank, FUN= update_binary_all,tsr = tsr,X_covar1 = X_covar1,X_covar2 = X_covar2,
               X_covar3 = X_covar3, Nsim = Nsim, linear = linear, cons = cons)
  
  log_lik = unlist(lapply(seq(length(upp)), function(x) tail(upp[[x]]$lglk,1)))
  if(is.null(X_covar1)) X_covar1 = whole_shape[1]
  if(is.null(X_covar2)) X_covar2 = whole_shape[2]
  if(is.null(X_covar3)) X_covar3 = whole_shape[3]
  
  
  BIC = unlist(lapply(seq(length(rank)), 
                      function(x) (prod(rank[[x]]) + sum((c(dim(X_covar1)[2],dim(X_covar2)[2],
                                                            dim(X_covar3[2]))-1)*
                                                           rank[[x]])) * log(prod(whole_shape))))
  BIC = -2*log_lik + BIC
  # BIC = rep(0,length(rank))
  # for(i in 1:length(rank)){
  #     BIC[i] =  (prod(rank[[i]]) + sum((c(dim(X_covar1)[2],dim(X_covar2)[2],whole_shape[3])-1)*
  #                                  rank[[i]]))*log(prod(whole_shape))
  #   }
  return(rank = rank[[which(BIC == min(BIC))]])
}



#### This function is to select the lambda in the penalty constrain version
sele_lambda = function(seed, lambda, ...){
  #lambda = as.list(lambda)
  re = lapply(lambda, FUN = conv_rate, seed = seed, ...)
  re = lapply(seq(length(re)), function(x) re[[x]]$RMSE)
  return(re)
}
# select_lambda = function(tsr,trueU,lambda, X_covar1, X_covar2, Nsim, linear = TRUE){
#   #lambda is a sequence of lambda options 
#   #tsr is an array
#   #trueU is an array
#   #have selected the rank
#   
#   d1 = dim(tsr)[1]; d2 = dim(tsr)[2]; d3 = dim(tsr)[3]
#   r1 = dim(trueU)[1] ; r2 = dim(trueU)[2] ; r3 = dim(trueU)[3]
#   
#   len = length(lambda)
#   rmse = c()
#   
#   for (i in 1:len) {
#     upp2 = update_binary(tsr, X_covar1, X_covar2, core_shape = c(r1,r2,r3), Nsim, linear = linear,
#                          cons = "penalty", lambda = lambda[i], alpha = 10*max(abs(trueU)), solver = "CG")
#     U_est2 = ttl(upp2$G,list(upp2$A,upp2$B,upp2$C),ms = c(1,2,3))@data
#     
#     rmse[i] =  sum((U_est2 - trueU)^2)/(d1*d2*d3)
#   }
#   
#   best_lambda = lambda[order(rmse)[1]]
#   return(best_lambda)
# }




###----  functions for simulation
#####---- This is the function used for generating data through different distribution
#         of core tensor in unsupervised setting
gene_data_all = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),
                        p1 = NULL,p2 = NULL, p3 = NULL,dis, gs_mean = 0,gs_sd = 10,
                        unf_a = 0,unf_b = 1,dup, signal){ 
  #dis can be "gaussian" and "uniform"
  #print(p3)
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  ####-------- generate data
  set.seed(seed)  # 24 # 37  #  347

  if(is.null(p1)) {
    A = as.matrix(randortho(d1)[,1:r1])     ## factor matrix
    X_covar1 = NULL ; W1 = NULL
  } else {
    X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 1),d1,p1)
    W1 = randortho(p1,type = c("orthonormal"))[,1:r1]
    A = X_covar1%*%W1
  }
  if(is.null(p2)) {
    B = as.matrix(randortho(d2)[,1:r2])     ## factor matrix
    X_covar2 = NULL ; W2 = NULL
  } else {
    X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd = 1),d2,p2)
    W2 = randortho(p2,type = c("orthonormal"))[,1:r2]
    B = X_covar2%*%W2
  }
  if(is.null(p3)) {
    C = as.matrix(randortho(d3)[,1:r3])     ## factor matrix
    X_covar3 = NULL ; W3 = NULL
  } else {
    X_covar3 = matrix(rnorm(d3*p3,mean = 0, sd = 1),d3,p3)
    W3 = randortho(p3,type = c("orthonormal"))[,1:r3]
    C = X_covar3%*%W3
  }
  
  ### G: core tensor
  if(dis == "gaussian"){
    G = as.tensor(array(data = rnorm(r1*r2*r3,mean = gs_mean,sd = gs_sd),dim = core_shape))
  }
  else if(dis == "uniform"){
    G = as.tensor(array(data = runif(r1*r2*r3,min = unf_a,max = unf_b),dim = core_shape))
  }
  
  ### U: ground truth
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  G=G/max(abs(U))*signal
  U=U/max(abs(U))*signal
  
  
  ### tsr:binary tensor
  tsr = lapply(seq(dup), function(x) array(rbinom(d1*d2*d3,1,
                                                  prob = as.vector( 1/(1 + exp(-U)))) ,
                                                  dim = c(d1,d2,d3)))
  
  # tsr = list()
  # for (i in 1:dup) {
  #   binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
  #   tsr[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data
  # }
  return(list(U = U,tsr = tsr,X_covar1 = X_covar1, X_covar2 = X_covar2,X_covar3 = X_covar3,
              W1 = W1,W2 = W2,W3 = W3, A=A,B=B,C=C,G=G))
}


# data = gene_data_all(seed = 24,whole_shape = rep(20,3),core_shape =  rep(3,3), p1 = 5,p2 =5 ,
#                      'gaussian', gs_mean = 0, gs_sd = 10, dup = 2,signal = 10)
# 
# data = gene_data_all(seed = 24,whole_shape = rep(20,3),core_shape =  rep(3,3), p1 = 5,p2 =5 ,
#                      dis = 'gaussian', gs_mean = 0, gs_sd = 10, dup = 2,signal = 10)
# 
# data = gene_data_all(seed = 24,whole_shape = rep(20,3),core_shape =  rep(3,3), p1 = 5,p2 =5 ,
#                      p3 = NULL,
#                      'gaussian', gs_mean = 0, gs_sd = 10, dup = 2,signal = 10)
# 


#####---- This is the function used for generating data through different distribution
#         of core tensor in supervised setting
#         dup: number of dupplicate generated from one ground truth
gene_data = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),p1,p2,dis,
                     gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1, dup){
  #dis can be "gaussian" and "uniform"
  
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  ####-------- generate data
  set.seed(seed)  # 24 # 37  #  347
  
  X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 1),d1,p1)
  X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd = 1),d2,p2)
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
  
  ### true coefficient
  C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))@data
  ### U: ground truth
  U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data
  
  ### tsr:binary tensor
  tsr = lapply(seq(dup), function(x) array(rbinom(d1*d2*d3,1,
                                                  prob = as.vector( 1/(1 + exp(-U)))) , dim = c(d1,d2,d3)))
  
  # tsr = list()
  # for (i in 1:dup) {
  #  binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
  #  tsr[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data}
  # 
  return(list(X_covar1 = X_covar1, X_covar2 = X_covar2,C_ts = C_ts, U = U,tsr = tsr))
}




conv_rate = function(seed,d,r, p1 = NULL, p2 = NULL,p3 = NULL , dis,gs_mean = 0,
                     gs_sd = 10,unf_a = 0,unf_b = 1,dup,signal, Nsim, linear = TRUE,
                     cons = 'vanilla' ,lambda = 1,solver = NULL){
  #cons can be "non","vanilla","penalty"
  rate = rep(0,length(d))
  RMSE = rep(0,length(d))
  for (i in 1:length(d)) {
    data = gene_data_all(seed,rep(d[i],3), rep(r[i],3), p1[i], p2[i],p3 = p3, dis, gs_mean, gs_sd,
                         unf_a, unf_b, dup,signal)
    X_covar1 = data$X_covar1
    X_covar2 = data$X_covar2
    C_ts = ttl(data$G,list(data$W1,data$W2,data$C),ms = c(1,2,3))@data
    U = data$U
    tsr = data$tsr
    RMSEi = rep(0,dup)
    for (j in 1:dup) {
      upp = update_binary(tsr = tsr[[j]], X_covar1 = X_covar1, X_covar2 = X_covar2, 
                          core_shape =  rep(r[i],3), Nsim, linear, cons, lambda = lambda, 
                          alpha = 10*max(abs(U)), solver = NULL)
      
      C_ts_est = ttl(upp$G,list(upp$W1,upp$W2,upp$C),ms = c(1,2,3))@data
      #U_est = ttl(upp$G,list(X_covar1%*%upp$W1,X_covar2%*%upp$W2,upp$C),ms = c(1,2,3))@data
      RMSEi[j] = sqrt(sum((C_ts_est - C_ts)^2))
      print(paste(j,"-th observation ---- when dimension is ",d[i],"-- rank is ",r[i]," ---------"))
      print(paste("p1 is ",p1[i],"---------  p2 is ",p2[i],"--------------------"))
    }
    RMSE[i] = mean(RMSEi)
    rate[i] = sqrt(r[i]^2*(d[i] + p1[i] + p2[i])/d[i]^2)
  }
  return(list(RMSE = RMSE, rate = rate))
}

# con = conv_rate(seed = 24,d = seq(20,70,10),r = rep(3,7), p1 = rep(5,7), p2 = rep(5,7),p3 = NULL, dis = 'gaussian',gs_mean = 0,gs_sd = 10, 
#                 dup=2, signal = 10,Nsim =50, linear = FALSE, cons = 'penalty' ,lambda = 1,
#                 solver = 'CG')
# 
