rm(list = ls())
###  previous tensor regression algorithms

##----------------         HOLRR
library(rTensor)
## This part follows algorithm 1 in Guillaume et al. 2016
## HOLRR: for N-order tensor, N \geq 3
HOLRR = function(X,tsr,core_shape, gamma = 0){
  # core_shape = (d0, d1, ... dp)
  # X_shape =  N , d0
  # tsr_shape = N, d1, d2, ...,dp
  p = length(core_shape) - 1
  tsr = as.tensor(tsr)
  Y_1 = k_unfold(tsr, m = 1)@data
  I = diag(dim(X)[2])
  A = solve(t(X)%*%X + gamma*I)%*%t(X)%*%Y_1%*%t(Y_1)%*%X
  if(core_shape[1] == dim(X)[2]){
    U0 = diag(core_shape[1])
  }else{
    U0 = eigen(A)$vectors[,1:core_shape[1]]
  }
  U = list(); U[[1]] = U0
  for(i in seq(2,length(core_shape))){
    if(core_shape[i] == dim(tsr)[i]){
      U[[i]] = diag(core_shape[i])
    }else{
      Y_i = k_unfold(tsr, m = i)@data
      U[[i]] = eigen(Y_i%*%t(Y_i))$vectors[,1:core_shape[i]]
    }
  }
  M = solve(t(U0)%*%(t(X)%*%X + gamma*I)%*%U0)%*%t(U0)%*%t(X)
  
  U_M = lapply(seq(length(core_shape)), function(x) t(U[[x]]))
  U_M[[1]] = M
  
  G = ttl(tsr, U_M, ms = seq(p+1))
  C_ts = ttl(G, U, ms = seq(p+1))
  return(list(G = G@data, U = U, C_ts = C_ts@data))
}


##------------------------   TPG + ITP
library(rTensor)
## This part follows algorithm 2 in Rose Yu et al. 2016
## ITP: for N-order tensor, N = 3
ITP = function(W, R, U){
  # R: scalar constraining tucker rank
  # W: Coefficient tensor: d1, d2, T
  # U: list containing factor matrices
  W = as.tensor(W)
  U1 = U[[1]]; U2 = U[[2]]; U3 = U[[3]]
  i = 1
  while(i <= R){
    u10 = U1[,i]; u20 = U2[,i]; u30 = U3[,i]
    for(j in 1:100){  ##  power iteration
      u11 = ttl(W,list(t(as.matrix(u20)), t(as.matrix(u30))), ms = c(2,3))@data
      u11 = as.vector(u11); u11 = u11/sum(u11^2)
      u21 = ttl(W,list(t(as.matrix(u11)), t(as.matrix(u30))), ms = c(1,3))@data
      u21 = as.vector(u21); u21 = u21/sum(u21^2)
      u31 = ttl(W,list(t(as.matrix(u11)), t(as.matrix(u21))), ms = c(1,2))@data
      u31 = as.vector(u31); u31 = u31/sum(u31^2)
      u10 = u11; u20 = u21; u30 = u31
      if((sum((u11 - u10)^2) + sum((u21 - u20)^2) + sum((u31 - u30)^2)) <= 0.05) break
    }
    U1[,i] = u11;  U2[,i] = u21;  U3[,i] = u31 
    i = i + 1
  }
  U = list(U1,U2,U3)
  UUT = lapply(1:3, function(x) U[[x]]%*%t(U[[x]]))
  W = ttl(W, UUT, m = c(1,2,3))@data
  return(list(W = W, U = U))
}

sketch = function(x){
  i = sample(1:length(x),1)
  x[i] = sample(c(1,-1),1,prob = c(0.5,0.5))
  return(x)
}

## This part follows algorithm 1 in Rose Yu et al. 2016
## TPG: for N-order tensor, N = 3
TPG = function(X, tsr, R, eta, iter, tol, seed = 1){
  # tsr: d3, d2, T
  # X: d3, d1
  # W: Coefficient tensor: d1, d2, T
  # eta: step size
  # iter: iteration times
  # tol: tolerance
  set.seed(seed)
  d3 = dim(X)[1]; Tt = dim(tsr)[3]; d2 = dim(tsr)[2]; d1 = dim(X)[2]
  l = sample(seq(d3-1))  # sketch to lower dim
  S = matrix(0,l, d3)
  S = as.matrix(sapply(as.data.frame(S), sketch))
  tsr = as.tensor(tsr)
  tsr = ttm(tsr,S,m = 1)
  X = S %*% X
  W = array(0, c(d1,d2,Tt));  W = as.tensor(W)
  gr = ttm(tsr - ttm(W,X,m = 1),t(X),m=1)
  W = W - eta*gr
  U = tucker(W, ranks = rep(R,3))$U
  W = W@data
  loss = 0
  for(i in 1:iter) {
    loss0 = loss
    gr = -ttm(tsr - ttm(as.tensor(W),X,m = 1),t(X),m=1)@data
    W = W - eta*gr
    ITPres = ITP(W, R, U)
    W = ITPres$W
    U = ITPres$U
    loss = sum((tsr@data - ttm(as.tensor(W),X,m = 1)@data)^2)
    cat("step", i, ":", loss, "\n")
    #print(W)
    if(loss <= tol | sum((loss - loss0)^2)<= 0.01 ){
      break
    }
  }
  return(W)
}
###----------------------------   HOPLS
## This part follows algorithm 1+2 --> 3 in Q. Zhao et al. 2012
## HOPLS: for N-order tensor, N \geq 3
library(rTensor)
HOPLS = function(X,tsr, R, Kn, tol){
  # R: number of latent vectors
  # Kn: vector, contains number of loadings in tsr (see fig 2 in Q. Zhao et al. 2012)
  # tol: tolerance, epsilon in algorithm
  E = list(); Fr = list() ## denote residuals of X, tsr
  t = list() ## latent vectors
  Q = list() ## loadings
  D = list() ## tensor of tsr (see fig 2 in Q. Zhao et al. 2012)
  E[[1]] = X; Fr[[1]] = as.tensor(tsr)
  for(r in 1:R){
    if(sum(E[[r]]^2) > tol & sum(Fr[[r]]@data^2) > tol){
      Cr = ttm(Fr[[r]],t(E[[r]]), m = 1)
      tckr = tucker(Cr, ranks = c(1,Kn))
      t_inte = E[[r]]%*%tckr$U[[1]]  ## pr = tckr$U[[1]]  vector
      t[[r]] = t_inte/sqrt(sum(t_inte^2))
      Q[[r]] = tckr$U ; Q[[r]][[1]] = t[[r]]
      Qrt = lapply(1:(1+length(Kn)), function(x) t(Q[[r]][[x]]))  ## get t(Q[[r]])
      #D_inte = ttm(Fr[[r]],t(t[[r]]), m = 1)
      #D[[r]] = ttl(D_inte,Qrt,ms = seq(2,length(Kn)+1))
      D[[r]] = ttl(Fr[[r]],Qrt,ms = seq(1,length(Kn)+1))
      E[[r+1]] = E[[r]] - t[[r]]%*%t(tckr$U[[1]])
      #inte = ttm(D[[r]], t[[r]], m = 1)
      #Fr[[r+1]] = Fr[[r]] - ttl(inte, Q[[r]], m = seq(2,length(Kn)+1))
      Fr[[r+1]] = Fr[[r]] - ttl(D[[r]], Q[[r]], m = seq(1,length(Kn)+1))
    }else break
  }
  return(list(Q = Q, D = D, t = t))
}
#HOPLS = HOPLS(X_covar1,tsr,10,c(3,3),0.01)


pre_err = function(HOPLS,tsr,Kn){
  pre = 0
  for(r in 1:R){
    inte = ttl(HOPLS$D[[r]], HOPLS$Q[[r]], m = seq(1,length(Kn)+1))
    pre = pre + inte
  }
  err = sum((pre@data - tsr)^2)
  return(err)
}
#pre_err(HOPLS,tsr,c(2,2))



