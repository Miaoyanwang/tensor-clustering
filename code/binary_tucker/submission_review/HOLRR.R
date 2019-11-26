## This part follows algorithm 1 in Guillaume et al. 2016
## HOLRR
library(rTensor)
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

  