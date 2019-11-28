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
HOPLS = HOPLS(X_covar1,tsr,10,c(3,3),0.01)


pre_err = function(HOPLS,tsr,Kn){
  pre = 0
  for(r in 1:R){
    inte = ttl(HOPLS$D[[r]], HOPLS$Q[[r]], m = seq(1,length(Kn)+1))
    pre = pre + inte
  }
  err = sum((pre@data - tsr)^2)
  return(err)
}
pre_err(HOPLS,tsr,c(2,2))













