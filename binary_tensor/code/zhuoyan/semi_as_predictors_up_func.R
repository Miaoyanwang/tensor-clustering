rm(list = ls())
library(rTensor)
library(pracma)
library(gtools)
library(rgl)
library("RColorBrewer")

#######################
# author: prof. Wang
glm_modify=function(y,x,start){
  
  ## initial coefficent
  ini_loglik=sum(log(inv.logit((2*y-1)*(x%*%start))))
  
  ## Option 1: glm fittig with default initilization
  fit1 = glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F))
  
  ## Option 2: glm with user specified initilization
  fit2= glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F),start=start)
  
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




#################  update

update_binary = function(tsr, Y_covar, core_shape, Nsim){
  
  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3] 
  p = dim(Y_covar)[2]
  
  
  ## get initialization
  
  C_ts_1 = glm_mat(Y_1,start = matrix(0,p,d2*d3),Y_covar)[[1]]
  C_ts = fold(C_ts_1, row_idx = 1, col_idx = c(2,3),modes = c(p,d2,d3))@data

  C_ts = as.tensor(C_ts)
  tckr = tucker(C_ts, ranks = core_shape)
  W = tckr$U[[1]] ; B = tckr$U[[2]] ; C = tckr$U[[3]]
  G = tckr$Z
  
  lglk = rep(0,4*Nsim)
  
  for(n in 1:Nsim){
    
    ###### update W
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    
    N_long = matrix(0,nrow = d1*d2*d3, ncol = p*r1)
    m=1
    for(j in 1:(d2*d3)){
      for(i in 1:d1){
        N_long[m,] =  kronecker_list(list(G_BC1[,j],Y_covar[i,]))
        m = m + 1
      }
    }
    
    
    mod_re = glm_modify(as.vector(Y_1), N_long, as.vector(W))
    coe = mod_re[[1]]
    W = matrix(coe, nrow = p, ncol = r1)
    lglk[4*n-3] = mod_re[[2]]
    print("W Done------------------")
    A = Y%*%W
    

    
    ##### update B
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
    
    #Y_2_and_Start = cbind(Y_2,B)
    #re = glm_mat(t(Y_2_and_Start),t(G_AC2))
    re = glm_mat(t(Y_2),start = t(B),t(G_AC2))
    
    B = t(re[[1]])
    lglk[4*n - 2] = re[[2]]
    ## orthogonal B*
    U = ttl(G,list(A,B,C),ms = c(1,2,3))
    tuk = tucker(U, ranks = core_shape)
    G = tuk$Z
    A = tuk$U[[1]]
    B = tuk$U[[2]]
    C = tuk$U[[3]]
    print("B Done------------------")
    
    ###### update C
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
    
    #Y_3_and_Start = cbind(Y_3,C)
    #re = glm_mat(t(Y_3_and_Start),t(G_AB3))
    re = glm_mat(t(Y_3),start = t(C),t(G_AB3))
    
    C = t(re[[1]])
    lglk[4*n - 1] = re[[2]]
    ## orthogonal C*
    U = ttl(G,list(A,B,C),ms = c(1,2,3))
    tuk = tucker(U, ranks = core_shape)
    G = tuk$Z
    A = tuk$U[[1]]
    B = tuk$U[[2]]
    C = tuk$U[[3]]
    print("C Done------------------")
    
    M_long = matrix(0,nrow = d1*d2*d3, ncol = r1*r2*r3)
    m=1
    ## update G
    for(k in 1:d3){
      for(j in 1:d2){
        for(i in 1:d1){
          M_long[m,] = kronecker_list(list(C[k,],B[j,],A[i,]))
          m = m + 1
        }
      }
    }
    
    
    
    mod_re = glm_modify(as.vector(tsr@data), M_long, as.vector(G@data))
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    lglk[4*n] = mod_re[[2]]
    print("G Done------------------")
    print(n)
    if(lglk[4*n] - lglk[4*(n-1)+1] <= 0.005) break
  }
  return(list(W = W,B = B,C = C,G = G,lglk = lglk))
}

up = update_binary(tensor,Y,c(4,4,5),5)
plot(up$lglk)

#### simulate

MSE = matrix(0,10,2)
colnames(MSE) = c('C','U')
for(n in 1:10){
  Y = matrix(rbinom(100,size = 1, prob = 0.5),20,5)
  W = randortho(5)[,1:2]
  B = randortho(20)[,1:2]
  C = randortho(20)[,1:2]
  G = as.tensor(array(data = rnorm(2*2*2),dim = c(2,2,2)))
  
  C_ts = ttl(G,list(W,B,C),ms = c(1,2,3))
  
  #U = ttl(G,list(Y%*%W,B,C),ms = c(1,2,3))@data
  U = ttm(C_ts,Y,m = 1)@data
  
  tsr = rbinom(20*20*20,1,prob = as.vector( 1/(1 + exp(-U)) ) )
  tsr = as.tensor(array(Y,dim = c(20,20,20)))@data
  
  up = update_binary(tsr,Y,c(2,2,2),10)
  #plot(up$lglk,ylim = c(min(up$lglk), max(up$lglk)))
  
  C_ts_hat = ttl(up$G,list(up$W,up$B,up$C),ms = c(1,2,3))
  MSE[n,1] = sum((C_ts@data - C_ts_hat@data)^2)
  
  U_hat = ttm(C_ts_hat,Y,m = 1)@data
  MSE[n,2] = sum((U - U_hat)^2)
  
  
}
MSE
plot(MSE[1:8,1],type = 'b')
plot(MSE[1:8,2],type = 'b')




























