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
  
  N_long = matrix(0,nrow = m*n, ncol = q1*q2)
  k=1
  for(j in 1:n){
    for(i in 1:m){
      N_long[k,] =  kronecker_list(list(X2[,j],X1[i,]))
      k = k + 1
    }
  }
  
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

update_binary = function(tsr, X_covar1, X_covar2, core_shape, Nsim){
  
  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3] 
  p_1 = dim(X_covar1)[2] ; p_2 = dim(X_covar2)[2]
  
  
  ## get initialization
  C_ts = glm_two_mat(tsr@data, X_covar1, t(X_covar2), ini = TRUE,linear=TRUE) ## add the linear model option for initilization


  # C_ts_1 = glm_mat(Y_1,start = matrix(0,p,d2*d3),X_covar)[[1]]
  # C_ts = fold(C_ts_1, row_idx = 1, col_idx = c(2,3),modes = c(p,d2,d3))@data
  
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
    
    
    # N_long = matrix(0,nrow = d1*d2*d3, ncol = p*r1)
    # m=1
    # for(j in 1:(d2*d3)){
    #   for(i in 1:d1){
    #     N_long[m,] =  kronecker_list(list(G_BC1[,j],X_covar[i,]))
    #     m = m + 1
    #   }
    # }
    # 
    # 
    # mod_re = glm_modify(as.vector(Y_1), N_long, as.vector(W))
    # coe = mod_re[[1]]
    # W = matrix(coe, nrow = p, ncol = r1)
    # lglk[4*n-3] = mod_re[[2]]
    
    
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
    
    # re = glm_mat(t(Y_2),start = t(B),t(G_AC2))
    # 
    # B = t(re[[1]])
    # lglk[4*n - 2] = re[[2]]
    
    ## orthogonal W2*
    C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))
    tuk = tucker(C_ts, ranks = core_shape)
    G = tuk$Z
    W1 = tuk$U[[1]]
    W2 = tuk$U[[2]]
    C = tuk$U[[3]]
    print("W2 Done------------------")
    
    
    # U = ttl(G,list(X_covar%*%W,B,C),ms = c(1,2,3))@data
    # U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    # 
    # L_W_hand[2*n]=sum(log(inv.logit((2*tsr@data-1)*U)))
    # 
    
    ##### calculate A
    A = X_covar1%*%W1
    
    ##### calculate B
    B = X_covar2%*%W2
    
    ###### update C
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
    
    #Y_3_and_Start = cbind(Y_3,C)
    #re = glm_mat(t(Y_3_and_Start),t(G_AB3))
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
    
    print(paste(n,"-th  iteration ---- when dimension is ",d1,"-- rank is ",r1," -----------------"))
    if(lglk[4*n] - lglk[4*n-1] <= 0.00005) break
    
    L_W_logLik[2*n + 1] = mod_re[[2]]  ### for checking W
    
  }
  return(list(W1 = W1,W2 = W2,C = C,G = G,lglk = lglk,L_W_logLik = L_W_logLik, L_W_hand = L_W_hand))
}


############  check on one simulation
set.seed(37)
X_covar1 = matrix(rbinom(100,size = 1, prob = 0.5),20,5)
X_covar2 = matrix(rbinom(100,size = 1, prob = 0.5),20,5)
W1 = randortho(5,type = c("orthonormal"))[,1:2]
W2 = randortho(5,type = c("orthonormal"))[,1:2]

C = randortho(20,type = c("orthonormal"))[,1:2]
G = as.tensor(array(data = rnorm(2*2*2,mean = 0,sd = 10),dim = c(2,2,2)))

C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))

U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data
#U = ttm(C_ts,X_covar,m = 1)@data

tsr = rbinom(20*20*20,1,prob = as.vector( 1/(1 + exp(-U)) ) )
tsr = as.tensor(array(tsr,dim = c(20,20,20)))@data

up = update_binary(tsr,X_covar1, X_covar2, c(2,2,2),10)

plot(up$lglk)
plot(up$lglk[up$lglk!=0])



####  codes below are useless, can be ignored

plot(up$L_W_logLik[2:20])
plot(up$L_W_hand[2:20])

C_ts_hat = ttl(up$G,list(up$W,up$B,up$C),ms = c(1,2,3))
sum((C_ts@data - C_ts_hat@data)^2)

#U_est = ttm(C_ts_hat,X_covar,m = 1)@data
U_est = ttl(up$G,list(X_covar%*%up$W,up$B,up$C),ms = c(1,2,3))@data
sum((U - U_est)^2)

hist(1/(1 + exp(-U)))

plot(as.vector(U),as.vector(U_est))
plot(as.vector(inv.logit(U)),as.vector(inv.logit(U_est)))
plot(as.vector(U[,,8]),as.vector(U_est[,,8]))
plot(as.vector(inv.logit(U[,,8])),as.vector(inv.logit(U_est[,,8])))


for(i in 1:20){
  png(paste('U_mode3_slice_10_siomoid',i,'.png'))
  plot(as.vector(inv.logit(U[,,i])),as.vector(inv.logit(U_est[,,i])))
  dev.off()
}


for(i in 1:20){
  png(paste('U_mode3_slice_10',i,'.png'))
  plot(as.vector(U[,,i]),as.vector(U_est[,,i]))
  dev.off()
}





#######################################################################
######################################################################
###----------- compare with unsupervised
set.seed(37)
X_covar = diag(20)
W = randortho(20)[,1:2]
B = randortho(20)[,1:2]
C = randortho(20)[,1:2]
G = as.tensor(array(data = rnorm(2*2*2,sd = 20),dim = c(2,2,2)))

C_ts = ttl(G,list(W,B,C),ms = c(1,2,3))

#U = ttl(G,list(Y%*%W,B,C),ms = c(1,2,3))@data
U = ttm(C_ts,Y,m = 1)@data

tsr = rbinom(20*20*20,1,prob = as.vector( 1/(1 + exp(-U)) ) )
tsr = as.tensor(array(tsr,dim = c(20,20,20)))@data

up2 = update_binary(tsr,X_covar,c(2,2,2),10)

plot(up2$lglk)


U_hat2 =  ttl(up2$G,list(X_covar%*%up2$W,up2$B,up2$C),ms = c(1,2,3))@data
sum((U - U_hat2)^2)


plot(as.vector(U),as.vector(U_hat2))
plot(as.vector(inv.logit(U)),as.vector(inv.logit(U_hat2)))

### unsupervised code
up2_un = update_binary(tsr,c(2,2,2),10)
plot(up2_un$lglk)
up2_un$lglk - up2$lglk


U_hat2_un =  ttl(up2$G,list(up2_un$A,up2_un$B,up2_un$C),ms = c(1,2,3))@data
sum((U - U_hat2_un)^2)
sum((U_hat2 - U_hat2_un)^2)



plot(as.vector(U_hat2),as.vector(U_hat2_un))
plot(as.vector(inv.logit(U_hat2)),as.vector(inv.logit(U_hat2_un)))


plot(as.vector(up2$W),as.vector(up2_un$A))
plot(as.vector(up2$B),as.vector(up2_un$B))
#plot(as.vector(inv.logit(up2$B)),as.vector(inv.logit(up2_un$B)))
plot(as.vector(up2$C),as.vector(up2_un$C))
plot(as.vector(up2$G@data),as.vector(up2_un$G@data))


#### simulate

MSE = matrix(0,20,2)
colnames(MSE) = c('C','U')
for(n in 1:8){
  Y = matrix(rbinom(100,size = 1, prob = 0.5),20,5)
  W = randortho(5)[,1:2]
  B = randortho(20)[,1:2]
  C = randortho(20)[,1:2]
  G = as.tensor(array(data = rnorm(2*2*2),dim = c(2,2,2)))
  
  C_ts = ttl(G,list(W,B,C),ms = c(1,2,3))
  
  #U = ttl(G,list(Y%*%W,B,C),ms = c(1,2,3))@data
  U = ttm(C_ts,Y,m = 1)@data
  
  tsr = rbinom(20*20*20,1,prob = as.vector( 1/(1 + exp(-U)) ) )
  tsr = as.tensor(array(tsr,dim = c(20,20,20)))@data
  
  up = update_binary(tsr,Y,c(2,2,2),10)
  #plot(up$lglk,ylim = c(min(up$lglk), max(up$lglk)))
  
  C_ts_hat = ttl(up$G,list(up$W,up$B,up$C),ms = c(1,2,3))
  MSE[n,1] = sum((C_ts@data - C_ts_hat@data)^2)
  
  U_hat = ttm(C_ts_hat,Y,m = 1)@data
  MSE[n,2] = sum((U - U_hat)^2)
  
  
}
MSE
plot(MSE[,1],type = 'b')
plot(MSE[,2],type = 'b')








