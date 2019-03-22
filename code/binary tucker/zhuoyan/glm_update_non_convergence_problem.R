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
  fit1 = glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = TRUE))

  ## Option 2: glm with user specified initilization
  fit2= glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = TRUE),start=start)

  ## report the result whichever gives the highest likelihood
  if(max(logLik(fit1),logLik(fit2))<ini_loglik) return (list(start, ini_loglik))
  else if(logLik(fit1)>logLik(fit2)) return(list(coef(fit1), logLik(fit1)))
  else return(list(coef(fit2), logLik(fit2)))
  
}

#######################


glm_f = function(y_and_start,X)  ## rbind y and initialized value
{
  n = dim(X)[1]
  R = dim(X)[2]
  y = y_and_start[1:n]
  start = y_and_start[(n+1):(n+R)]
  
  # ml = glm(y ~ -1 + X, family = binomial(link = 'logit'),
  #         control = glm.control(maxit = 100, trace = TRUE),start = start)
  # coe = ml$coefficients
  # lglk = logLik(ml)
  re = glm_modify(y, X, start)
  coe = re[[1]]
  lglk = re[[2]]
  return(list(coe,lglk))
}

#    form U = X*Beta
glm_mat = function(Y_and_Start,X){
  R = dim(X)[2]   ## R in note
  p = dim(Y_and_Start)[2]
  la = sapply(as.data.frame(Y_and_Start),glm_f,X)
  re = t(matrix(unlist(la), nrow = p, byrow=T))
  beta = re[1:R,]
  lglk = sum(re[R + 1,])
  return(list(beta,lglk))
}
############################################################
########################################################################

# glm_f = function(y,X,start)  ## rbind y and initialized value
# {
#   ml = glm(y ~ -1 + X, family = binomial(link = 'logit'),
#                     control = glm.control(maxit = 100, trace = F),start = start)
#   coe = ml$coefficients
#   lglk = logLik(ml)
#   return(list(coe,lglk))
# }
# 
# #    form U = X*Beta
# glm_mat = function(Y,X,start){
#   
#   Y
# 
#   la = mapply(glm_f,Y,X,start)
#   re = t(matrix(unlist(la), nrow = p, byrow=T))
#   beta = re[1:R,]
#   lglk = sum(re[R + 1,])
#   return(list(beta,lglk))
# }

#####  some test
set.seed(37)
X = matrix(rnorm(30,5,9),10,3)
Y = matrix(rbinom(20,1,0.5),10,2)
start = matrix(rnorm(6,1,0.5),3,2)



Y_and_Start = rbind(Y,matrix(rnorm(6,1,0.5),3,2))

y_and_start = Y_and_Start[,1]
glm_f(Y_and_Start[,1],X)

y_and_start = Y_and_Start[,2]
glm_f(Y_and_Start[,2],X)


re = glm_mat(Y_and_Start,X)
re[1]
re[[2]]
####



update_binary = function(ts, core_shape, Nsim){
  ## get initialization
  ts1 = 10*(2*ts - 1)
  ts1 = as.tensor(ts1)
  ts = as.tensor(ts)
  tckr = tucker(ts1, ranks = core_shape)
  A = tckr$U[[1]] ; B = tckr$U[[2]] ; C = tckr$U[[3]]
  G = tckr$Z
  d1 = dim(ts)[1] ; d2 = dim(ts)[2] ; d3 = dim(ts)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3] 
  Y_1 = unfold(ts, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(ts, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(ts, row_idx = 3, col_idx = c(1,2))@data
  
  lglk = rep(0,4*Nsim)
  for(n in 1:Nsim){
    
    ###### update A
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    
    Y_1_and_Start = cbind(Y_1,A)
    re = glm_mat(t(Y_1_and_Start),t(G_BC1))
    
    A = t(re[[1]])
    lglk[4*n - 3] = re[[2]]
    ## orthogonal A*
    U = ttl(G,list(A,B,C),ms = c(1,2,3))
    tuk = tucker(U, ranks = core_shape)
    G = tuk$Z
    A = tuk$U[[1]]
    B = tuk$U[[2]]
    C = tuk$U[[3]]
    print("A Done------------------")
    
    ##### update B
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
    
    Y_2_and_Start = cbind(Y_2,B)
    re = glm_mat(t(Y_2_and_Start),t(G_AC2))

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
    
    Y_3_and_Start = cbind(Y_3,C)
    re = glm_mat(t(Y_3_and_Start),t(G_AB3))
    
    
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
    
    # ml = glm(as.vector(ts@data) ~ -1 + M_long, family = binomial(link = 'logit'),
    #           control = glm.control(maxit = 150,trace = TRUE),start = as.vector(G@data))
    # 
    # G = as.tensor(array(data = ml$coefficients,dim = core_shape))
    # lglk[4*n] = logLik(ml)
    
    mod_re = glm_modify(as.vector(ts@data), M_long, as.vector(G@data))
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    lglk[4*n] = mod_re[[2]]
    
    print(n)
    if(lglk[4*n] - lglk[4*(n-1)+1] <= 0.005) break
  }
  return(list(A = A,B = B,C = C,G = G,lglk = lglk))
}

### apply on real data

library(R.matlab)
data = readMat('dnations.mat')
relname = unlist(data$relnnames)
att = unlist(data$attnames)
con = unlist(data$countrynames)
data = data$R
data[is.na(data)] = 0
data[1,1,]


up = update_binary(data,c(5,5,6),30)

plot(up$lglk,type = 'b')
up$lglk

kmA = kmeans(up$A,centers = 5)
kmB = kmeans(up$B,centers = 5)
kmC = kmeans(up$C,centers = 6)

table(kmeans(up$A,centers = 4)$cluster)

mem_A = model.matrix(~ -1 + factor(kmA$cluster))
rownames(mem_A) = con
mem_A = rbind(mem_A[mem_A[,1] == 1,],mem_A[mem_A[,2] == 1,],
      mem_A[mem_A[,3] == 1,],mem_A[mem_A[,4] == 1,],mem_A[mem_A[,5] == 1,])

mem_B = model.matrix(~ -1 + factor(kmB$cluster))
rownames(mem_B) = con
mem_B = rbind(mem_B[mem_B[,1] == 1,],mem_B[mem_B[,2] == 1,],
              mem_B[mem_B[,3] == 1,],mem_B[mem_B[,4] == 1,],mem_B[mem_B[,5] == 1,])

mem_C = model.matrix(~ -1 + factor(kmC$cluster))
rownames(mem_C) = relname
mem_C = rbind(mem_C[mem_C[,1] == 1,],mem_C[mem_C[,2] == 1,],mem_C[mem_C[,3] == 1,],
              mem_C[mem_C[,4] == 1,],mem_C[mem_C[,5] == 1,],mem_C[mem_C[,6] == 1,])



U_est = ttl(up$G,list(mem_A,mem_B,mem_C),ms = c(1,2,3))@data

plot_tensor(U_est)

U_est1 = U_est[,,1]
rownames(U_est1) = rownames(mem_A)
colnames(U_est1) = rownames(mem_B)



library(stats)
heatmap(U_est1,Rowv = NA, Colv = NA)
heatmap(U_est[,,11],Rowv = NA, Colv = NA)

kmC$cluster

################  simulation
simu = function(Nsim,core_shape, whole_shape){
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  
  MSE = seq(Nsim)
  U_true = list()
  U_hat = list()
  for(i in 1:Nsim){
    ### generate data
    G = as.tensor(array(data = rnorm(r1*r2*r3),dim = core_shape))
    A = randortho(d1)[,1:r1]
    B = randortho(d2)[,1:r2]
    C = randortho(d3)[,1:r3]
    
    U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    
    Y = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
    Y = as.tensor(array(Y,dim = c(d1,d2,d3)))@data
    
    up = update_binary(Y,c(r1,r2,r3),5)
    U_est = ttl(up$G,list(up$A,up$B,up$C),ms = c(1,2,3))@data
    
    MSE[i] = sum((U_est - U)^2)
    U_true[[i]] = U
    U_hat[[i]] = U_est
  }
  
  return(list(MSE = MSE,U_true = U_true, U_hat = U_hat))
}
simulation = simu(10,c(3,4,5),c(20,30,40))
plot(simulation[[1]],type = 'b')

U_true = simulation$U_true
U_hat = simulation$U_hat[7]


simulation$MSE
