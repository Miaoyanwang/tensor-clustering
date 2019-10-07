####  Unsupervised checking
# This code follows the manual of [Algorithm:  Unsupervised Binary Tensor Factorization]

rm(list = ls())
library(rTensor)
library(pracma)
library(gtools)
library(rgl)
library("RColorBrewer")

########-------------------------------------------
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
# This part follows Sec 1.1 to implement matrix GLM
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
########################################################################

#################  update
# This part follows Sec 1.3 to update core tensor and factor matrix

update_binary_raw = function(ts, core_shape, Nsim){
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
    
    re = glm_mat(t(Y_1),start = t(A),t(G_BC1))
    
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
    
    ##### update G
    M_long = matrix(0,nrow = d1*d2*d3, ncol = r1*r2*r3)
    m=1
    for(k in 1:d3){
      for(j in 1:d2){
        for(i in 1:d1){
          M_long[m,] = kronecker_list(list(C[k,],B[j,],A[i,]))
          m = m + 1
        }
      }
    }
    
    
    
    mod_re = glm_modify(as.vector(ts@data), M_long, as.vector(G@data))
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    lglk[4*n] = mod_re[[2]]
    print("G Done------------------")
    print(paste(n,"-th  iteration -----------------------------"))
    if(lglk[4*n] - lglk[4*n-1] <= 0.00005) break
    
  }
  return(list(A = A,B = B,C = C,G = G,lglk = lglk))
}




########------------   simulatioin for 1 time
whole_shape = rep(20,3) ; core_shape = rep(3,3)
d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]

####-------- generate data
set.seed(24)  # 24 # 37  #  347
A = randortho(d1)[,1:r1]     ## factor matrix
B = randortho(d2)[,1:r2]     ## factor matrix
C = randortho(d3)[,1:r3]     ## factor matrix

### G: core tensor
G = as.tensor(array(data = rnorm(r1*r2*r3,mean = 0,sd = 10),dim = core_shape))
###
#G = as.tensor(array(data = runif(r1*r2*r3,0,1),dim = core_shape))

### U: ground truth
U = ttl(G,list(A,B,C),ms = c(1,2,3))@data

#hist(U); quantile(U,c(0,0.01,0.25,0.75,0.99,1))
#hist(inv.logit(U)); quantile(inv.logit(U),c(0,0.01,0.25,0.75,0.99,1))

### ts:binary tensor
ts = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
ts = as.tensor(array(ts,dim = c(d1,d2,d3)))@data


# ### the logLik of true U and data Y with noise -------------------------------------
# sum(log(inv.logit((2*ts-1)*U)))
# # p=plogis(U)
# # sum(Y*log(p)+(1-Y)*log(1-p))



upp = update_binary_raw(ts,c(r1,r2,r3),100)

U_est = ttl(upp$G,list(upp$A,upp$B,upp$C),ms = c(1,2,3))@data
### U_est: estimated ground truth
sum(log(inv.logit((2*ts-1)*U_est)))


upp$lglk

plot(upp$lglk)
plot(upp$lglk[upp$lglk!=0])



sum((U_est - U)^2)/prod(whole_shape)  ### MSE
plot(as.vector(U),as.vector(U_est));abline(0,1)   ## plot U vs U_est
plot(as.vector(inv.logit(U)),as.vector(inv.logit(U_est)));abline(0,1)  ## plot U vs U_est in logical scale





#plot(as.vector(U[,,8]),as.vector(U_est[,,8]));abline(0,1)
#plot(as.vector(inv.logit(U[,,8])),as.vector(inv.logit(U_est[,,8])));abline(0,1)




# #########---------------   identify  spatial location of outliers
# 
# which(as.vector(U_est > 3000))
# 
# U_est[4,2,5]           ####  spatial location of outliers
# as.vector(U_est)[1624]
# 
# U_est[10,2,5]
# as.vector(U_est)[1630]
# 
# 
# 
# sum(ts[1:15,1:5,2:8])
# 
# which(U_est > 3000)

























