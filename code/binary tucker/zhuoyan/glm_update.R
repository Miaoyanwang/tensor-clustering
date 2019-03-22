rm(list = ls())
library(rTensor)
library(pracma)
library(gtools)

library(rgl)
library("RColorBrewer")

glm_f = function(y,X)
{
  ml = glm(y ~ -1 + X, family = binomial(link = 'logit'),
           control = glm.control(maxit = 100, trace = TRUE))
  coe = ml$coefficients
  lglk = logLik(ml)
  return(list(coe,lglk))
}

#    form U = X*Beta
glm_mat = function(Y,X){
  beta_row = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  la = sapply(as.data.frame(Y),glm_f,X)
  re = t(matrix(unlist(la), nrow = p, byrow=T))
  beta = re[1:beta_row,]
  lglk = sum(re[beta_row + 1,])
  return(list(beta,lglk))
}

# ### test
# set.seed(37)
# X = matrix(rnorm(30,5,9),10,3)
# Y = matrix(rbinom(20,1,0.5),10,2)
# 
# re = glm_mat(Y,X)
# re[1]
# 
# # df = as.data.frame(Y)
# # la = sapply(as.data.frame(Y),glm_f,X)
# # ma = t(matrix(unlist(la), nrow=2, byrow=T))
# # colnames(df)
# # 
# ml = glm(Y[,2] ~ 0 + X, family = binomial(link = 'logit'))
# ml$coefficients
# ml$deviance
# 
# ####

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
  
  lglk = seq(1,4*Nsim)
  for(n in 1:Nsim){
  
    ###### update A
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    re = glm_mat(t(Y_1),t(G_BC1))
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
    re = glm_mat(t(Y_2),t(G_AC2))
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
    re = glm_mat(t(Y_3),t(G_AB3))
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
  }
  return(list(A = A,B = B,C = C,G = G,lglk = lglk))
}


up = update_binary(data,c(4,4,5),2)

####  function to plot tensor
# author: prof.Wang
marker = list(color = brewer.pal(9, "RdBu"))

### function to plot 3d array
plot_tensor=function(tensor){
  
  position=positionfun(dim(tensor))$position
  quan=c(quantile(tensor,(0:8)/8),max(tensor))
  col=tensor
  for(i in 1:9){
    col[(tensor>=quan[i])&(tensor<quan[i+1])]=marker$color[i]
  }
  
  plot3d(position[,1],position[,2],position[,3],col=col,alpha=0.3,size=5,xlab="",ylab="",zlab="")
}
## function to find the array index
positionfun=function(d){
  x=rep(1,d[2])%x%rep(1,d[3])%x%c(1:d[1])
  y=rep(1,d[3])%x%(c(1:d[2])%x%rep(1,d[1]))
  z=c(1:d[3])%x%rep(1,d[1])%x%rep(1,d[2])
  position=cbind(x,y,z)
  return(list("position"=position))
}

######  simulation

simu = function(Nsim,core_shape, whole_shape){
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  
  MSE = seq(Nsim)
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
  }
  
  return(MSE = MSE)
}
MSE = simu(10,c(3,4,5),c(20,30,40))
plot(MSE,type = 'b')


plot_tensor(U)
rgl.snapshot("Ui.png")
plot_tensor(U_est)
rgl.snapshot("U_t.png")

#################################  draft

set.seed(37)

A = simu(1,c(3,4,5),c(20,30,40))



set.seed(37)
G = as.tensor(array(data = rnorm(60),dim = c(3,4,5)))
A = randortho(20)[,1:3]
B = randortho(30)[,1:4]
C = randortho(40)[,1:5]

bl = sort(sample(seq(20),2,replace = F))
A[,] = 0
A[1:bl[1],1] = 1
A[(bl[1] + 1):bl[2],2] = 1
A[(bl[2] + 1):20,3] = 1

bl = sort(sample(seq(30),3,replace = F))
B[,] = 0
B[1:bl[1],1] = 1
B[(bl[1] + 1):bl[2],2] = 1
B[(bl[2] + 1):bl[3],3] = 1
B[(bl[3] + 1):30,4] = 1


bl = sample(seq(40),4,replace = F)
C[,] = 0
C[1:bl[1],1] = 1
C[(bl[1] + 1):bl[2],2] = 1
C[(bl[2] + 1):bl[3],3] = 1
C[(bl[3] + 1):bl[4],4] = 1
C[(bl[4] + 1):40,5] = 1



G@modes

U = ttl(G,list(A,B,C),ms = c(1,2,3))@data

Y = rbinom(20*30*40,1,prob = as.vector( 1/(1 + exp(-U)) ) )
Y = as.tensor(array(Y,dim = c(20,30,40)))@data


up = update_binary(Y,c(3,4,5),5)
plot(up$lglk)

plot_tensor(U)

kmA = kmeans(up$A,centers = 3)
kmB = kmeans(up$B,centers = 4)
kmC = kmeans(up$C,centers = 5)

mem_A = model.matrix(~ -1 + factor(kmA$cluster))
#rownames(mem_A) = con
mem_A = mem_A[,c(1,3,2)]

mem_B = model.matrix(~ -1 + factor(kmB$cluster))
#rownames(mem_B) = con
mem_B = mem_B[,c(2,1,4,3)]

mem_C = model.matrix(~ -1 + factor(kmC$cluster))
#rownames(mem_C) = relname
mem_C = mem_C[,c(2,3,4,1,5)]



U_est = ttl(up$G,list(mem_A,mem_B,mem_C),ms = c(1,2,3))@data

plot_tensor(U_est)


sum((ttl(up$G,list(up$A,up$B,up$C),ms = c(1,2,3))@data - U)^2)




