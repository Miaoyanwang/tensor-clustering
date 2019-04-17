###### Semi-supervised binary tensor decomposition ######
###### pipline #####

rm(list = ls())
setwd("/Users/March/Desktop/Tensor/code/Binarytensor/semi-supervised/")

source("upgrade_covariate.R")

dnoations = readMat("/Users/March/Desktop/Tensor/code/Binarytensor/data/dnations.mat")

if(!require(pracma)){
  install.packages("pracma")
  stopifnot(require(pracma))
}

semi_supervised_binarydc_ttm = function(tensor,A,k,maxiter){
  #tensor is a binary tensor with dim d1 x d2 x d3
  #A is a covariate matrix with dim d1 x p
  #k is the rank of core tensor ,k = c(k1 , k2 , k3) 
  #maxiter is the max iteration approved
  
  #tensor should be an array
  if(is.array(tensor) == F){
    print("tensor should be an array")
    return()
  }
  
  #replace NaN to 0
  tensor[is.na(tensor)] = 0
  A[is.na(A)] = 0
  #set tensor as a tensor class
  tensor = as.tensor(tensor)
  
  #initializiation
  #y[,j,k] = X \beta[,j,k]
  
  y_initial = unfold(tensor,row_idx = 1,col_idx =c(2,3))

  beta = mapply(glm_modify,y = as.data.frame(y_initial@data),
                start = as.data.frame(matrix(0,nrow = p,ncol = dim(tensor)[2]*dim(tensor)[3])),
                MoreArgs = list(x = A))
  
  B = fold(beta, row_idx = 1 , col_idx = c(2,3),
           modes = c(dim(A)[2],dim(tensor)[2],dim(tensor)[3]))
  
  initial_tucker = tucker(B, ranks =  k)
  C = initial_tucker$Z
  N1 = as.matrix(initial_tucker$U[[1]])
  N2 = as.matrix(initial_tucker$U[[2]])
  N3 = as.matrix(initial_tucker$U[[3]])
  
  likelihood = cal_likelihood(tensor,C,N1,N2,N3,A)#the first likelihood
  modes = list("U" = C,"M1"=N1,"M2" =N2,"M3"=N3)
  
  for (iter in 2:maxiter) {
    modes_new = upgrade_glm_covariate(tensor,modes$U,modes$M1,modes$M2,modes$M3,A)
    likelihood = c(likelihood,modes_new$likelihood)
    print(iter)
    if((likelihood[(1+4*(iter-1))]-likelihood[(4*(iter-2)+1)])<= 0.01){
      #compare the full upgrade likelihood
      modes = modes_new
      break
    }
    modes = modes_new
  }
  
  return(list("modes"=modes,"likelihood"=likelihood ))
}


#Simulation

k1 = 3
k2 = 3
k3 = 3
k = c(k1,k2,k3)

p = 5

d1 = 20
d2 = 20
d3 = 20

sim = 10
MSE = rep(0,sim)
like_list = list()
like_orgin = list()
for (i in 1:sim) {
  #generate covariate X
  #set.seed(0426)
  #X = matrix(rnorm(p*d1,sd = 15) , nrow = d1 , ncol = p)
  X =matrix(rbinom(100,size = 1, prob = 0.5),20,5)
  #X = diag(1,nrow = d1,ncol = p)
  N1 = randortho(p,type = "orthonormal")[,1:k1]
  N2 = randortho(d2,type = "orthonormal")[,1:k2]
  N3 = randortho(d3,type = "orthonormal")[,1:k3]
  C = as.tensor(array(rnorm(k1*k2*k3,sd = 20),dim = c(k1,k2,k3)))
  true_B = ttl(C, list_mat = list(N1,N2,N3),ms = c(1,2,3))
  
  true_prob = ttm(true_B,X,m = 1)@data
  
  tsr = rbinom(20*20*20,1,prob = as.vector( 1/(1 + exp(-true_prob)) ) )
  tsr = as.tensor(array(tsr,dim = c(20,20,20)))@data
  
  test = semi_supervised_binarydc_ttm(tsr,X,k,20)
  like_list[[i]] = test$likelihood
  #plot(test$likelihood)
  like_orgin[[i]] = cal_likelihood(as.tensor(tsr),C,N1,N2,N3,X)
  modes = test$modes
  hat_B = ttl(modes$U,list_mat = list(modes$M1, modes$M2,modes$M3),ms = c(1,2,3))
  MSE[i] = mean((true_B@data - hat_B@data)^2)
}


par(mfrow = c(2,2))
plot(like_list[[10]])
abline(h = like_orgin[[10]],col = 2)
plot(like_list[[3]])
abline(h = like_orgin[[3]],col = 2)
plot(like_list[[6]])
abline(h = like_orgin[[6]],col = 2)
plot(like_list[[5]])
abline(h = like_orgin[[5]],col = 2)

