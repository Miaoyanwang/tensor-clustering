### upgrade function####
#rm(list = ls())
#setwd("/Users/March/Desktop/Tensor/code/Binarytensor/")

#packages 
library(R.matlab)
library(rTensor)
library(tensr)
library(tensor)

#upgrading functions
mode_upgrade_glm = function(tensor,U_0,M1_0,M2_0,M3_0){
  #U_0,M1_0,M2_0,M3_0 are initial core tensor and factor matrices
  #M1_0 \in R^{d_1 \times k_1}
  #class of tensor should be array
  
  #get the dimensions
  k1 = dim(U_0)[1]
  k2 = dim(U_0)[2]
  k3 = dim(U_0)[3]
  
  d1 = dim(M1_0)[1]
  d2 = dim(M2_0)[1]
  d3 = dim(M3_0)[1]
  
  #upgrade core tensor
  
  #first find y of glm #find predictor
  y = rep(0,d1*d2*d3)
  x = matrix(0,nrow = d1*d2*d3,ncol = k1*k2*k3)
  m=1
  for (i in 1:d1) {
    for (j in 1:d2) {
      for (k in 1:d3) {
        y[m] = tensor[i,j,k]
        x[m,] =kronecker(kronecker(M1_0[i,],M2_0[j,]),M3_0[k,]) 
        m = m+1
      }
    }
  }
  
  #y0 = unfold(as.tensor(tensor),row_idx = 1,col_idx = c(3,2))
  #row indices from M1, column indices from M3 and M1
  #y = as.vector(t(y0))
  
  # what y looks like
  #[1,1,1],...,[1,1,d3],[1,2,1],...,[1,2,d3],...,[1,d2,d3],
  #[2,1,1],...,[2,1,d3],[1,2,1],...,[1,2,d3],...,[2,d3,d3],...
  
  glforu = glm(y~0+x,family=binomial(link="logit"),control = list(maxit = 50))
  glforu_coef = glforu$coefficients
  U = matrix(0,nrow = k1*k2*k3)
  dim(U) = c(k1,k2,k3)
  n=1
  for (i in 1:k1) {
    for (j in 1:k2) {
      for (k in 1:k3) {
        U[i,j,k] =  glforu$coefficients[n]
        n=n+1
      }
    }
  }
  
  ###upgrade M_3###
  
  UM1 = amprod(U, M1_0, 1)
  UM1M2 = amprod(UM1,M2_0,2)
  
  x3 = UM1M2[1,,]
  for (a in 2:d1) {
    x3 = rbind(x3,UM1M2[a,,])
  }
  
  M3 = matrix(0,nrow = d3,ncol = k3) 
  for (i in 1:d3) {
    y = tensor[1,,3]
    for (j in 2:d1) {
      y = c(y,tensor[j,,i])
    }
    gm = glm(y~0+x3,family=binomial(link="logit"),control = list(maxit = 50))
    M3[i,]=gm$coefficients
  }

  ###upgrade M_2###
  
  UM1M3 = amprod(UM1,M3,3)
  
  x2 = UM1M3[1,,]
  for (a in 2:d1) {
    x2 = cbind(x2,UM1M3[a,,])
  }
  
  M2 = matrix(0,nrow = d2,ncol = k2) 
  for (i in 1:d2) {
    y = tensor[1,i,]
    for (j in 2:d1) {
      y = c(y,tensor[j,i,])
    }
    gm = glm(y~0+t(x2),family=binomial(link="logit"),control = list(maxit = 50))
    M2[i,]=gm$coefficients
  }
  
  ###upgrade M_1###
  
  UM2 = amprod(U,M2,2)
  UM2M3 = amprod(UM2,M3,3)
  
  x1 = UM2M3[,1,]
  for (a in 2:d2) {
    x1 = cbind(x1,UM2M3[,a,])
  }
  
  M1 = matrix(0,nrow = d1,ncol = k1) 
  
  for (i in 1:d1) {
    y = tensor[i,1,]
    for (j in 2:d2) {
      y = c(y,tensor[i,j,])
    }
    gm = glm(y~0+t(x1),family=binomial(link="logit"),control = list(maxit = 50))
    M1[i,]=gm$coefficients
  }
  
  return(list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3))
}

cal_likelihood = function(tensor,U,M1,M2,M3){#they should be array or matrix
  T1 = amprod(U,M1,1)
  T2 = amprod(T1,M2,2)
  est_tensor = amprod(T2,M3,3)
  first_term = innerProd(as.tensor(tensor),as.tensor(est_tensor))
  second_term = sum(log(1+exp(est_tensor)))
  like = first_term-second_term
  
  return(like)
}
