### upgrade function####
#rm(list = ls())
#setwd("/Users/March/Desktop/Tensor/code/Binarytensor/")

#packages 
library(R.matlab)
library(rTensor)
library(tensr)
library(tensor)
library(gtools)

#upgrading functions
mode_upgrade_glm = function(tensor,U_0,M1_0,M2_0,M3_0){
  #U_0,M1_0,M2_0,M3_0 are initial core tensor and factor matrices
  #M1_0 \in R^{d_1 \times k_1}
  
  #get the dimensions
  k1 = dim(U_0)[1]
  k2 = dim(U_0)[2]
  k3 = dim(U_0)[3]
  
  d1 = dim(M1_0)[1]
  d2 = dim(M2_0)[1]
  d3 = dim(M3_0)[1]
  
  #set likelihood 
  likelihood = rep(0,4)
  
  #find problem matrices
  mode_p1 = NULL
  mode_p2 = NULL
  mode_p3 = NULL
  mode_p4 = NULL

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
  if(is.array(U_0)==FALSE){
    U_0 = U_0@data #tensor class is a S4 class, need to use @ to get the data
  }
  
  if(is.array(U_0)==FALSE){
    print("not array")
  }
  
  #find the start coefficient
  start_u = rep(0,k1*k2*k3)
  a=1
  for (i in 1:k1) {
    for (j in 1:k2) {
      for (k in 1:k3) {
        start_u[a]= U_0[i,j,k]
        a = a+1
      }
    }
  }
  
  glforu = glm(y~0+x,family=binomial(link="logit"),control = list(maxit = 100))
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
  
  likelihood[1] = cal_likelihood(tensor,U,M1_0,M2_0,M3_0)
  
  if(likelihood[1] == -Inf){#if give inf then try another strategy
    glforu = glm(y~0+x,family=binomial(link="logit"),control = list(maxit = 100)
                 ,start = start_u)
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
    
    likelihood[1] = cal_likelihood(tensor,U,M1_0,M2_0,M3_0)
  }
  
  if(likelihood[1] == -Inf){ #after the strategy if the likelihood is still inf
    print("still inf")
    mode_p1 = list("U"=U,"M1"=M1_0,"M2" =M2_0,"M3"=M3_0)
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
    y = tensor[1,,i]
    for (j in 2:d1) {
      y = c(y,tensor[j,,i])
    }
    gm = glm(y~0+x3,family=binomial(link="logit"),control = list(maxit = 50))
    M3[i,]=gm$coefficients
  }

  likelihood[2] = cal_likelihood(tensor,U,M1_0,M2_0,M3)
  
  if(likelihood[2] == -Inf){
    M3 = matrix(0,nrow = d3,ncol = k3)
    for (i in 1:d3) {
      y = tensor[1,,i]
      for (j in 2:d1) {
        y = c(y,tensor[j,,i])
      }
      gm = glm(y~0+x3,family=binomial(link="logit"),control = list(maxit = 50)
               ,start = M3_0[i,])
      M3[i,]=gm$coefficients
    }
    
    likelihood[2] = cal_likelihood(tensor,U,M1_0,M2_0,M3)
  }
  
  if(likelihood[2] == -Inf){
    print("still inf")
    mode_p2 = list("U"=U,"M1"=M1_0,"M2" =M2_0,"M3"=M3)
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
  
  likelihood[3] = cal_likelihood(tensor,U,M1_0,M2,M3)
  
  if(likelihood[3] == -Inf){
    M2 = matrix(0,nrow = d2,ncol = k2)
    for (i in 1:d2) {
      y = tensor[1,i,]
      for (j in 2:d1) {
        y = c(y,tensor[j,i,])
      }
      gm = glm(y~0+t(x2),family=binomial(link="logit"),control = list(maxit = 50)
               ,start = M2_0[i,])
      M2[i,]=gm$coefficients
    }
    likelihood[3] = cal_likelihood(tensor,U,M1_0,M2,M3)
  }
  
  if(likelihood[3] == -Inf){
    print("still inf")
    mode_p3 = list("U"=U,"M1"=M1_0,"M2" =M2,"M3"=M3)
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
  
  likelihood[4] = cal_likelihood(tensor,U,M1,M2,M3)
  if(likelihood[4] == -Inf){
    M1 = matrix(0,nrow = d1,ncol = k1)

    for (i in 1:d1) {
      y = tensor[i,1,]
      for (j in 2:d2) {
        y = c(y,tensor[i,j,])
      }
      gm = glm(y~0+t(x1),family=binomial(link="logit"),control = list(maxit = 50)
               ,start = M1_0[i,])
      M1[i,]=gm$coefficients
    }
    
    likelihood[4] = cal_likelihood(tensor,U,M1,M2,M3)
  }
  if(likelihood[4] == -Inf){
    print("still inf")
    mode_p4 = list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3)
  }
  
  #,"mode_p" = list("p1" = mode_p1,"p2" = mode_p2,"p3" = mode_p3,"p4" = mode_p4)
  return(list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3,"likelihood" = likelihood
              ,"mode_p" = list("p1" = mode_p1,"p2" = mode_p2,"p3" = mode_p3,"p4" = mode_p4)))
}

cal_likelihood = function(tensor,U,M1,M2,M3){
  if(is.array(U) ==TRUE) U =as.tensor(U)
  T1 = ttm(U,M1,1)#must use class tensor
  T2 = ttm(T1,M2,2)
  est_tensor = ttm(T2,M3,3)
  #first_term = innerProd(as.tensor(tensor),est_tensor)
  #second_term = sum(log(1+exp(est_tensor@data)))
  #like = first_term-second_term
  like=sum(log(inv.logit((2*tensor-1)*est_tensor@data)))
  return(like)
}

###### below are the process I tried to find the problem #####
###### so below are messy and unnecessary contents #######
###### I even consider whether it is because of the overflow of the computer###
###### But initialization and correlation are still the main problems ####

# U = problem$p1$U
# M1=problem$p1$M1
# M2=problem$p1$M2
# M3=problem$p1$M3
# sum(exp(est_tensor@data))#this would be inf
# 
# library(Rmpfr)
# big_est = mpfr(est_tensor@data,precBits = 6)
# second_term = sum(log(1+exp(big_est)))
# second_term
# 
# ###another way to calculate the likelihood 
# p = exp(est_tensor@data)/(1+exp(est_tensor@data))
# p
# first_term1 = innerProd(as.tensor(tensor),as.tensor(p))
# #不太行 exp（e+13）真的会inf
# 
# ###check correlation###
# d1 = 14 
# d2 = 14
# d3 = 56
# k1 = 3
# k2 = 3
# k3 = 5
# y = rep(0,d1*d2*d3)
# x = matrix(0,nrow = d1*d2*d3,ncol = k1*k2*k3)
# 
# m=1
# 
# for (i in 1:d1) {
#   for (j in 1:d2) {
#     for (k in 1:d3) {
#       y[m] = tensor[i,j,k]
#       x[m,] =kronecker(kronecker(mode_t$M1[i,],mode_t$M2[j,]),mode_t$M3[k,]) 
#       m = m+1
#     }
#   }
# }
# dim(x)
# cor(x)#the correlation of x is not very big
# 
# glforu = glm(y~0+x,family=binomial(link="logit"),control = list(maxit = 100))
# logLik(glforu)
# 
# ######搞不定了 我要尝试initialize了
# #上一个u
# 
# start_u = rep(0,k1*k2*k3)
# n=1
# for (i in 1:k1) {
#   for (j in 1:k2) {
#     for (k in 1:k3) {
#       start_u[n] = mode_t$U[i,j,k] 
#       n=n+1
#     }
#   }
# }
# 
# glforu0 = glm(y~0+x,family=binomial(link="logit"),control = list(maxit = 100)
#               ,start = start_u)
# 
# u_new = matrix(0,nrow = k1*k2*k3)
# dim(u_new) = c(k1,k2,k3)
# n=1
# for (i in 1:k1) {
#   for (j in 1:k2) {
#     for (k in 1:k3) {
#       u_new[i,j,k] =  glforu0$coefficients[n]
#       n=n+1
#     }
#   }
# }
# 
# cal_likelihood(tensor,U , mode_t$M1,mode_t$M2,mode_t$M3)
# cal_likelihood(tensor,u_new , M1,M2,M3)#是initialize 的问题
