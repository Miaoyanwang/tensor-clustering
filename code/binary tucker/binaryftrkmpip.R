########pipline######
########### binary tensor factorization and cluster(kmeans)#########
#rm(list = ls())
#setwd("/Users/March/Desktop/Tensor/code/Binarytensor/")

source("upgrade_new.R")

######packages#######
library(R.matlab)
library(rTensor)
library(tensorr)

binary_tensor_fctr_km = function(tensor,k,maxiter){#k is for cp factorization 
  #and k can be a vector is we use trucker decomposition
  ##########replace NaN by 0######
  tensor[is.na(tensor)] = 0
  
  #########find the initial matrices#########
  #########take the elements as continuous variables#######
  #######zoom y to 10*(2y-1)#######
  tensor_i = tensor
  tensor_i = as.tensor(10*(2*tensor_i-1))
  
  ######The model would be#######
  ######logit(p_ijk) = C x M1 x M2 x M3
  ######initialize C,M1,M2,M3#######
  ###here use cp factorization#####
  inital_dep = cp(tnsr = tensor_i,num_components = k)
  M1_0 = as.matrix(inital_dep$U[[1]])
  M2_0 = as.matrix(inital_dep$U[[2]])
  M3_0 = as.matrix(inital_dep$U[[3]])
  ##initial C would be identity tensor###
  C_0 = matrix(0,nrow = k^3)
  dim(C_0) = c(k,k,k)
  for (i in 1:k) {####special initialize for cp factorization
    C_0[i,i,i] =1
  }
  
  likelihood = rep(0,maxiter)
  
  likelihood[1] = cal_likelihood(tensor,C_0,M1_0,M2_0,M3_0)
  modes = list("U"=C_0,"M1"=M1_0,"M2" =M2_0,"M3"=M3_0)
  
  for (iter in 2:maxiter) {
    modes_new = mode_upgrade_glm(tensor,modes$U,modes$M1,modes$M2,modes$M3)
    likelihood[iter] = cal_likelihood(tensor,modes_new$U,modes_new$M1
                                      ,modes_new$M2,modes_new$M3)
    print(iter)
    if(likelihood[iter]<=likelihood[iter-1]){
      likelihood = likelihood[1:iter]
      break
    }
    
    modes = modes_new
  }
  
  M1_km = kmeans(modes$M1,centers = k)#if use trucker, there would be k1
  M2_km = kmeans(modes$M2,centers = k)
  M3_km = kmeans(modes$M3,centers = k)
  
  return(list("modes"=modes,"M1_CL" = M1_km$cluster,"M2_CL" = M2_km$cluster
              ,"M3_CL" = M3_km$cluster,"likelihood"=likelihood))
}

dnoations = readMat("dnations.mat")
dnoations_tr = dnoations$R
test  =binary_tensor_fctr_km(dnoations_tr,5,10)


