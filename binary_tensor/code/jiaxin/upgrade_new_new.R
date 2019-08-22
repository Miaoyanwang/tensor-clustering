#**Setting for the workpath and clean the environment**

#rm(list = ls())
#setwd("/Users/March/Desktop/Tensor/code/Binarytensor/")


#**Set a global variable to control the log-likelihood **

global_max = 1e+14


#**Set the packages **

if(!require("R.matlab")){
  install.packages("R.matlab")
  stopifnot(require("R.matlab"))
}
if(!require("rTensor")){
  install.packages("rTensor")
  stopifnot(require("rTensor"))
}
if(!require("gtools")){
  install.packages("gtools")
  stopifnot(require("gtools"))
}


#**The function to calculate the log-likelihood, in which the input U should be a tensor**

cal_likelihood = function(tensor,U,M1,M2,M3){
  #U must be tensor 
  #tensor must be tensor
  est_tensor = ttl(U,list_mat = list(M1,M2,M3),ms = c(1,2,3))
  
  temp=(2*tensor@data-1)*est_tensor@data
  like=sum(log(inv.logit(temp[temp>=0])))+sum(temp[temp<0]-log(1+exp(temp[temp<0])))
  
  return(like)
}


#**Normalize the factor matrices**

normalize_tucker = function(U,M1,M2,M3,k){#k is a vector of rank
  #U must be tensor
  est_tensor = ttl(U,list_mat = list(M1,M2,M3),ms = c(1,2,3))
  
  normal = tucker(est_tensor,ranks = k)
  
  U = normal$Z
  M1 = as.matrix(normal$U[[1]])
  M2 = as.matrix(normal$U[[2]])
  M3 = as.matrix(normal$U[[3]])
  
  return(list("U" = U,"M1" = M1,"M2"=M2,"M3"=M3))
}



#**The function to calculate the coefficients given by three different methods 
#of glm and choose the best one with max log-likelihood**

glm_modify=function(y,x,start){
  ## initial coefficent
  ini_loglik=sum(log(inv.logit((2*y-1)*(x%*%start))))
  
  ## Option 1: glm fittig with default initilization
  fit1 = glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50))
  
  ## Option 2: glm with user specified initilization
  fit2= glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50),start=start)
  
  ## report the result whichever gives the highest likelihood
  if(max(logLik(fit1),logLik(fit2))<ini_loglik) return (start)
  
  else if(logLik(fit1)>logLik(fit2)) return(coef(fit1))
  
  else return(coef(fit2))
}


#**The pipline function to upgrade the core tensor and three factor matrices**

upgrade_glm_tucker = function(tensor,U,M1,M2,M3){
  #tensor must be a tensor
  #U should be a tensor
  #M1,M2,M3 are matrices
  
  ### get dimension
  k1 = dim(U)[1]
  k2 = dim(U)[2]
  k3 = dim(U)[3]
  
  d1 = dim(M1)[1] #M1 has k1 columns
  d2 = dim(M2)[1]
  d3 = dim(M3)[1]
  
  #set likelihood 
  likelihood = rep(0,4)
  
  #find problem matrices
  mode_p1 = NULL
  mode_p2 = NULL
  mode_p3 = NULL
  mode_p4 = NULL
  
  ### upgrading U
  #find x
  x = matrix(0,nrow = d1*d2*d3,ncol = k1*k2*k3)
  m=1
  for (k in 1:d3) {
    for (j in 1:d2) {
      for (i in 1:d1) {
        #y[m] = tensor[i,j,k]
        x[m,] =kronecker(kronecker(M3[k,],M2[j,]),M1[i,]) 
        m = m+1
      }
    }
  }
  
  coef_u = glm_modify(as.vector(tensor@data),x,as.vector(U@data))
  U = as.tensor(array(coef_u,dim = c(k1,k2,k3)))
  likelihood[1] = cal_likelihood(tensor,U,M1,M2,M3)
  
  if(abs(likelihood[1])  >= global_max ){ #after modify the glm still gives -inf
    print("still inf")
    mode_p1 = list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3)
  }
  
  normal = normalize_tucker(U,M1,M2,M3,c(k1,k2,k3))
  U = normal$U#return a tensor
  M1 = normal$M1
  M2 = normal$M2
  M3 = normal$M3
  
  ######  upgrading M3
  y3 = unfold(tensor,row_idx = c(2,1),col_idx = 3)@data
  
  UM1M2 = ttl(U,list_mat = list(M1,M2),c(1,2))
  x3 = unfold(UM1M2,row_idx = c(2,1),col_idx = 3)@data
  
  M3 = mapply(glm_modify,y = as.data.frame(y3),start = as.data.frame(t(M3))
              ,MoreArgs = list(x = x3))
  M3 = t(M3)
  
  likelihood[2] = cal_likelihood(tensor,U,M1,M2,M3)
  
  if(abs(likelihood[2])  >= global_max ){ #after modify the glm still gives -inf
    print("still inf")
    mode_p2 = list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3)
  }
  
  normal = normalize_tucker(U,M1,M2,M3,c(k1,k2,k3))
  U = normal$U#return a tensor
  M1 = normal$M1
  M2 = normal$M2
  M3 = normal$M3
  
  ####### upgrading M2
  y2 = unfold(tensor,row_idx = c(3,1),col_idx = 2)@data
  
  UM1M3 = ttl(U,list_mat = list(M1,M3),c(1,3))
  x2 = unfold(UM1M3,row_idx = c(3,1),col_idx = 2)@data
  
  M2 = mapply(glm_modify,y = as.data.frame(y2),start = as.data.frame(t(M2))
              ,MoreArgs = list(x = x2))
  M2 = t(M2)
  
  likelihood[3] = cal_likelihood(tensor,U,M1,M2,M3)
  
  if(abs(likelihood[3])  >= global_max ){ #after modify the glm still gives -inf
    print("still inf")
    mode_p3 = list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3)
  }
  
  normal = normalize_tucker(U,M1,M2,M3,c(k1,k2,k3))
  U = normal$U#return a tensor
  M1 = normal$M1
  M2 = normal$M2
  M3 = normal$M3
  
  ###### upgrading M1
  y1 = unfold(tensor,row_idx = c(3,2),col_idx = 1)@data
  
  UM2M3 = ttl(U,list_mat = list(M2,M3),c(2,3))
  x1 = unfold(UM2M3,row_idx = c(3,2),col_idx = 1)@data
  
  M1 = mapply(glm_modify,y = as.data.frame(y1),start = as.data.frame(t(M1))
              ,MoreArgs = list(x = x1))
  M1 = t(M1)
  
  likelihood[4] = cal_likelihood(tensor,U,M1,M2,M3)
  
  if(abs(likelihood[4])  >= global_max ){ #after modify the glm still gives -inf
    print("still inf")
    mode_p4 = list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3)
  }
  
  normal = normalize_tucker(U,M1,M2,M3,c(k1,k2,k3))
  U = normal$U#return a tensor
  M1 = normal$M1
  M2 = normal$M2
  M3 = normal$M3
  
  return(list("U"=U,"M1"=M1,"M2" =M2,"M3"=M3,"likelihood" = likelihood
              ,"mode_p" = list("p1" = mode_p1,"p2" = mode_p2
                               ,"p3" = mode_p3,"p4" = mode_p4)))
  
}
