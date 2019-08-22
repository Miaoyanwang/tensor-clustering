########pipline######
########### binary tensor factorization and cluster(kmeans)#########
#rm(list = ls())
#setwd("/Users/March/Desktop/Tensor/code/Binarytensor/")

source("upgrade_new.R")

######packages#######
library(R.matlab)
library(rTensor)
library(tensr)

binary_tensor_fctr_km = function(tensor,k,maxiter){
  #k is number for cp factorization 
  #and k can be a vector is we use trucker decomposition
  ##########replace NaN by 0######
  tensor[is.na(tensor)] = 0
  
  #########find the initial matrices#########
  #if use cp
  #########take the elements as continuous variables#######
  #######zoom y to 10*(2y-1)#######
  # tensor_i = tensor
  # tensor_i = as.tensor(10*(2*tensor_i-1))
  
  ######The model would be#######
  ######logit(p_ijk) = C x M1 x M2 x M3
  ######initialize C,M1,M2,M3#######

  # inital_dep = cp(tnsr = tensor_i,num_components = k)
  # M1_0 = as.matrix(inital_dep$U[[1]])
  # M2_0 = as.matrix(inital_dep$U[[2]])
  # M3_0 = as.matrix(inital_dep$U[[3]])
  # ##initial C would be identity tensor###
  # C_0 = matrix(0,nrow = k^3)
  # dim(C_0) = c(k,k,k)
  # for (i in 1:k) {####special initialize for cp factorization
  #   C_0[i,i,i] =1
  # }
  # 
  
  #here use tucker
  inital_tucker = tucker(as.tensor(10*(2*tensor-1)),ranks = k)
  C_0 = inital_tucker$Z
  M1_0 = as.matrix(inital_tucker$U[[1]])
  M2_0 = as.matrix(inital_tucker$U[[2]])
  M3_0 = as.matrix(inital_tucker$U[[3]])
  likelihood = cal_likelihood(tensor,C_0,M1_0,M2_0,M3_0)
  likelihood_normal = cal_likelihood(tensor,C_0,M1_0,M2_0,M3_0)
  modes = list("U"=C_0,"M1"=M1_0,"M2" =M2_0,"M3"=M3_0)
  
  mode_problem = NULL
  for (iter in 2:maxiter) {
    modes_new = mode_upgrade_glm(tensor,modes$U,modes$M1,modes$M2,modes$M3)
    likelihood = c(likelihood,modes_new$likelihood)
    likelihood_normal = c(likelihood_normal,modes_new$likelihood_normal)
    print(iter)
    if(likelihood[(1+4*(iter-1))]<=likelihood[(4*(iter-2)+1)]){
      #compare the full upgrade likelihood
      mode_problem = modes_new
      #record the problem matrices and core tensor
      #and the problem matrices and core tensor for each step can be 
      #found by mode_problem$mode_p
      #which record the problem matries and tensor if the likelihood is inf
      break
    }
    modes = modes_new
  }
  
  M1_km = kmeans(modes$M1,centers = k[1])#if use trucker, there would be k1
  M2_km = kmeans(modes$M2,centers = k[2])
  M3_km = kmeans(modes$M3,centers = k[3])
  
  return(list("modes"=modes,"M1_CL" = M1_km$cluster,"M2_CL" = M2_km$cluster
              ,"M3_CL" = M3_km$cluster,"likelihood"=likelihood
              ,"mode_problem" = mode_problem
              ,"likelihood_normal" = likelihood_normal))
}

dnoations = readMat("dnations.mat")
dnoations_tr = dnoations$R
test_222  =binary_tensor_fctr_km(dnoations_tr,c(2,2,2),50)
test_223  =binary_tensor_fctr_km(dnoations_tr,c(2,2,3),50)
test_333  =binary_tensor_fctr_km(dnoations_tr,c(3,3,3),50)
test = binary_tensor_fctr_km(dnoations_tr,c(3,3,4),50)

# tensor = dnoations_tr
# tensor[is.na(tensor)] = 0
# inital_tucker = tucker(as.tensor(tensor),ranks = c(4,4,4))
# C_0 = inital_tucker$Z
# M1_0 = as.matrix(inital_tucker$U[[1]])
# M2_0 = as.matrix(inital_tucker$U[[2]])
# M3_0 = as.matrix(inital_tucker$U[[3]])
# cal_likelihood(tensor,C_0,M1_0,M2_0,M3_0)

###why likelihood would drop to negative infinity???
problem = test$mode_problem$mode_p$p1#p3 means the third upgrading step sth wrong
cor(problem$M3)#the third step corresponds to upgrade M2 or B in our document

mode_t = test$modes
#mode_t 是出现inf情况前那次upgrading最后的更新结果 
#在这里就是第三次带入mode_upgrade_glm 的那些矩阵和tensor

test1 = mode_upgrade_glm(tensor,mode_t$U,mode_t$M1,mode_t$M2,mode_t$M3)
#test1 就是第三次的更新结果
#以下在看更新U时glm出了什么问题

k1 =3
k2 = 3
k3 =4
d1 = 14
d2 = 14
d3 = 56

tensor = dnoations_tr
tensor[is.na(tensor)] = 0

y = rep(0,d1*d2*d3)
x = matrix(0,nrow = d1*d2*d3,ncol = 3*3*4)

#scale 好像会改变数值的结果
#不过这里可以看到其实这些矩阵的数值没有很奇怪
# M1_0 = scale(mode_t$M1,scale = F)
# M2_0 = scale(mode_t$M2,scale = F)
# M3_0 = scale(mode_t$M3,scale = T)
M1_0 = mode_t$M1
M2_0 = mode_t$M2
M3_0 = mode_t$M3
cor(M1_0)#这三个矩阵的列也都近似正交
cor(M2_0)
cor(M3_0)

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
cor(x)#correlation 都很小
head(x)

start_u = rep(0,k1*k2*k3)
a=1
for (i in 1:k1) {
  for (j in 1:k2) {
    for (k in 1:k3) {
      start_u[a]= mode_t$U[i,j,k]#以之前的U的值作为glm搜索的起始值
      a = a+1
    }
  }
}

glforu = glm(y~0+x,family=binomial(link="logit"),control = list(maxit = 200)
             ,start = start_u )#报错了 并且回归系数的绝对值都很大

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

U
mode_t$U#更新前的u的数值还是很正常
cal_likelihood(tensor,U,mode_t$M1,mode_t$M2,mode_t$M3)


###如果更新m3的时候出了问题 检验更新m3的predictor是不是相关

# UM1 = amprod(problem$U, problem$M1, 1)
# UM1M2 = amprod(UM1,problem$M2,2)
# 
# x3 = UM1M2[1,,]
# for (a in 2:d1) {
#   x3 = rbind(x3,UM1M2[a,,])
# }
# cor(x3)
# dim(x3)
# cal_likelihood(tensor,problem$U,problem$M1,problem$M2,M3_y)
# 
# M3_y = matrix(0,nrow = d3,ncol = 4)
# for (i in 1:d3) {
#   y = tensor[1,,i]
#   for (j in 2:d1) {
#     y = c(y,tensor[j,,i])
#   }
#   gm = glm(y~0+x3,family=binomial(link="logit"),control = list(maxit = 50)
#            ,start = test$modes$M3[i,])
#   M3_y[i,]=gm$coefficients
# }

###如果更新m2的时候出了问题 检验更新m2的predictor是不是相关

UM1 = amprod(problem$U,problem$M1,1)
UM1M3 = amprod(UM1,problem$M3,3)

x2 = UM1M3[1,,]
for (a in 2:d1) {
  x2 = cbind(x2,UM1M3[a,,])
}
cor(t(x2)) 

#要转置因为要把tensor切片成矩阵的话 
#r会默认mode数小的是行，mode数大的是列



