rm(list = ls())
library(rTensor)
library(pracma)
library(gtools)
library(rgl)
library("RColorBrewer")
###  update glm without any contrain
#######################
# author: prof. Wang
glm_modify=function(y,x,start){
  
  ## initial coefficent
  ini_loglik=sum(log(inv.logit((2*y-1)*(x%*%start))))
  
  ## Option 1: glm fittig with default initilization
  fit1 = glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F))
  
  ## Option 2: glm with user specified initilization
  fit2= glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50,trace = F),start=start)
  
  ## report the result whichever gives the highest likelihood
  if(max(logLik(fit1),logLik(fit2))<ini_loglik) return (list(start, ini_loglik))
  else if(logLik(fit1)>logLik(fit2)) return(list(coef(fit1), logLik(fit1)))
  else return(list(coef(fit2), logLik(fit2)))
  
}

### matrix form glm
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

###########
########  some test
set.seed(37)
X = matrix(rnorm(30,5,9),10,3)
Y = matrix(rbinom(20,1,0.5),10,2)
start = matrix(rnorm(6,1,0.5),3,2)



glm_f(Y[,1],start = start[,1],X)
glm_f(Y[,2],start = start[,2],X)


re = glm_mat(Y,start = start,X)
re[1]
re[[2]]
####

################################################################################################
############  update A,B with constrain
###########################################  define loss
loss = function(beta,z,X,lambda,y,D){
  p=plogis(X %*% beta)
  L=-z*log(p)-(1-z)*log(1-p)
  LwR2=sum(L)+lambda*sum((y - D%*%beta)^2)
  return(c(LwR2))
}

loss_gr <- function(beta,z,X,lambda,y,D){
  p=plogis(X %*% beta)
  v=t(X) %*% (p - z)
  regu_loss = 2*lambda*t(D)%*%(D%*%beta - y)
  return(c(v)+regu_loss)
}

######## a little test
library(mlbench)
set.seed(37)
d=mlbench.2dnormals(100,2)
X=d$x
z=ifelse(d$classes==1,1,0)
lambda=1
y = rnorm(20)
D = matrix(rnorm(40),20,2)

optim(par = runif(2),loss,loss_gr,z = z,X = X,lambda = lambda,y = y,D = D
      ,method="BFGS")
##################################### define opt  10 steps each iteration
opt = function(beta_initial, z = z,y = y,X = X,lambda = lambda,D = D
               ,method="BFGS"){
  res = optim(par = beta_initial,loss,loss_gr,z = z,X = X,lambda = lambda,y = y,D = D
              ,method="BFGS",control = list(maxit = 50,trace = F))
  coe = res$par
  loss = res$value
  return(list(coe,loss))
}
### a little test
res = opt(beta_initial = c(0.2,5),z = z,y = y,X = X,lambda = lambda,D = D
    ,method="BFGS")
unlist(res)
##################################################  matrix form GLM
#    form U = X*Beta
glm_mat_regu = function(Beta,Z,Y,X,lambda,D){
  R = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  ma = mapply(opt, as.data.frame(Beta),as.data.frame(Z), as.data.frame(Y),MoreArgs = list(X,lambda,D))
  re = t(matrix(unlist(ma), nrow = p, byrow=T))
  beta = re[1:R,]
  loss = sum(re[R + 1,])
  return(list(beta,loss))
}
## a little test
Y = cbind(y,y)
Z = cbind(z,z)
Beta = matrix(c(0.2,5,0.2,5),2,2)
glm_mat_regu(Beta,Z,Y,X,lambda,D)

################################################################################################
############  update W,U with regularizer
### standardized
library(glmnet)
glmnet_regu = function(y,X,lambda = rep(0.005))  
{
  re = glmnet(X,y,alpha = 1,lambda = lambda,intercept = FALSE,standardize = FALSE)
  coe = re$beta
  coe = coe@x    ### extrace coefficient
  #lglk = re[[2]]
  return(coe)
}

#    form U = X*Beta
glmnet_mat_regu = function(Y,X,lambda = rep(0.005)){
  R = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  ma = mapply(glmnet_regu, as.data.frame(Y),MoreArgs = list(X,lambda))
  re = t(matrix(unlist(ma), nrow = p, byrow=T))
  beta = re[1:R,]
  #lglk = sum(re[R + 1,])
  return(beta)
}

## a little test
set.seed(37)
X=matrix(rnorm(100*5),100,5)
y=rnorm(100)
fit1=glmnet(X,y,lambda = 0.005,intercept = FALSE)
coef(fit1)

fit1 = glmnet_regu(y,X,lambda = 0.005)
fit1

Y = cbind(y,y)
fit2 = glmnet_mat_regu(Y,X,lambda = 0.005)
fit2






















