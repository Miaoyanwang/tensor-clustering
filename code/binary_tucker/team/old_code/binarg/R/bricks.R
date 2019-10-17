library(rTensor)
library(pracma)
library(gtools)
library(rgl)
library("RColorBrewer")

loss = function(beta,y,X,lambda,alpha){
  U = X %*% beta
  p=plogis(X %*% beta)
  p[p==1] = p[p==1] - 0.0001
  p[p==0] = p[p==0] + 0.0001
  # if U is too large/small exp(U)/(1+exp(U)) would give value 1/0
  # that leads to the penalty function can not be intialized
  L=-y*log(p)-(1-y)*log(1-p)
  L2 =sum(L) - lambda * sum(log(1 - (U / alpha)^2))  ## object function
  #L2 =sum(L) + lambda * sum(log((U / alpha)^2))  ## object function
  return(c(L2))
}

## This part follows Sec 2.1.2 to compute Gradient of object function
loss_gr <- function(beta,y,X,lambda,alpha){
  U = X %*% beta
  p=plogis(X %*% beta)
  L_g=t(X) %*% (p - y)
  penal_g = 2 * lambda * t(X) %*% (U/(alpha^2 - U^2))
  #penal_g = 2 * lambda * t(X) %*% (1/U)
  return(c(L_g) + c(penal_g))
}

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

###########---------  GLM on two modes
##---  This section follows 1.3 of [Algorithm: Semi-Supervised Binary Tensor Factorization]
## add the linear model option for initilization
##---  function version
glm_two = function(Y, X1, X2, ini = FALSE, start, linear=FALSE, lm = FALSE # used for supervised scale down
                   ){ ## Y_size = m * n
  #X2 = t(X2)
  # logit(E(Y)) = X1 %*% coe %*% X2
  m = dim(Y)[1] ; n = dim(Y)[2]
  q1 = dim(X1)[2] ; q2 = dim(X2)[1]

  N_long = kronecker_list(list(t(X2),X1))

  if(ini == TRUE){
    coe_start = rnorm(q1*q2)
  }
  else {coe_start = as.vector(start)}

  if(linear==TRUE){
    mod_re = lm(10*(2*as.vector(Y)-1)~-1+N_long)
    coe = matrix(mod_re[[1]], nrow = q1, ncol = q2)
    lglk= logLik(mod_re)
    return(list(coe = coe, lglk = lglk))

  }

  if(lm==TRUE){
    mod_re = lm(as.vector(Y)~-1+N_long)
    coe = matrix(mod_re[[1]], nrow = q1, ncol = q2)
    lglk= logLik(mod_re)
    return(list(coe = coe, lglk = lglk))

  }

  mod_re = glm_modify(as.vector(Y), N_long, coe_start)
  coe = mod_re[[1]]
  coe = matrix(coe, nrow = q1, ncol = q2)
  lglk= mod_re[[2]]
  return(list(coe = coe, lglk = lglk))
}

#####---- This function is a parallel version of GLM on two modes, used for initialization
## add the linear model option for initilization
glm_two_mat = function(Y, X1, X2, ini = TRUE, start = NULL,linear, lm = FALSE){
  Yl = lapply(seq(dim(Y)[3]), function(x) Y[ , , x])
  re = lapply(Yl, glm_two, X1, X2, ini = ini,linear=linear, lm = lm)

  coe = lapply(seq(length(re)), function(x) re[[x]]$coe) ## extract coe
  coe = array(unlist(coe), dim = c(dim(X1)[2],dim(X2)[1],dim(Y)[3]))  ## form coe
  return(coe)
}

#####---- This is the function used for generating data through different distribution
#         of core tensor in unsupervised setting

gene_data_un = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),dis,
gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1, dup){ 
    #dis can be "gaussian" and "uniform"
    
    d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
    r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
    ####-------- generate data
    set.seed(seed)  # 24 # 37  #  347
    A = as.matrix(randortho(d1)[,1:r1])    ## factor matrix
    B = as.matrix(randortho(d2)[,1:r2])     ## factor matrix
    C = as.matrix(randortho(d3)[,1:r3])     ## factor matrix
    
    
    ### G: core tensor
    if(dis == "gaussian"){
        G = as.tensor(array(data = rnorm(r1*r2*r3,mean = gs_mean,sd = gs_sd),dim = core_shape))
    }
    else if(dis == "uniform"){
        G = as.tensor(array(data = runif(r1*r2*r3,min = unf_a,max = unf_b),dim = core_shape))
    }
    
    ### U: ground truth
    U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    
    ### ts:binary tensor
    ts = list()
    for (i in 1:dup) {
        binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
        ts[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data
    }
    
    return(list(U = U,ts = ts,A=A,B=B,C=C,G=G))
}


#####---- This is the function used for generating data through different distribution
#         of core tensor in supervised setting
#         dup: number of dupplicate generated from one ground truth
gene_data = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),p1,p2,dis,
                     gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1, dup){
  #dis can be "gaussian" and "uniform"

  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  ####-------- generate data
  set.seed(seed)  # 24 # 37  #  347

  X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 10),d1,p1)
  X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd = 10),d2,p2)
  W1 = randortho(p1,type = c("orthonormal"))[,1:r1]
  W2 = randortho(p2,type = c("orthonormal"))[,1:r2]

  C = randortho(d3,type = c("orthonormal"))[,1:r3]

  ### G: core tensor
  if(dis == "gaussian"){
    G = as.tensor(array(data = rnorm(r1*r2*r3,mean = gs_mean,sd = gs_sd),dim = core_shape))
  }
  else if(dis == "uniform"){
    G = as.tensor(array(data = runif(r1*r2*r3,min = unf_a,max = unf_b),dim = core_shape))
  }

  ### true coefficient
  C_ts = ttl(G,list(W1,W2,C),ms = c(1,2,3))@data
  ### U: ground truth
  U = ttl(G,list(X_covar1%*%W1,X_covar2%*%W2,C),ms = c(1,2,3))@data

  ### tsr:binary tensor
  tsr = list()
  for (i in 1:dup) {
    binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
    tsr[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data
  }

  return(list(X_covar1 = X_covar1, X_covar2 = X_covar2,C_ts = C_ts, U = U,tsr = tsr))
}
