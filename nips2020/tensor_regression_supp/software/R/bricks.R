#library(rTensor)
library(pracma)
#library(gtools)
library(MASS)
library(speedglm)
#library(fastglm)
#library(rgl)
#library("RColorBrewer")



#######################
######----------   functions for update

## This part used for penalty (Conjugate gradient method)

loss = function(beta,y,X,lambda,alpha,dist){
  U = X %*% beta
  L= - loglike(y,U,dist)
  if(max(abs(U))>alpha) return(Inf)
  else{
    L2 =L - lambda * sum(log(1 - (U/ alpha)^2))  ## object function
    return(c(L2))
  }
}

loglike=function(data,linearpre,dist){
  if(dist=="binary"){
    p=plogis(linearpre)## log-likelihood
    L=dbinom(data,1,p,log=TRUE)
    return(sum(L))
  }
  else if (dist=="normal"){
    sigma_est=mean((data-linearpre)^2)
    L=dnorm(data,linearpre,sqrt(sigma_est),log=TRUE)
    sum(L)##log-likelihood
  }
  else if (dist=="poisson"){
    lambda=exp(linearpre)
    L=dpois(data,lambda, log = TRUE) ## log-likelihood for Poisson
    return(sum(L))
  }
}

loss_gr <- function(beta,y,X,lambda,alpha,dist){
  U = X %*% beta
  if(dist=="binary") p=plogis(U)
  else if (dist=="normal") p=U
  else if (dist=="poisson") p=exp(U)
  L_g=t(X) %*% (p - y)
  penal_g = 2 * lambda * t(X) %*% (U/(alpha^2 - U^2))
  return(c(L_g) + c(penal_g))
}


glm_modify=function(y,x,dist){
  if(dist=="binary"){
    fit1 =suppressWarnings(speedglm(y~-1+x,family=binomial(link="logit")))

    return(list(coef(fit1), logLik(fit1)))
  }
  else if (dist=="normal"){
    fit1 =speedlm(y~-1+x)

    return(list(coef(fit1), logLik(fit1)))
  }
  else if (dist=="poisson"){
    fit1 =suppressWarnings(speedglm(y~-1+x,family=poisson(link="log")))
    return(list(coef(fit1), logLik(fit1)))
  }
}

############################################################

#    form U = X*Beta ## in parallel
glm_mat = function(Y,X,dist){
  R = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  ma = mapply(glm_modify, y =  as.data.frame(Y),MoreArgs = list(X),dist=dist)
  re = t(matrix(unlist(ma), nrow = p, byrow=T))
  beta = re[1:R,]

  lglk = sum(re[R + 1,])
  return(list(t(beta),lglk))
}

###########---------  GLM on two modes
glm_two = function(Y, X1, X2, dist){
  ## Y_size = m * n
  # logit(E(Y)) = X1 %*% coe %*% X2
  m = dim(Y)[1] ; n = dim(Y)[2]
  q1 = dim(X1)[2] ; q2 = dim(X2)[1]

  N_long = kronecker_list(list(t(X2),X1))

  mod_re=glm_modify(as.vector(Y),N_long,dist)
  coe = mod_re[[1]]
  coe = matrix(coe, nrow = q1, ncol = q2)
  lglk= mod_re[[2]]
  return(list(coe=coe,lglk=lglk))
}



update_core=function(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist){

  M_long = kronecker_list(list(C,B,A))
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  violate=0

  if (cons == 'penalty'){
    if(max(abs(U))>=alpha){
      G=G/max(abs(U))*(alpha-0.01)
      U = U/max(abs(U))*(alpha-0.01)
      message("Violate constrain ------------------")
      violate = 1
    }
    mod_re = optim(par = as.vector(G@data),loss,loss_gr,y = as.vector(tsr@data),X =M_long, lambda = lambda, alpha = alpha, method = solver,dist=dist)
    coe = mod_re$par
    G = as.tensor(array(data = coe,dim = c(core_shape)))
    lglk = -mod_re$value
  }else {
    mod_re = glm_modify(as.vector(tsr@data), M_long,dist=dist)
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    lglk = loglike(tsr@data,U,dist)

    if(cons== 'vanilla' & max(abs(U))>=alpha){
      G=G/max(abs(U))*(alpha-0.01)
      U = U/max(abs(U))*(alpha-0.01)
      message("Violate constrain ------------------")
      violate = 1
    }
    lglk = loglike(tsr@data,U,dist)}


  return(list(G=G,violate=violate,lglk=lglk))
}



U_to_mean=function(U,dist){
  if(dist=="normal") return(U)
  else if (dist=="binary") return(plogis(U))
  else if (dist=="poisson") return(exp(U))
}
###
###----  This function shows how we select rank
##   recommend to use non-constrain verison to select rank



massive_glm=function(response,X,dist){
  d1=dim(response)[1]
  d2=dim(response)[2]
  coe=array(0,dim=c(d1,d2,dim(X)[2]))
  if(dist=="normal"){
    for(i in 1:d1){
      for(j in 1:d2){
        fit=lm(response[i,j,]~-1+X)
        coe[i,j,]=coef(fit)
      }
    }
  }
  else if(dist=="binary"){
    for(i in 1:d1){
      for(j in 1:d2){
        fit=suppressWarnings(glm(response[i,j,]~-1+X,family=binomial("logit")))
        coe[i,j,]=coef(fit)
      }
    }
  }
  else if(dist=="poisson"){
    for(i in 1:d1){
      for(j in 1:d2){
        fit=suppressWarnings(glm(response[i,j,]~-1+X,family=poisson("log")))
        coe[i,j,]=coef(fit)
      }
    }
  }
  return(coe)
}

