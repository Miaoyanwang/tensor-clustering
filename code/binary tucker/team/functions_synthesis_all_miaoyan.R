library(rTensor)
library(pracma)
library(gtools)
library(MASS)
library(speedglm)
#library(fastglm)
#library(rgl)
library("RColorBrewer")

#######################
######----------   functions for update

## This part used for penalty (Conjugate gradient method)

KL=function(U_est,U_true,dist,sigma_est=1,sigma_true=1){
    d=prod(dim(U_est))
    if(dist=="normal"){
      KL=d*(0.5*log(sigma_true)-0.5*log(sigma_est)+sigma_est/(2*sigma_true)-0.5)+0.5/sigma_true*sum((U_est-U_true)^2)
    }else if (dist=="binary"){
        p_est=plogis(U_est)
        p_true=plogis(U_true)
        KL=sum(p_est*log(p_est/p_true)+(1-p_est)*log((1-p_est)/(1-p_true)))
    }else if (dist=="poisson"){
        lambda_est=exp(U_est)
        lambda_true=exp(U_true)
        KL=sum(lambda_est*log(lambda_est/lambda_true)+lambda_true-lambda_est)
    }
return(KL)
}

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

## This part follows Sec 2.1.2 (?) to compute Gradient of object function 
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
          fit1 =speedglm(y~-1+x,family=binomial(link="logit"))
          
        return(list(coef(fit1), logLik(fit1)))
      }
      else if (dist=="normal"){
          fit1 =speedlm(y~-1+x)
          
          return(list(coef(fit1), logLik(fit1)))
      }
      else if (dist=="poisson"){
          fit1 =speedglm(y~-1+x,family=poisson(link="log"))
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
##---  This section follows 1.3 of [Algorithm: Semi-Supervised Binary Tensor Factorization]
##---  function version
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
            print("Violate constrain ------------------")
            violate = 1
        }
        mod_re = optim(par = as.vector(G@data),loss,loss_gr,y = as.vector(tsr@data),X =M_long, lambda = lambda, alpha = alpha, method = solver,dist=dist)
        coe = mod_re$par
        G = as.tensor(array(data = coe,dim = core_shape))
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
            print("Violate constrain ------------------")
            violate = 1
        }
        lglk = loglike(tsr@data,U,dist)}
 
        
    return(list(G=G,violate=violate,lglk=lglk))
}

#################  main update function ###--------   
update_binary_all = function(tsr,X_covar1 = NULL, X_covar2 = NULL,X_covar3 = NULL, core_shape, Nsim=20, cons, lambda = 0.1, alpha = 1, solver = NULL,dist){
 
                                 
  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  
  ###  check whether unsupervised on each mode
  un_m1 = FALSE ; un_m2 = FALSE ; un_m3 = FALSE 
  if(is.null(X_covar1)|(identical(X_covar1,diag(d1)))) {X_covar1 = diag(d1) ; un_m1 = TRUE}
  if(is.null(X_covar2)|(identical(X_covar2,diag(d2)))) {X_covar2 = diag(d2) ; un_m2 = TRUE}
  if(is.null(X_covar3)|(identical(X_covar3,diag(d3)))) {X_covar3 = diag(d3) ; un_m3 = TRUE}
  p1 = dim(X_covar1)[2] ; p2 = dim(X_covar2)[2] ; p3 = dim(X_covar3)[2]
  
    
  if(dist=="binary"){
      tsr.transform=as.tensor(2*tsr@data-1)
  }else if(dist=="poisson"){
      tsr.transform=as.tensor(log(tsr@data+0.1))###?? new initilization
 }else if (dist=="normal"){
    tsr.transform=tsr
}

  C_ts=ttl(tsr.transform,list(ginv(X_covar1),ginv(X_covar2),ginv(X_covar3)),ms=c(1,2,3))
  
  tckr = tucker(C_ts, ranks = core_shape)
  W1 = tckr$U[[1]] ; W2 = tckr$U[[2]] ; W3 = tckr$U[[3]] ## tucker factors 
  G = tckr$Z
  A = X_covar1%*%W1
  B = X_covar2%*%W2
  C = X_covar3%*%W3
  
  core=update_core(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist)
  G=core$G
  lglk=core$lglk
  violate=core$violate
  
  for(n in 1:Nsim){
    ## parameter from previous step
    
    W10 = W1 ; W20 = W2 ; W30 = W3 ; G0=G; A0=A;B0=B;C0=C;lglk0=tail(lglk,1);

    ###### update W1
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    
    if(un_m1) {re = glm_mat(t(Y_1),t(G_BC1),dist=dist) ## no covariate
    } else {re = glm_two(Y = Y_1, X1 = X_covar1, X2 = G_BC1, dist=dist)}
    
    W1 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W1*
    qr_res=qr(W1)
    W1=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),1)
    print("W1 Done------------------")
    
    ##### calculate A
    A = X_covar1%*%W1;
    ##### update W2
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
    
    if(un_m2) {re = glm_mat(t(Y_2),t(G_AC2),dist=dist)
    } else {re = glm_two(Y_2, X_covar2, G_AC2, dist=dist)}
    
    W2 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W2* 
    qr_res=qr(W2)
    W2=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),2)
    print("W2 Done------------------")
    
    ##### calculate B
    B = X_covar2%*%W2;
    
    
    ###### update W3
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
    
    if(un_m3) {re = glm_mat(t(Y_3),t(G_AB3),dist=dist)
    } else {re = glm_two(Y_3, X_covar3, G_AB3,dist=dist)}

    W3 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W3*
    qr_res=qr(W3)
    W3=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),3)
    print("W3 Done------------------")
    
    ##### calculate C
    C = X_covar3%*%W3;
    
    #########-----------------------------------------------
    ###  obtain core tensor under constraint 
    core=update_core(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist)
    G=core$G
    lglk=c(lglk,core$lglk)
    violate=c(violate,core$violate)
    
    print("G Done------------------")
    
    print(paste(n,"-th  iteration -- when dimension is",d1,d2,d3,"- rank is ",r1,r2,r3," -----------------"))

if (tail(lglk,1)-lglk0 <= 0.0005 & tail(lglk,1)>= lglk0 ) break
   else if (tail(lglk,1)-lglk0 < 0) {
       W1 = W10 ; W2 = W20 ; W3 = W30; G=G0; lglk=lglk[-c((length(lglk)-3):length(lglk))]; 
       A=A0;B=B0;C=C0;
       break
    } 
    
  }
  
U=ttl(G,list(A,B,C),ms = c(1,2,3))@data

  sigma_est=mean((tsr@data-U_to_mean(U,dist))^2)

  return(list(W1 = W1,W2 = W2,W3 = W3,G = G,U=U, C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data,lglk = lglk, sigma=sigma_est,violate = violate))
}


U_to_mean=function(U,dist){
    if(dist=="normal") return(U)
    else if (dist=="binary") return(plogis(U))
    else if (dist=="poisson") return(exp(U))
}
###
###----  This function shows how we select rank 
##   recommend to use non-constrain verison to select rank


sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,
                     rank1 = c(3,5), rank2=c(3,5),rank3=c(4,5),Nsim,cons = 'non',dist){
                         #rank = 1:3
   whole_shape=dim(tsr)                      
   p=rep(0,3)
   if(is.null(X_covar1)) p[1]=whole_shape[1] else p[1]=dim(X_covar1)[2]
   if(is.null(X_covar2)) p[2]=whole_shape[2] else p[2]=dim(X_covar2)[2]
   if(is.null(X_covar3)) p[3]=whole_shape[3] else p[3]=dim(X_covar3)[2]
  
                         
  rank = expand.grid(rank1,rank2,rank3)
  rank=rank[which(apply(rank,1,function(x){max(x)<= sqrt(prod(x))})==1),]
  rank=rank[which(apply(rank,1,function(x){sum(x<=p)==3})),]
  rank_matrix=rank
  rank=as.matrix(rank)
  
  whole_shape = dim(tsr)
  rank = lapply(1:dim(rank)[1], function(x) rank[x,]) ## turn rank to a list
  upp = lapply(rank, FUN= update_binary_all,tsr = tsr,X_covar1 = X_covar1,X_covar2 = X_covar2,X_covar3 = X_covar3, Nsim = Nsim, cons = cons,dist=dist)
 
  lglk= unlist(lapply(seq(length(upp)), function(x) tail(upp[[x]]$lglk,1)))
  BIC = unlist(lapply(seq(length(rank)), function(x) (prod(rank[[x]]) + sum((p-1)*rank[[x]])) * log(prod(whole_shape))))
  BIC = -2*lglk + BIC
  rank_matrix=cbind(rank_matrix,lglk,BIC)
  
  return(list(rank = rank[[which(BIC == min(BIC))]],result=rank_matrix))
}



#### This function is to select the lambda in the penalty constrain version
sele_lambda = function(seed, lambda, ...){
  #lambda = as.list(lambda)
  re = lapply(lambda, FUN = conv_rate, seed = seed, ...)
  re = lapply(seq(length(re)), function(x) re[[x]]$RMSE)
  return(re)
}

# re = sele_lambda(seed = 24, lambda = c(0.1,1,10,100,1000),  d = 20,r = 3, p1 = 5, p2 = 5, p3 = 5, dis = 'gaussian',gs_mean = 0,gs_sd = 1,unf_a = 0,unf_b = 1, 
#                  dup = 1, signal = 10, Nsim = 50, linear = TRUE, cons = 'penalty' ,
#                  solver = 'CG')



###----  functions for simulation
#####---- This is the function used for generating data through different distribution
#         of core tensor in  semi-supervised setting
## p is the dimension of the covaraite. p = 0 represents the case without covaraites
gene_data_all = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),p=c(3,3,0),dist, dup, signal){ 
    
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  p1 = p[1]; p2 = p[2]; p3 = p[3];
  ####-------- generate data
  set.seed(seed)  # 24 # 37  #  347
  X_covar1 = X_covar2 = X_covar3 = NULL  

if(p1<=0){
    X_covar1=diag(1,d1)
    p1=d1
}else{
    X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 1/sqrt(d1)),d1,p1)
}
W1 =as.matrix(randortho(p1)[,1:r1])   
A = X_covar1%*%W1 ## factor matrix
 
if(p2<=0){
     X_covar2=diag(1,d2)
     p2=d2
 }else{
     X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd =1/sqrt(d2)),d2,p2)
 }
 W2 = as.matrix(randortho(p2)[,1:r2])   
 B = X_covar2%*%W2 ## factor matrix
 
if(p3<=0){
     X_covar3=diag(1,d3)
     p3=d3
}else{
     X_covar3 = matrix(rnorm(d3*p3,mean = 0, sd = 1/sqrt(d3)),d3,p3)
}
W3 = as.matrix(randortho(p3)[,1:r3])   
C= X_covar3%*%W3 ## factor matrix 
 
  ### G: core tensor
  G = as.tensor(array(data = rnorm(r1*r2*r3,mean = 0,sd = 1),dim = core_shape))
 
  
  ### U: linear predictor
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  G=G/max(abs(U))*signal ## rescale subject to entrywise constraint
  U=U/max(abs(U))*signal
  
  C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data ## coefficient
  
  ### tsr:binary tensor
  if(dist=="binary"){
      tsr = lapply(seq(dup), function(x) array(rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)))),dim = c(d1,d2,d3)))}
  else if (dist=="normal"){
       tsr = lapply(seq(dup), function(x) array(rnorm(d1*d2*d3,U),dim = c(d1,d2,d3)))}
  else if (dist=="poisson"){
      tsr = lapply(seq(dup), function(x) array(rpois(d1*d2*d3,exp(U)),dim = c(d1,d2,d3)))}

  
  return(list(tsr = tsr,X_covar1 = X_covar1, X_covar2 = X_covar2,X_covar3 = X_covar3,
              W1 = W1,W2 = W2,W3 = W3, G=G, U=U,C_ts=C_ts))
}





conv_rate = function(seed,signal=10,Nsim=20,cons="no",lambda = 1,alpha=10,solver ="GC",c_range,dist="binary",dup=10,d_range,p_range,match_dp=TRUE){
  #cons can be "non","vanilla","penalty"
  if(match_dp==TRUE){
      n=m=nrow(d_range)
      s=nrow(c_range)
      error=array(0,dim=c(n,s,dup))
      error_KL=array(0,dim=c(n,s,dup))
      for(i in 1:n){
              for(k in 1:s){
                  data = gene_data_all(seed, d_range[i,], c_range[k,],p_range[i,],dist, dup, signal)
                  X_covar1=data$X_covar1
                  X_covar2=data$X_covar2
                  X_covar3=data$X_covar3
                  for(l in 1:dup){
result=update_binary_all(data$tsr[[l]],X_covar1,X_covar2,X_covar3,c_range[k,],Nsim,cons,lambda,alpha,solver,dist)
                      error[i,k,l]=mean((result$U-data$U)^2)
                      error_KL[i,k,l]=KL(result$U,data$U,dist,sigma_est=result$sigma)
                      print(paste(l,"-th replicate---- when dimension is ",d_range[i,],"-- covariate is ",p_range[i,]," --------core is ",c_range[k,]))
                      
                      
  }
              }
      }
      return(list(error,error_KL))  
  }
  else if (match_dp==FALSE){    
  n=nrow(d_range)
  m=nrow(p_range)
  s=nrow(c_range)
  error=array(0,dim=c(n,m,s,dup))
  error_KL=array(0,dim=c(n,m,s,dup))
  for(i in 1:n){
      for(j in 1:m){
          for(k in 1:s){
          data = gene_data_all(seed, d_range[i,], c_range[k,],p_range[j,],dist, dup, signal)
          X_covar1=data$X_covar1
          X_covar2=data$X_covar2
          X_covar3=data$X_covar3
          for(l in 1:dup){
result=update_binary_all(data$tsr[[l]],X_covar1,X_covar2,X_covar3,c_range[k,],Nsim,cons,lambda,alpha,solver,dist)
error[i,j,k,l]=mean((result$U-data$U)^2)
error_KL[i,j,k,l]=KL(result$U,data$U,dist,sigma_est=result$sigma)
print(paste(l,"-th replicate---- when dimension is ",d_range[i,],"-- covariate is ",p_range[j,]," --------core is ",c_range[k,]))
          }
          }
      }
  }      
 return(list(error,error_KL))  
  }
}
 
