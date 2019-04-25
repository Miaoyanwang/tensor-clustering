rm(list=ls())
require("tensorsparse")

########################################
###   Table 4    #######################
########################################


set.seed(1)
cat("The first situation:\n")
n = 50; p = 50; q = 50; k = 3; r = 5; l = 4
our_mse = 1:3; cp_mse = 1:3
for (i in 1:3){
  cat("When the noise is", c(4,8,12)[i], ":\n")
  for (iter in 1:50){
    error = c(4,8,12)[i]
    try = get.data(n,p,q,k,r,l,error=error)
    krl = choosekrl_bic(try$x,2:5,2:5,2:5)
    k = krl$estimated_krl[1]; r = krl$estimated_krl[2]; l = krl$estimated_krl[3]
    lambda = chooseLambda(try$x,k,r,l,lambda=sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05))
    our_result = classify2(try$x,k,r,l,lambda$lambda)
    our_mse[i] = sparse.evaluate(our_result,try,CER=FALSE,show=FALSE)$mse
    cp_result = cp_kmeans(try$x,2:5,2:5,2:5)
    cp_mse[i] = sparse.evaluate(cp_result,try,CER=FALSE,show=FALSE)$mse
  }
  cat("The mse of our method is:", round(mean(our_mse),4),"(", round(sd(our_mse),4),");\n")
  cat("the mse of cp k-means is:", round(mean(cp_mse), 4),"(", round(sd(cp_mse),4), ").\n")
}


set.seed(2)
cat("The second situation:\n")
n = 50; p = 50; q = 50
our_mse = 1:3; cp_mse = 1:3
for (i in 1:3){
  cat("When the noise is", c(4,8,12)[i], ":\n")
  for (iter in 1:50){
    error = c(4,8,12)[i]
    try = get.data(n,p,q,multiplicative = 2, error=error)
    krl = choosekrl_bic(try$x,2:5,2:5,2:5)
    k = krl$estimated_krl[1]; r = krl$estimated_krl[2]; l = krl$estimated_krl[3]
    lambda = chooseLambda(try$x,k,r,l,lambda=sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05))
    our_result = classify2(try$x,k,r,l,lambda$lambda)
    our_mse[i] = sparse.evaluate(our_result,try,CER=FALSE,show=FALSE)$mse
    cp_result = cp_kmeans(try$x,2:5,2:5,2:5,)
    cp_mse[i] = sparse.evaluate(cp_result,try,CER=FALSE,show=FALSE)$mse
  }
  cat("The mse of our method is:", round(mean(our_mse),4),"(", round(sd(our_mse),4),");\n")
  cat("the mse of cp k-means is:", round(mean(cp_mse), 4),"(", round(sd(cp_mse),4), ").\n")
}


set.seed(3)
cat("The third situation:\n")
n = 50; p = 50; q = 50
for (i in 1:3){
  cat("When the noise is", c(4,8,12)[i], ":\n")
  for (iter in 1:50){
    our_mse = 1:50; cp_mse = 1:50
    error = c(4,8,12)[i]
    try = get.data(n,p,q,multiplicative = 1, error=error)
    krl = choosekrl_bic(try$x,2:5,2:5,2:5)
    k = krl$estimated_krl[1]; r = krl$estimated_krl[2]; l = krl$estimated_krl[3]
    lambda = chooseLambda(try$x,k,r,l,lambda=sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05))
    our_result = classify2(try$x,k,r,l,lambda$lambda)
    our_mse[i] = sparse.evaluate(our_result,try,CER=FALSE,show=FALSE)$mse
    cp_result = cp_kmeans(try$x,2:5,2:5,2:5)
    cp_mse[i] = sparse.evaluate(cp_result,try,CER=FALSE,show=FALSE)$mse
  }
  cat("The mse of our method is:", round(mean(our_mse),4),"(", round(sd(our_mse),4),");\n")
  cat("the mse of cp k-means is:", round(mean(cp_mse), 4),"(", round(sd(cp_mse),4), ").\n")
}







