rm(list=ls())
require("tensorsparse")

########################################
###   Table 3    #######################
########################################

n=40;p=40;q=40;k=3;r=5;l=4;iteration=50
out = list()
for (i in 1:3){
  result = c()
  for (iter in 1: iteration){
    set.seed(iter)
    error = c(1,4,8)[i]
    data = get.data(n,p,q,k,r,l,error)
    krl = choosekrl_bic(data$x,2:6,2:6,2:6)$estimated_krl
    ev = sparse.evaluate(classify2(data$x,krl[1],krl[2],krl[3]),data,show=FALSE)
    result = c(result,ev$cerC,ev$cerD,ev$cerE)
  }
  result = matrix(result,ncol=3,byrow=TRUE)
  meanre = apply(result, 2, mean)
  sdre = apply(result, 2, sd)
  cat("the cerC is", meanre[1], "(", sdre[1], "), the cerD is", meanre[2], "(", sdre[2], "), the cerE is", meanre[3], "(", sdre[3], ").\n")
}

n=50;p=50;q=50;k=3;r=5;l=4;iteration=50
out = list()
for (i in 1:3){
  result = c()
  for (iter in 1: iteration){
    set.seed(iter)
    error = c(1,4,8)[i]
    data = get.data(n,p,q,k,r,l,error)
    krl = choosekrl_bic(data$x,2:6,2:6,2:6)$estimated_krl
    ev = sparse.evaluate(classify2(data$x,krl[1],krl[2],krl[3]),data,show=FALSE)
    result = c(result,ev$cerC,ev$cerD,ev$cerE)
  }
  result = matrix(result,ncol=3,byrow=TRUE)
  meanre = apply(result, 2, mean)
  sdre = apply(result, 2, sd)
  cat("the cerC is", meanre[1], "(", sdre[1], "), the cerD is", meanre[2], "(", sdre[2], "), the cerE is", meanre[3], "(", sdre[3], ").\n")
}
