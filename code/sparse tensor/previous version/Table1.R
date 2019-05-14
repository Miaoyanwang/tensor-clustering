rm(list=ls())
require(tensorsparse)

########################################
###   Table 1    #######################
########################################

set.seed(1)
n=40;p=40;q=40;k=3;r=5;l=4;iteration=50
out = list()
for (i in 1:3){
  error = c(4,8,12)[i]
  out[[i]] = simulation(n,p,q,k,r,l,error,iteration=iteration)
}

set.seed(2)
n=40;p=40;q=80;k=3;r=5;l=4;iteration=50
out1 = list()
for (i in 1:3){
  error = c(4,8,12)[i]
  out1[[i]] = simulation(n,p,q,k,r,l,error,iteration=iteration)
}

set.seed(3)
n=40;p=40;q=40;k=4;r=4;l=4;iteration=50
out2 = list()
for (i in 1:3){
  error = c(4,8,12)[i]
  out2[[i]] = simulation(n,p,q,k,r,l,error,iteration=iteration)
}

set.seed(4)
n=40;p=40;q=80;k=4;r=4;l=4;iteration=50
out3 = list()
for (i in 1:3){
  error = c(4,8,12)[i]
  out3[[i]] = simulation(n,p,q,k,r,l,error,iteration=iteration)
}
