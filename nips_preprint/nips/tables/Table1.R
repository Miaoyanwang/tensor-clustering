rm(list=ls())
require(tensorsparse)

########################################
###   Table 1    #######################
########################################

n=40;p=40;q=40;k=3;r=5;l=4;iteration=50
out = list()
for (i in 1:4){
  error = (i-1)*5
  out[[i]] = simulation(n,p,q,k,r,l,error,iteration=iteration)
}

n=50;p=50;q=50;k=3;r=5;l=4;iteration=50
out1 = list()
for (i in 1:4){
  error = (i-1)*5
  out1[[i]] = simulation(n,p,q,k,r,l,error,iteration=iteration)
}

