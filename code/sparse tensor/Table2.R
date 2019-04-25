rm(list=ls())
require(tensorsparse)

########################################
###   Table 2    #######################
########################################

set.seed(1)
n=40;p=40;q=40;k=3;r=5;l=4;iteration=50
for (i in 1:3){
  error = c(4,8,12)[i]
  out = sim.choosekrl(n,p,q,k,r,l,error,iteration)
  res<-Calculate(c(k,r,l),out)
  reskr<-Calculatekrl(out)
  cat("The krl choosed by bic:\n")
  cat("The correct rate is", res,", meank=",reskr$meank, "(",reskr$sdek,"), meanr=",reskr$meanr, "(", reskr$sder, "), meanl=", reskr$meanl, "(", reskr$sdel, ").\n")
}

set.seed(2)
n=40;p=40;q=40;k=3;r=4;l=2;iteration=50
for (i in 1:3){
  error = c(4,8,12)[i]
  out = sim.choosekrl(n,p,q,k,r,l,error,iteration)
  res<-Calculate(c(k,r,l),out)
  reskr<-Calculatekrl(out)
  cat("The krl choosed by bic:\n")
  cat("The correct rate is", res,", meank=",reskr$meank, "(",reskr$sdek,"), meanr=",reskr$meanr, "(", reskr$sder, "), meanl=", reskr$meanl, "(", reskr$sdel, ").\n")
}

set.seed(3)
n=40;p=45;q=50;k=3;r=4;l=2;iteration=50
for (i in 1:3){
  error = c(4,8,12)[i]
  out = sim.choosekrl(n,p,q,k,r,l,error,iteration)
  res<-Calculate(c(k,r,l),out)
  reskr<-Calculatekrl(out)
  cat("The krl choosed by bic:\n")
  cat("The correct rate is", res,", meank=",reskr$meank, "(",reskr$sdek,"), meanr=",reskr$meanr, "(", reskr$sder, "), meanl=", reskr$meanl, "(", reskr$sdel, ").\n")
}
