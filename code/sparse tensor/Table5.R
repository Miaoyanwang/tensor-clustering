rm(list=ls())
require("tensorsparse")

########################################
###   Table 5    #######################
########################################

set.seed(1)
n=40;p=40;q=40;k=3;r=5;l=4;error=4;sparserate=0.5;iteration = 50
lambda_bar = mean(sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05)*error,error))
lambda_range=c(0,100,200,lambda_bar)
sparsityrate = 1:iteration
correctzerorate = 1:iteration
correctonerate = 1:iteration
totalincorrectrate = 1:iteration
cerC = 1:iteration
cerD = 1:iteration
cerE = 1:iteration
for(lambda in lambda_range){
  for(iter in 1:iteration){
    set.seed(iter)
    raw = get.data(n,p,q,k,r,l,error,TRUE,sparserate)
    bires = classify2(raw$x,k,r,l,lambda)
    result = sparse.evaluate(bires,raw,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
    cerC[iter] = result$cerC
    cerD[iter] = result$cerD
    cerE[iter] = result$cerE
  }
  cat("when lambda =", lambda, ":\n")
  cat("sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "total correct rate=", round(1-mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), "),\n",
      "cerC=", round(mean(cerC),4), "(", round(sd(cerC),4), "),\n",
      "cerD=", round(mean(cerD),4), "(", round(sd(cerD),4), "),\n",
      "cerE=", round(mean(cerC),4), "(", round(sd(cerE),4), ").\n"
  )
}

set.seed(2)
n=40;p=40;q=40;k=3;r=5;l=4;error=8;sparserate=0.5;iteration = 50
lambda_bar = mean(sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05)*error,error))
lambda_range=c(0,100,200,lambda_bar)
sparsityrate = 1:iteration
correctzerorate = 1:iteration
correctonerate = 1:iteration
totalincorrectrate = 1:iteration
cerC = 1:iteration
cerD = 1:iteration
cerE = 1:iteration
for(lambda in lambda_range){
  for(iter in 1:iteration){
    set.seed(iter)
    raw = get.data(n,p,q,k,r,l,error,TRUE,sparserate)
    bires = classify2(raw$x,k,r,l,lambda)
    result = sparse.evaluate(bires,raw,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
    cerC[iter] = result$cerC
    cerD[iter] = result$cerD
    cerE[iter] = result$cerE
  }
  cat("when lambda =", lambda, ":\n")
  cat("sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "total correct rate=", round(1-mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), "),\n",
      "cerC=", round(mean(cerC),4), "(", round(sd(cerC),4), "),\n",
      "cerD=", round(mean(cerD),4), "(", round(sd(cerD),4), "),\n",
      "cerE=", round(mean(cerC),4), "(", round(sd(cerE),4), ").\n"
  )
}


set.seed(3)
n=40;p=40;q=40;k=3;r=5;l=4;error=8;sparserate=0.8;iteration = 50
lambda_bar = mean(sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt(2*(n*p*q)/(k*r*l))*seq(2,6,by=0.1),error))
lambda_range=c(0,100,200,lambda_bar)
sparsityrate = 1:iteration
correctzerorate = 1:iteration
correctonerate = 1:iteration
totalincorrectrate = 1:iteration
cerC = 1:iteration
cerD = 1:iteration
cerE = 1:iteration
for(lambda in lambda_range){
  #lambda = lambda_bar
  for(iter in 1:iteration){
    set.seed(iter)
    raw = get.data(n,p,q,k,r,l,error,TRUE,sparserate)
    bires = classify2(raw$x,k,r,l,lambda)
    result = sparse.evaluate(bires,raw,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
    cerC[iter] = result$cerC
    cerD[iter] = result$cerD
    cerE[iter] = result$cerE
  }
  cat("when lambda =", lambda, ":\n")
  cat("sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "total correct rate=", round(1-mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), "),\n",
      "cerC=", round(mean(cerC),4), "(", round(sd(cerC),4), "),\n",
      "cerD=", round(mean(cerD),4), "(", round(sd(cerD),4), "),\n",
      "cerE=", round(mean(cerC),4), "(", round(sd(cerE),4), ").\n"
  )
}
