rm(list=ls())
require("tensorsparse")

########################################
###   Table 4    #######################
########################################

set.seed(1)
n=40;p=40;q=40;k=4;r=4;l=4;error=4;sparserate=0.5;iteration = 50
lambda.sel = sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05)*error,error)
lambda_bar = mean(lambda.sel)
lambda_range=c(lambda_bar)
lambda_sd = sd(lambda.sel)
cat("The selected lambda is", lambda_bar, "(", lambda_sd, ").\n")
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
    result = sparse.evaluate(bires,raw,TRUE,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
    cerC[iter] = result$cerC
    cerD[iter] = result$cerD
    cerE[iter] = result$cerE
  }
  cat("when lambda =", lambda, ":\n")
  cat("estimated sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "sparsity error rate=", round(mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), "),\n",
      "cerC=", round(mean(cerC),4), "(", round(sd(cerC),4), "),\n",
      "cerD=", round(mean(cerD),4), "(", round(sd(cerD),4), "),\n",
      "cerE=", round(mean(cerC),4), "(", round(sd(cerE),4), ").\n"
  )
}

set.seed(2)
n=40;p=40;q=40;k=4;r=4;l=4;error=8;sparserate=0.5;iteration = 50
lambda.sel = sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05)*error,error)
lambda_bar = mean(lambda.sel)
lambda_range=c(lambda_bar)
lambda_sd = sd(lambda.sel)
cat("The selected lambda is", lambda_bar, "(", lambda_sd, ").\n")
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
    result = sparse.evaluate(bires,raw,TRUE,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
    cerC[iter] = result$cerC
    cerD[iter] = result$cerD
    cerE[iter] = result$cerE
  }
  cat("when lambda =", lambda, ":\n")
  cat("estimated sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "sparsity error rate rate=", round(mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), "),\n",
      "cerC=", round(mean(cerC),4), "(", round(sd(cerC),4), "),\n",
      "cerD=", round(mean(cerD),4), "(", round(sd(cerD),4), "),\n",
      "cerE=", round(mean(cerC),4), "(", round(sd(cerE),4), ").\n"
  )
}


set.seed(3)
n=40;p=40;q=40;k=4;r=4;l=4;error=8;sparserate=0.8;iteration = 50
lambda.sel = sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05)*error,error)
lambda_bar = mean(lambda.sel)
lambda_range=c(lambda_bar)
lambda_sd = sd(lambda.sel)
cat("The selected lambda is", lambda_bar, "(", lambda_sd, ").\n")
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
    result = sparse.evaluate(bires,raw,TRUE,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
    cerC[iter] = result$cerC
    cerD[iter] = result$cerD
    cerE[iter] = result$cerE
  }
  cat("when lambda =", lambda, ":\n")
  cat("estimated sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "sparsity error rate rate=", round(mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), "),\n",
      "cerC=", round(mean(cerC),4), "(", round(sd(cerC),4), "),\n",
      "cerD=", round(mean(cerD),4), "(", round(sd(cerD),4), "),\n",
      "cerE=", round(mean(cerC),4), "(", round(sd(cerE),4), ").\n"
  )
}
