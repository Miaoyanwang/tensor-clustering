rm(list=ls())
require("tensorsparse")

########################################
###   Table 5    #######################
########################################

n=50;p=50;q=50;k=3;r=5;l=4;error=1;sparserate=0.5;rangekrl=3:6
iteration = 50
lambda_bar = mean(sim.chooseLambda(n,p,q,k,r,l,sparserate,iteration,sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05),error))
lambda_range=c(0,5,10,15,lambda_bar)
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
    selectedkrl = choosekrl_bic(raw$x,rangekrl,rangekrl,rangekrl)$estimated_krl
    bires = classify2(raw$x,selectedkrl[1,1],selectedkrl[1,2],selectedkrl[1,3],lambda)
    result = sparse.evaluate(bires,raw,FALSE)
    sparsityrate[iter] = result$sparsityrate
    correctzerorate[iter] = result$correctzerorate
    correctonerate[iter] = result$correctonerate
    totalincorrectrate[iter] = result$totalincorrectrate
  }
  cat("when lambda =", lambda, ":\n")
  cat("sparsity rate=", round(mean(sparsityrate),4), "(", round(sd(sparsityrate),4), "),\n",
      "correct zero rate=", round(mean(correctzerorate),4), "(", round(sd(correctzerorate),4), "),\n",
      "correct one rate=", round(mean(correctonerate),4), "(", round(sd(correctonerate),4), "),\n",
      "total incorrect rate=", round(mean(totalincorrectrate),4), "(", round(sd(totalincorrectrate),4), ").\n",
  )
}
