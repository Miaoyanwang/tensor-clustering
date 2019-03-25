rm(list=ls())
require("sparseBCnew")
require("tensorsparse")
if(!require("clues")){
  install.packages("clues")
  stopifnot(require("clues"))
}


n = 200; p = 200; q = 1; k = 4; r = 5; l = 1
data = get.data(n,p,q,k,r,l,error=1,sort=TRUE,sparse.percent=0.7,mumin=-4,mumax=4)
test = data$x
truth = data$truthX
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs
mus = data$mus
binaryX = data$binaryX

lambda_0 =c(floor(n*p*q/k/r/l))
lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.1)


method = "L0"
verlam = chooseLambda(test,k,r,l,lambda=lambda,method=method)

plot(lambda,verlam$BIC,pch=16,col="pink")
cat("the lambda tensorsparse choose is", verlam$lambda,".\n")


splam = sparseBC.BIC(test[,,1],k,r,lambda = lambda,method=method)
points(lambda,splam$BIC,pch=17,col="lightblue")
cat("the lambda sparseBC choose is ", splam$lambda, ".\n")

#the result of tensor
sim = label2(test,4,5,1,threshold=1e-10,lambda=verlam$lambda,sim.times=5,trace=FALSE)
judgeX = sim$judgeX

plot_tensor(test)
Sys.sleep(1.5)

plot_tensor(truth)
Sys.sleep(1.5)

plot_tensor(judgeX)


cat("The result of with the lambda tensorsparse choose:\n")
sparse.evaluate(sim,data)

#result of sparseBC
sim = label2(test,4,5,1,threshold=1e-10,lambda=splam$lambda,sim.times=5,trace=FALSE)
judgeX = sim$judgeX

plot_tensor(test)
Sys.sleep(1.5)

plot_tensor(truth)
Sys.sleep(1.5)
plot_tensor(judgeX)


cat("The result of with the lambda sparseBC choose:\n")
sparse.evaluate(sim,data)

#compare the label result
tc = classify2(test,4,5,1,threshold=1e-10,lambda=3000,trace=FALSE)
mc = sparseBC(test[,,1],4,5,lambda=3000,center = FALSE)

similarityC<-adjustedRand(tc$Cs,mc$Cs,randMethod=c("Rand"))
similarityD<-adjustedRand(tc$Ds,mc$Ds,randMethod=c("Rand"))
cat("The similarity of two clustering result is ",similarityC,",",similarityD,".\n")

tcBIC = tensor.calculateBIC(test,tc,"L0")
mcBIC = CalculateBIC(test[,,1],mc,"L0")
cat("The BIC calculated by tensorsparse is", tcBIC[1], ", the BIC calculated by sparseBC is", mcBIC, ".\n")




