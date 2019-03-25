rm(list=ls())
require("sparseBCnew")
require("tensorsparse")
if(!require("clues")){
  install.packages("clues")
  stopifnot(require("clues"))
}


n = 200; p = 200; q = 1; k = 4; r = 5; l = 1
data = get.data(n,p,q,k,r,l,error=1,sort=TRUE,sparse.percent=0.7,mumin=-5,mumax=5)
test = data$x
truth = data$truthX
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs
mus = data$mus
binaryX = data$binaryX

lambda_0 =c(floor(n*p*q/k/r/l))


method = "L1"
if (method == "L0") lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.05)
if (method == "L1") lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.05)

verlam = chooseLambda(test,k,r,l,lambda=lambda,method=method)
cat("the lambda tensorsparse choose is", verlam$lambda,".\n")
cat("the nonzeromus calculated by tensorsparse is ", verlam$nonzeromus, ".\n")
cat("=================================================================\n")

splam = sparseBC.BIC(test[,,1],k,r,lambda = lambda,method=method)
plot(lambda,verlam$BIC,pch=16,col="pink",ylim=c(min(c(verlam$BIC,splam$BIC)),max(c(verlam$BIC,splam$BIC))))
points(lambda,splam$BIC,pch=17,col="lightblue")
cat("the lambda sparseBC choose is ", splam$lambda, ".\n")
cat("the nonzeromus calculated by sparseBC is ", splam$nonzeromus, ".\n")
cat("=================================================================\n")
#the result of tensor
sim = label2(test,4,5,1,threshold=1e-10,lambda=verlam$lambda,sim.times=5,trace=FALSE)
judgeX = sim$judgeX

plot_tensor(test)
Sys.sleep(0.5)

plot_tensor(truth)
Sys.sleep(0.5)

plot_tensor(judgeX)


cat("The result of with the lambda tensorsparse choose:\n")
sparse.evaluate(sim,data)
cat("=================================================================\n")
#result of sparseBC
sim = label2(test,4,5,1,threshold=1e-10,lambda=splam$lambda,sim.times=5,trace=FALSE)
judgeX = sim$judgeX

plot_tensor(test)
Sys.sleep(0.5)

plot_tensor(truth)
Sys.sleep(0.5)
plot_tensor(judgeX)


cat("The result of with the lambda sparseBC choose:\n")
sparse.evaluate(sim,data)
cat("=================================================================\n")
#compare the label result
tc = classify2(test,4,5,1,threshold=1e-10,lambda=verlam$lambda,trace=FALSE,method=method)
mc = sparseBC(test[,,1],4,5,lambda=verlam$lambda,center = FALSE,method=method)

similarityC<-adjustedRand(tc$Cs,mc$Cs,randMethod=c("Rand"))
similarityD<-adjustedRand(tc$Ds,mc$Ds,randMethod=c("Rand"))
cat("The similarity of two clustering result is ",similarityC,",",similarityD,".\n")

tcBIC = tensor.calculateBIC(test,tc,method=method)
mcBIC = CalculateBIC(test[,,1],mc,method=method)
cat("The BIC calculated by tensorsparse is", tcBIC[1], ", the BIC calculated by sparseBC is", mcBIC, ".\n")
cat("=================================================================\n")
tc = classify2(test,4,5,1,threshold=1e-10,lambda=splam$lambda,trace=FALSE,method=method)
mc = sparseBC(test[,,1],4,5,lambda=splam$lambda,center = FALSE,method=method)

similarityC<-adjustedRand(tc$Cs,mc$Cs,randMethod=c("Rand"))
similarityD<-adjustedRand(tc$Ds,mc$Ds,randMethod=c("Rand"))
cat("The similarity of two clustering result is ",similarityC,",",similarityD,".\n")

tcBIC = tensor.calculateBIC(test,tc,method=method)
mcBIC = CalculateBIC(test[,,1],mc,method=method)
cat("The BIC calculated by tensorsparse is", tcBIC[1], ", the BIC calculated by sparseBC is", mcBIC, ".\n")
cat("=================================================================\n")
#=================================
lambda = 4000
tc = classify2(test,4,5,1,lambda=lambda,trace=FALSE,method=method)
mc = sparseBC(test[,,1],4,5,lambda=lambda,center = FALSE,method=method)

similarityC<-adjustedRand(tc$Cs,mc$Cs,randMethod=c("Rand"))
similarityD<-adjustedRand(tc$Ds,mc$Ds,randMethod=c("Rand"))
cat("When lambda =", lambda, ":\n")
cat("The similarity of two clustering result is ",similarityC,",",similarityD,".\n")

tcBIC = tensor.calculateBIC(test,tc,method=method)
mcBIC = CalculateBIC(test[,,1],mc,method=method)
cat("The BIC calculated by tensorsparse is", tcBIC[1], ", the BIC calculated by sparseBC is", mcBIC, ".\n")
cat("=================================================================\n")


