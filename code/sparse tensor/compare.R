source('tensorsparse.R')
source('plot.R')
if(!require("sparseBC")){
  install.packages("sparseBC")
  stopifnot(require("sparseBC"))
}



n=30;p=30;q=1;k=3;r=3;l=1
data = get.data(n,p,q,k,r,l,error=2,sort=TRUE)
test = data$x
truth = data$truthX
plot_tensor(test)
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs

#better method to get labels>>>>>>>>>>>>>>>>>>>>>>>
sim = label2(test,k,r,l,threshold=5e-2,lambda=0,sim.times=5,trace=FALSE)
judgeX = sim$judgeX
#true distribution of mu
plot_tensor(truth)
Sys.sleep(1.5)
#plot_tensor_2(reorderClusters(truth,truthCs,truthDs,truthEs))
#the input data matrix
plot_tensor(test)
Sys.sleep(1.5)
cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
cat("The CER(clustering error rate) is ",cerC,",",cerD,",",cerE,".\n")
#the result of classifying
plot_tensor(judgeX)



bicluster = sparseBC(matrix.x,k,r,0)

tencluster = sim.chooseLambda(200,200,1,4,5,1,iteration = 5,lambda=c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200),4)

verlam = chooseLambda(tensor.x,4,5,1,lambda=c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200))
plot(lambda,verlam$BIC)
plot(1:length(lambda),verlam$logRSS,pch=16,col="pink",ylim = c(-50,24))
points(1:length(lambda),verlam$nonzeromus,pch=17,col="lightblue")

splam = sparseBC.BIC(matrix.x,4,5,lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200))
plot(lambda,splam$BIC)

set.seed(11)
test = get.data(200,200,1,4,5,1,1,FALSE,0.5,TRUE)$truthX
verlam = chooseLambda(test,4,5,1,lambda=c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200))
plot(lambda,verlam$BIC,pch=16)

splam = sparseBC.BIC(test[,,1],4,5,lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200))
points(lambda,splam$BIC,pch=17)

tl = label2(test,4,5,1,100)
ml = sparseBC(test[,,1],4,5,100)
print(adjustedRand(ml$Cs,tl$Cs))
print(adjustedRand(ml$Ds,tl$Ds))
