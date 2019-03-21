source('tensorsparse.R')
source('plot.R')
require("sparseBCnew")

#simulate the data matrix>>>>>>>>>>>>>>>>>>>>>>>>>
n=30;p=30;q=30;k=3;r=3;l=3
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


#get the mean and sd of cer in 50 simulations>>>>>>>>>>>>>>>>>>>>>>>>
n=30;p=30;q=30;k=3;r=3;l=3;error=1;lambda=0;iteration=50
simulation(n,p,q,k,r,l,error,lambda,iteration)


#one method to selecting k,r,l>>>>>>>>>>>>>>>>>>>>>>
#find the bottleneck
n=30;p=30;q=30;k=3;r=3;l=3
data = get.data(n,p,q,k,r,l,error=2,sort=TRUE)
test = data$x
truth = data$truthX
plot_tensor(test)
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs
#source(file="simulation.R", keep.source=TRUE)
Rprof(filename="Rprof.out", line.profiling=TRUE)
range.k = 2:4; range.r = 2:4; range.l = 2:4
sparse.choosekrl(test,range.k,range.r,range.l,trace=TRUE)
Rprof(NULL)
print(summaryRprof(filename="Rprof.out", lines="show"))


#evaluate the sparse.choosekrl()>>>>>>>>>>>>>>>>>
out <- sim.choosekrl(20,20,20,2,2,4,error=1,sim.times=5)
res<-Calculate(c(2,2,4),out)
reskr<-Calculatekrl(out)
cat("The correct rate is", res,", the mean krl and sd of krl we choose is", reskr, ".\n")

#test the chooseLambda()>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
n=200;p=200;q=1;k=4;r=5;l=1
set.seed(1)
data = get.data(n,p,q,k,r,l,error=2,sort=TRUE,sparse.percent=0.5)
test = data$x
truth = data$truthX
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs
mus = data$mus
binaryX = data$binaryX
plot_tensor(binaryX)

lambda_summary = chooseLambda(test,k,r,l);print(lambda_summary)
print(lambda_summary)
lambda = lambda_summary$lambda

sim<-label2(test,k,r,l,lambda)

cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
cat("The CER(clustering error rate) is ",cerC,",",cerD,",",cerE,".\n")


plot_tensor(test)
Sys.sleep(1.5)
plot_tensor(truth)
Sys.sleep(1.5)
judgeX = sim$judgeX
plot_tensor(judgeX)

restensor<-(abs(judgeX)>1e-10)*1
totalzero<-sum(restensor==0)/(n*p*q)
correctzero<-sum(restensor[which(binaryX==0)]==0)/sum(binaryX==0)
correctone<-sum(restensor[which(binaryX==1)]==1)/(n*p-sum(binaryX==0))
totalincorrect<-sum(abs(restensor-binaryX))/(n*p)
cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect,".")

