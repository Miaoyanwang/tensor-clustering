source('tensorsparse.R')
source('plot.R')


#simulate the data matrix>>>>>>>>>>>>>>>>>>>>>>>>>
n=30;p=30;q=30;k=4;r=3;l=2
set.seed(123)
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
out <- sim.choosekrl(20,20,20,2,2,4,error=1)
res<-Calculate(c(2,2,4),out)
reskr<-Calculatekrl(out)

#test the chooseLambda()>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#need to redefine the cer!!!=====================!!!!!=================
n=20;p=30;q=40;k=2;r=3;l=4
set.seed(123)
data = get.data(n,p,q,k,r,l,error=2,sort=TRUE,sparse.percent=0.5)
test = data$x
truth = data$truthX
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs
mus = data$mus

lambda_summary = chooseLambda(test,k,r,l);print(lambda_summary)
lambda = lambda_summary$lambda
lambda_summary2 = chooseLambda2(test,k,r,l);print(lambda_summary2)
lambda2 = lambda_summary2$minimum


sim = label2(test,k,r,l,threshold=5e-2,lambda=lambda,sim.times=5,trace=FALSE)
judgeX = sim$judgeX
correct_zeros = sum(c(judgeX==0) & c(truth == 0))/sum(c(truth == 0))
cat("Correct zeros rate is", correct_zeros, ".\n")
wrong_zeros = sum(c(judgeX==0) & c(truth != 0))/sum(truth!=0)
cat("Wrong zeros rate is", wrong_zeros, ".\n")
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


sim = label2(test,k,r,l,threshold=5e-2,lambda=lambda2,sim.times=5,trace=FALSE)
judgeX = sim$judgeX
correct_zeros = sum(c(judgeX==0) & c(truth == 0))/sum(c(truth == 0))
cat("Correct zeros rate is", correct_zeros, ".\n")
wrong_zeros = sum(c(judgeX==0) & c(truth != 0))/sum(truth!=0)
cat("Wrong zeros rate is", wrong_zeros, ".\n")

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


lambda_summary3 = chooseLambda3(test,k,r,l);print(lambda_summary3)
lambda3 = lambda_summary3$minimum

sim = label3(test,k,r,l,threshold=5e-2,lambda=lambda3,sim.times=5,trace=FALSE)
judgeX = sim$judgeX
correct_zeros = sum(c(judgeX==0) & c(truth == 0))/sum(c(truth == 0))
cat("Correct zeros rate is", correct_zeros, ".\n")
wrong_zeros = sum(c(judgeX==0) & c(truth != 0))/sum(truth!=0)
cat("Wrong zeros rate is", wrong_zeros, ".\n")

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


lambda_summary4 = chooseLambda4(test,k,r,l);print(lambda_summary4)
lambda4 = lambda_summary4$minimum

sim = label4(test,k,r,l,threshold=5e-2,lambda=lambda4,sim.times=5,trace=FALSE)
judgeX = sim$judgeX
correct_zeros = sum(c(judgeX==0) & c(truth == 0))/sum(c(truth == 0))
cat("Correct zeros rate is", correct_zeros, ".\n")
wrong_zeros = sum(c(judgeX==0) & c(truth != 0))/sum(truth!=0)
cat("Wrong zeros rate is", wrong_zeros, ".\n")

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




lambda=0.1;n=30;p=30;q=30;k=5;r=3;l=6;sparse.percent=0.3
lambda_sim = sim.lambda_lasso(lambda,n,p,q,k,r,l,sparse.percent = sparse.percent)
print(lambda_sim)
