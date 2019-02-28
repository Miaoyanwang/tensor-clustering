source('tensorsparse.R')

#simulate the data matrix>>>>>>>>>>>>>>>>>>>>>>>>>
n=30;p=30;q=30;k=4;r=3;l=2
data = get.data(n,p,q,k,r,l,error=2,sort=TRUE)
test = data$x
truth = data$truthX
source('plot.R')
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
#plot_tensor(reorderClusters(truth,truthCs,truthDs,truthEs))
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
n=50;p=50;q=50;k=3;r=3;l=3;error=3;lambda=0;iteration=50
simulation(n,p,q,k,r,l,error,lambda,iteration)


#one method to selecting k,r,l>>>>>>>>>>>>>>>>>>>>>>
range.k = 2:4; range.r = 2:4; range.l = 2:4
sparse.choosekrl(test,range.k,range.r,range.l,trace=TRUE)

#evaluate the sparse.choosekrl()>>>>>>>>>>>>>>>>>
out <- sim.choosekrl(20,40,30,2,4,3)
res<-Calculate(c(2,4,3),out)
reskr<-Calculatekrl(out)




