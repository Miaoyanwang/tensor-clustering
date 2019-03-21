source('tensorsparse.R')
source('plot.R')
require("sparseBCnew")

CalculateBIC <- function(x,bires,method="new"){
  if(method=="new"){
    RSS=log(sum((x-bires$mus)^2))*nrow(x)*ncol(x)
    df=log(nrow(x)*ncol(x))*sum(bires$Mus!=0)
    return(RSS+df)
  } else {
    mat <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    Cs <- bires$Cs
    Ds <- bires$Ds
    for(i in unique(Cs)){
      for(j in unique(Ds)){
        if(bires$Mus[i,j]!=0){
          mat[Cs==i,Ds==j] <- mean(x[Cs==i,Ds==j])
        }
      }
    }
    mat[abs(bires$mus)<=1e-8] <- mean(x[abs(bires$mus)<=1e-8])
    return(log(sum((x-mat)^2))*nrow(x)*ncol(x) + log(nrow(x)*ncol(x))*sum(bires$Mus!=0))
  }
}



set.seed(11)
n = 200; p = 200; q = 1; k = 4; r = 5; l = 1
data = get.data(n,p,q,k,r,l,error=1,sort=TRUE,sparse.percent=0.5,mumin=-4,mumax=4)
test = data$x
truth = data$truthX
truthCs = data$truthCs
truthDs = data$truthDs
truthEs = data$truthEs
mus = data$mus
binaryX = data$binaryX

lambda_0 =c(floor(n*p*q/k/r/l))
lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.1)


method = "new"
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

restensor<-(abs(judgeX)>1e-10)*1
totalzero<-sum(restensor==0)/(n*p*q)
correctzero<-sum(restensor[which(binaryX==0)]==0)/sum(binaryX==0)
correctone<-sum(restensor[which(binaryX==1)]==1)/(n*p*q-sum(binaryX==0))
totalincorrect<-sum(abs(restensor-binaryX))/(n*p*q)
cat("The result of with the lambda tensorsparse choose:\n")
cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect,".\n")

#result of sparseBC
sim = label2(test,4,5,1,threshold=1e-10,lambda=splam$lambda,sim.times=5,trace=FALSE)
judgeX = sim$judgeX

plot_tensor(test)
Sys.sleep(1.5)

plot_tensor(truth)
Sys.sleep(1.5)
plot_tensor(judgeX)

restensor<-(abs(judgeX)>1e-10)*1
totalzero<-sum(restensor==0)/(n*p)
correctzero<-sum(restensor[which(binaryX==0)]==0)/sum(binaryX==0)
correctone<-sum(restensor[which(binaryX==1)]==1)/(n*p-sum(binaryX==0))
totalincorrect<-sum(abs(restensor-binaryX))/(n*p)
cat("The result of with the lambda sparseBC choose:\n")
cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect,".\n")

#compare the label result
tc = classify2(test,4,5,1,threshold=1e-10,lambda=3000,trace=FALSE)
mc = sparseBC(test[,,1],4,5,lambda=3000,center = FALSE)

similarityC<-adjustedRand(tc$Cs,mc$Cs,randMethod=c("Rand"))
similarityD<-adjustedRand(tc$Ds,mc$Ds,randMethod=c("Rand"))
cat("The similarity of two clustering result is ",similarityC,",",similarityD,".\n")

tcBIC = tensor.calculateBIC(test,tc,"original")
mcBIC = CalculateBIC(test[,,1],mc,"original")
cat("The BIC calculated by tensorsparse is", tcBIC[1], ", the BIC calculated by sparseBC is", mcBIC, ".\n")




