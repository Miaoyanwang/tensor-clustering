##########################################################
# R-code for Table 2
# K=2 R=4 
##########################################################
set.seed(1)
library(sparseBC)
Do <- function(n,p,k,r){
classification<-list()
for(a in 1:50){
cat("Starting", a, fill=TRUE)
  mus<-runif(k*r,-3,3)
  mus<-matrix(c(mus),nrow=k,ncol=r,byrow=FALSE)
  truthCs<-sample(1:k,n,rep=TRUE)
  truthDs<-sample(1:r,p,rep=TRUE)
  x<-matrix(rnorm(n*p,mean=0,sd=2),nrow=n,ncol=p)
  for(i in 1:max(truthCs)){
     for(j in 1:max(truthDs)){ 
         x[truthCs==i, truthDs==j] <- x[truthCs==i, truthDs==j] + mus[i,j]
     }
  }
  x<-x-mean(x)
  classification[[a]]<-sparseBC.choosekr(x,1:5,1:5,0,0.1)$estimated_kr
}
return(classification)
}

# Calculate accuracy of the estimated clusters
Calculate<-function(true,results){	
real<-matrix(true,ncol=2)
percent<-0
for(i in 1:length(results)){
  if(nrow(results[[i]])>1){
    for(a in 1:nrow(results[[i]])){
#  	  cat("iteration",a,fill=TRUE)
      if(sum(results[[i]][a,]==real)==2){
        percent<-percent+(1/nrow(results[[i]]))
      }
    }
  }
  else if(nrow(results[[i]])<2){
    if(sum(results[[i]]==real)==2){
       percent<-percent+1
      }	
  }
}
return(percent)
}

#Calculate Average k and r values
Calculatekr<-function(results){
#k<-0
#r<-0
k<-rep(NA,length(results))
r<-rep(NA,length(results))
for(i in 1:length(results)){
  if(nrow(results[[i]]>1)){
  	tempk<-0
  	tempr<-0
    for(a in 1:nrow(results[[i]])){
      tempk<-tempk+(results[[i]][a,1]/nrow(results[[i]]))
      tempr<-tempr+(results[[i]][a,2]/nrow(results[[i]]))
    }
    k[i]<-tempk
    r[i]<-tempr
  }
  else if(nrow(results[[i]]<2)){
    k[i]<-results[[i]][1,1]
    r[i]<-results[[i]][1,2]
  }
}
return(list(meank=mean(k),meanr=mean(r),sdek=sd(k)/sqrt(length(k)),sder=sd(r)/sqrt(length(r))))
#return(list(averagek=k/length(results),averager=r/(length(results))))
}

out <- Do(100,100,2,4)
res<-Calculate(c(2,4),out)
reskr<-Calculatekr(out)

out2 <- Do(100,500,2,4)
res2<-Calculate(c(2,4),out2)
reskr2<-Calculatekr(out2)

out3 <- Do(500,100,2,4)
res3<-Calculate(c(2,4),out3)
reskr3<-Calculatekr(out3)

out4 <- Do(500,500,2,4)
res4<-Calculate(c(2,4),out4)
reskr4<-Calculatekr(out4)

tab<-rbind(c("np","kr","accuracy","averagek","averager","sek","ser"),c("n=100,p=100","k=2 r=4",res*2,reskr[[1]],reskr[[2]],reskr[[3]],reskr[[4]]),c("n=100,p=500","k=2 r=4",res2*2,reskr2[[1]],reskr2[[2]],reskr2[[3]],reskr2[[4]]),c("n=500,p=100","k=2 r=4",res3*2,reskr3[[1]],reskr3[[2]],reskr3[[3]],reskr3[[4]]), c("n=500,p=500","k=2 r=4",res4*2,reskr4[[1]],reskr4[[2]],reskr4[[3]],reskr4[[4]]))
write.table(tab,file="Choosekrresults.txt",sep=",")
