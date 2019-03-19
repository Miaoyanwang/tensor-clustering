source("MVNLambdaKnown.R")
n=200
p=200
sd=2
iter=50
lambda=c(15,30,50,75,100,150,200,250,300,350,400,450)
matrixbicluster4<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=lambda,standarddeviation=sd,max.iter=50)

tab<-rbind(c("Method","rowCER","colCER","CorrectZero","CorrectNonZero","Sparsity","SparsityError"),c(paste("Matrix Known",matrixbicluster4$selectedlambda,sep=","),paste(round(matrixbicluster4$cerrow[1],digits=4),round(matrixbicluster4$cerrow[2],digits=4),sep=","),paste(round(matrixbicluster4$cercol[1],digits=4),round(matrixbicluster4$cercol[2],digits=4),sep=","),paste(round(matrixbicluster4$correctzero[1],digits=4),round(matrixbicluster4$correctzero[2],digits=4),sep=","),paste(round(matrixbicluster4$correctone[1],digits=4),round(matrixbicluster4$correctone[2],digits=4),sep=","),paste(round(matrixbicluster4$totalzero[1],digits=4),round(matrixbicluster4$totalzero[2],digits=4),sep=","),paste(round(matrixbicluster4$totalincorrect[1],digits=4),round(matrixbicluster4$totalincorrect[2],digits=4),sep=",")))


write(tab,file="n200p200Known.txt",sep=",")
