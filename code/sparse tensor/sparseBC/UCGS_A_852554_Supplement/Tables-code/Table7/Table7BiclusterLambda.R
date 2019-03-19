source("BiclusterLambda.R")
n=200
p=200
sd=2
iter=50
lambda=c(0,50,100,150,200,250,300,400,500,600,700,800,900,1000)
bicluster<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=lambda,standarddeviation=sd)

tab<-rbind(c("Method","rowCER","colCER","CorrectZero","CorrectNonZero","Sparsity","SparsityError"),c(paste("Biclustering auto",mean(bicluster$selectedlambda),sep=","),paste(round(bicluster$cerrow[1],digits=4),round(bicluster$cerrow[2],digits=4),sep=","),paste(round(bicluster$cercol[1],digits=4),round(bicluster$cercol[2],digits=4),sep=","),paste(round(bicluster$correctzero[1],digits=4),round(bicluster$correctzero[2],digits=4),sep=","),paste(round(bicluster$correctone[1],digits=4),round(bicluster$correctone[2],digits=4),sep=","),paste(round(bicluster$totalzero[1],digits=4),round(bicluster$totalzero[2],digits=4),sep=","),paste(round(bicluster$totalincorrect[1],digits=4),round(bicluster$totalincorrect[2],digits=4),sep=",")))


write(tab,file="n200p200bc.txt",sep=",")
