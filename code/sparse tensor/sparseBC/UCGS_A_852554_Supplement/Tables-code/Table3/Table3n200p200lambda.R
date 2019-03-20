#########################
# Run AutoBIC
#########################
source("Table3Functionslambda.R")
n=200
p=200
sd=4
iter=50
lambda<-c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)

bicluster5<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=lambda,standarddeviation=sd)

tab<-rbind(c("Method","rowCER","colCER","CorrectZero","CorrectNonZero","Sparsity","SparsityError"),c(paste("Biclustering Auto Lambda",round(mean(bicluster5$selectedlambda),digits=4),sep=","),paste(round(bicluster5$cerrow[1],digits=4),round(bicluster5$cerrow[2],digits=4),sep=","),paste(round(bicluster5$cercol[1],digits=4),round(bicluster5$cercol[2],digits=4),sep=","),paste(round(bicluster5$correctzero[1],digits=4),round(bicluster5$correctzero[2],digits=4),sep=","),paste(round(bicluster5$correctone[1],digits=4),round(bicluster5$correctone[2],digits=4),sep=","),paste(round(bicluster5$totalzero[1],digits=4),round(bicluster5$totalzero[2],digits=4),sep=","),paste(round(bicluster5$totalincorrect[1],digits=4),round(bicluster5$totalincorrect[2],digits=4),sep=",")))

write(tab,file="autobicn200p200.txt",sep=",")
