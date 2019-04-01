source("Table4Function.R")
n=200
p=200
sd=4
iter=50
kmeans<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=0,method="kmeans",rank=1,standarddeviation=sd)

bicluster1<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=0,method="biclustering",rank=1,standarddeviation=sd)

bicluster2<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=200,method="biclustering",rank=1,standarddeviation=sd)

bicluster3<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=500,method="biclustering",rank=1,standarddeviation=sd)

bicluster4<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1000,method="biclustering",rank=1,standarddeviation=sd)

bicluster5<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1000,method="biclusteringauto",rank=1,standarddeviation=sd)

plaid<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1200,method="plaid2",rank=1,standarddeviation=sd)

ssvd1<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1200,method="ssvd",rank=1,standarddeviation=sd)

ssvd2<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1200,method="ssvd",rank=2,standarddeviation=sd)

ssvd3<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1200,method="ssvd",rank=3,standarddeviation=sd)

ssvd4<-Do(n=n,p=p,k=4,r=5,iteration=iter,lambda=1200,method="ssvd",rank=4,standarddeviation=sd)

tab<-rbind(c("Method","rowCER","colCER","CorrectZero","CorrectNonZero","Sparsity","SparsityError"),c("Kmeans",paste(round(kmeans$cerrow[1],digits=4),round(kmeans$cerrow[2],digits=4),sep=","),paste(round(kmeans$cercol[1],digits=4),round(kmeans$cercol[2],digits=4),sep=","),"-","-","-","-"),c("Biclustering Lambda=0",paste(round(bicluster1$cerrow[1],digits=4),round(bicluster1$cerrow[2],digits=4),sep=","),paste(round(bicluster1$cercol[1],digits=4),round(bicluster1$cercol[2],digits=4),sep=","),"-","-","-","-"),c("Biclustering Lambda=200",paste(round(bicluster2$cerrow[1],digits=4),round(bicluster2$cerrow[2],digits=4),sep=","),paste(round(bicluster2$cercol[1],digits=4),round(bicluster2$cercol[2],digits=4),sep=","),paste(round(bicluster2$correctzero[1],digits=4),round(bicluster2$correctzero[2],digits=4),sep=","),paste(round(bicluster2$correctone[1],digits=4),round(bicluster2$correctone[2],digits=4),sep=","),paste(round(bicluster2$totalzero[1],digits=4),round(bicluster2$totalzero[2],digits=4),sep=","),paste(round(bicluster2$totalincorrect[1],digits=4),round(bicluster2$totalincorrect[2],digits=4),sep=",")),c("Biclustering Lambda=500",paste(round(bicluster3$cerrow[1],digits=4),round(bicluster3$cerrow[2],digits=4),sep=","),paste(round(bicluster3$cercol[1],digits=4),round(bicluster3$cercol[2],digits=4),sep=","),paste(round(bicluster3$correctzero[1],digits=4),round(bicluster3$correctzero[2],digits=4),sep=","),paste(round(bicluster3$correctone[1],digits=4),round(bicluster3$correctone[2],digits=4),sep=","),paste(round(bicluster3$totalzero[1],digits=4),round(bicluster3$totalzero[2],digits=4),sep=","),paste(round(bicluster3$totalincorrect[1],digits=4),round(bicluster3$totalincorrect[2],digits=4),sep=",")),c("Biclustering Lambda=1000",paste(round(bicluster4$cerrow[1],digits=4),round(bicluster4$cerrow[2],digits=4),sep=","),paste(round(bicluster4$cercol[1],digits=4),round(bicluster4$cercol[2],digits=4),sep=","),paste(round(bicluster4$correctzero[1],digits=4),round(bicluster4$correctzero[2],digits=4),sep=","),paste(round(bicluster4$correctone[1],digits=4),round(bicluster4$correctone[2],digits=4),sep=","),paste(round(bicluster4$totalzero[1],digits=4),round(bicluster4$totalzero[2],digits=4),sep=","),paste(round(bicluster4$totalincorrect[1],digits=4),round(bicluster4$totalincorrect[2],digits=4),sep=",")),c(paste("Biclustering Auto Lambda",round(mean(bicluster5$selectedlambda),digits=4),sep=","),paste(round(bicluster5$cerrow[1],digits=4),round(bicluster5$cerrow[2],digits=4),sep=","),paste(round(bicluster5$cercol[1],digits=4),round(bicluster5$cercol[2],digits=4),sep=","),paste(round(bicluster5$correctzero[1],digits=4),round(bicluster5$correctzero[2],digits=4),sep=","),paste(round(bicluster5$correctone[1],digits=4),round(bicluster5$correctone[2],digits=4),sep=","),paste(round(bicluster5$totalzero[1],digits=4),round(bicluster5$totalzero[2],digits=4),sep=","),paste(round(bicluster5$totalincorrect[1],digits=4),round(bicluster5$totalincorrect[2],digits=4),sep=",")),c("Plaid","-","-",paste(round(plaid$correctzero[1],digits=4),round(plaid$correctzero[2],digits=4),sep=","),paste(round(plaid$correctone[1],digits=4),round(plaid$correctone[2],digits=4),sep=","),paste(round(plaid$totalzero[1],digits=4),round(plaid$totalzero[2],digits=4),sep=","),paste(round(plaid$totalincorrect[1],digits=4),round(plaid$totalincorrect[2],digits=4),sep=",")),c("ssvd1","-","-",paste(round(ssvd1$correctzero[1],digits=4),round(ssvd1$correctzero[2],digits=4),sep=","),paste(round(ssvd1$correctone[1],digits=4),round(ssvd1$correctone[2],digits=4),sep=","),paste(round(ssvd1$totalzero[1],digits=4),round(ssvd1$totalzero[2],digits=4),sep=","),paste(round(ssvd1$totalincorrect[1],digits=4),round(ssvd1$totalincorrect[2],digits=4),sep=",")),c("ssvd2","-","-",paste(round(ssvd2$correctzero[1],digits=4),round(ssvd2$correctzero[2],digits=4),sep=","),paste(round(ssvd2$correctone[1],digits=4),round(ssvd2$correctone[2],digits=4),sep=","),paste(round(ssvd2$totalzero[1],digits=4),round(ssvd2$totalzero[2],digits=4),sep=","),paste(round(ssvd2$totalincorrect[1],digits=4),round(ssvd2$totalincorrect[2],digits=4),sep=",")),c("ssvd3","-","-",paste(round(ssvd3$correctzero[1],digits=4),round(ssvd3$correctzero[2],digits=4),sep=","),paste(round(ssvd3$correctone[1],digits=4),round(ssvd3$correctone[2],digits=4),sep=","),paste(round(ssvd3$totalzero[1],digits=4),round(ssvd3$totalzero[2],digits=4),sep=","),paste(round(ssvd3$totalincorrect[1],digits=4),round(ssvd3$totalincorrect[2],digits=4),sep=",")),c("ssvd4","-","-",paste(round(ssvd4$correctzero[1],digits=4),round(ssvd4$correctzero[2],digits=4),sep=","),paste(round(ssvd4$correctone[1],digits=4),round(ssvd4$correctone[2],digits=4),sep=","),paste(round(ssvd4$totalzero[1],digits=4),round(ssvd4$totalzero[2],digits=4),sep=","),paste(round(ssvd4$totalincorrect[1],digits=4),round(ssvd4$totalincorrect[2],digits=4),sep=",")))


















write(tab,file="Sparsen200p200sd4.txt",sep=",")
