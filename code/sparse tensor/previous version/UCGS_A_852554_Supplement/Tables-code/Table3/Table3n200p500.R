source("Table3Functions.R")
out<-Do(200,500,4,5,200,400,800,50)
tab<-rbind(c("Method","rowCER","colCER","Sparsity"),c("Kmeans",paste(out[[1]][[1]],out[[1]][[2]],sep=","),paste(out[[6]][[1]],out[[6]][[2]],sep=","),0), c("Bicluster0",paste(out[[2]][[1]],out[[2]][[2]],sep=","),paste(out[[7]][[1]],out[[7]][[2]],sep=","),0), c("BiCluster1",paste(out[[3]][[1]],out[[3]][[2]],sep=","),paste(out[[8]][[1]],out[[8]][[2]],sep=","),paste(out[[11]][[1]],out[[11]][[2]],sep=",")),c("BiCluster2",paste(out[[4]][[1]],out[[4]][[2]],sep=","),paste(out[[9]][[1]],out[[9]][[2]],sep=","),paste(out[[12]][[1]],out[[12]][[2]],sep=",")),c("BiCluster3",paste(out[[5]][[1]],out[[5]][[2]],sep=","),paste(out[[10]][[1]],out[[10]][[2]],sep=","),paste(out[[13]][[1]],out[[13]][[2]],sep=",")))
write(tab,file="n200p500.txt",sep=",")

