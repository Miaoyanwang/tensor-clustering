#################################################
# May 7 2013
# Table 7 for MVN  bicluster unknown
# We split the 50 iterations and run it separately since it is computational intensive
# The R-code presented here is for one iteration
#################################################
# Simulation for MVN using BIC Lambda
#source("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/Functions.R")
#source("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/MatrixBiclustermodified.R")
#source("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/MVNLambdaUnknown.R")
source("MVNLambdaUnknown.R")
library(glasso)

args<-commandArgs(TRUE)
print(args)
iteration <- as.numeric(args[[1]])
n=200
p=200
sd=2
lambda=c(0,10,20,30,40,50,60,70)
res1<-Do(n=n,p=p,k=4,r=5,iteration=iteration,lambda=lambda,standarddeviation=sd,max.iter=50)

#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/cerrow.",iteration,sep=""),res1$cerrow,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/cercol.",iteration,sep=""),res1$cercol,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/totalzero.",iteration,sep=""),res1$totalzero,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/correctzero.",iteration,sep=""),res1$correctzero,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/correctone.",iteration,sep=""),res1$correctone,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/totalincorrect.",iteration,sep=""),res1$totalincorrect,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

#write.table(file=paste("/home/students/keanming/Dropbox/Sparse_Biclustering/April_2013_Simulations/MVN_MSE/May2013sd2unc/AutoBICUnknown/n200p200/selectlambda",iteration,sep=""),res1$selectlambda,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)






