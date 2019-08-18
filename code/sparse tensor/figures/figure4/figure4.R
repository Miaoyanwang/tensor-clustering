rm(list=ls())
require("tensorsparse")
require(ggplot2)

########################################
###   Figure 4    #######################
########################################


set.seed(1)
cat("The first situation:\n")
n = 50; p = 50; q = 50; k = 4; r = 4; l = 4; s = 3
our_cerC = 1:50; cp_cerC = 1:50
our_cerD = 1:50; cp_cerD = 1:50
our_cerE = 1:50; cp_cerE = 1:50
our_cerC_mean = 1:3;our_cerD_mean = 1:3;our_cerE_mean = 1:3
cp_cerC_mean = 1:3;cp_cerD_mean = 1:3;cp_cerE_mean = 1:3
our_cerC_sd = 1:3;our_cerD_sd = 1:3;our_cerE_sd = 1:3
cp_cerC_sd = 1:3;cp_cerD_sd = 1:3;cp_cerE_sd = 1:3
for (i in 1:3){
  cat("When the noise is", c(0,10,20)[i], ":\n")
  for (iter in 1:50){
    error = c(0,10,20)[i]
    try = get.data(n,p,q,k,r,l,error=error,multiplicative = s)
    our_result = classify2(try$x,k,r,l)
    our_ev = sparse.evaluate(our_result,try,CER=TRUE,show=FALSE)
    our_cerC[iter] = our_ev$cerC;our_cerD[iter] = our_ev$cerD;our_cerE[iter] = our_ev$cerE
    
    cp_result = cp_kmeans(try$x,k,r,l,multiplicative = s)
    cp_ev = sparse.evaluate(cp_result,try,CER=TRUE,show=FALSE)
    cp_cerC[iter] = cp_ev$cerC;cp_cerD[iter] = cp_ev$cerD;cp_cerE[iter] = cp_ev$cerE
  }
  our_cerC_mean[i] = mean(our_cerC); our_cerD_mean[i] = mean(our_cerD); our_cerE_mean[i] = mean(our_cerE)
  cp_cerC_mean[i] = mean(cp_cerC); cp_cerD_mean[i] = mean(cp_cerD); cp_cerE_mean[i] = mean(cp_cerE)
  our_cerC_sd[i] = sd(our_cerC); our_cerD_sd[i] = sd(our_cerD); our_cerE_sd[i] = sd(our_cerE)
  cp_cerC_sd[i] = sd(cp_cerC); cp_cerD_sd[i] = sd(cp_cerD); cp_cerE_sd[i] = sd(cp_cerE)
  cat("The cerC of our method is:", round(mean(our_cerC),4),"(", round(sd(our_cerC),4),");\n")
  cat("The cerD of our method is:", round(mean(our_cerD),4),"(", round(sd(our_cerD),4),");\n")
  cat("The cerE of our method is:", round(mean(our_cerE),4),"(", round(sd(our_cerE),4),");\n")
  cat("the cerC of cp k-means is:", round(mean(cp_cerC), 4),"(", round(sd(cp_cerC),4), ");\n")
  cat("the cerD of cp k-means is:", round(mean(cp_cerD), 4),"(", round(sd(cp_cerD),4), ");\n")
  cat("the cerE of cp k-means is:", round(mean(cp_cerE), 4),"(", round(sd(cp_cerE),4), ").\n")
}
figure4_1 = data.frame(cerC_mean = c(our_cerC_mean,cp_cerC_mean), 
                     method = c("us","us","us","CPD k-means","CPD k-means","CPD k-means"),
                     noise = c(0,10,20,0,10,20))

pdf("multidata.pdf",width=5,height=3)
ggplot(data=figure4_1, aes(x=noise,y=cerC_mean))+geom_line(aes(color=method))
dev.off()

set.seed(2)
cat("The second situation:\n")
n = 50; p = 50; q = 50; k = 4; r = 4; l = 4
our_cerC = 1:50; cp_cerC = 1:50
our_cerD = 1:50; cp_cerD = 1:50
our_cerE = 1:50; cp_cerE = 1:50
our_cerC_mean = 1:3;our_cerD_mean = 1:3;our_cerE_mean = 1:3
cp_cerC_mean = 1:3;cp_cerD_mean = 1:3;cp_cerE_mean = 1:3
our_cerC_sd = 1:3;our_cerD_sd = 1:3;our_cerE_sd = 1:3
cp_cerC_sd = 1:3;cp_cerD_sd = 1:3;cp_cerE_sd = 1:3
for (i in 1:3){
  cat("When the noise is", c(0,10,20)[i], ":\n")
  for (iter in 1:50){
    error = c(0,10,20)[i]
    try = get.data(n,p,q,k,r,l, error=error)
    our_result = classify2(try$x,k,r,l)
    our_ev = sparse.evaluate(our_result,try,CER=TRUE,show=FALSE)
    our_cerC[iter] = our_ev$cerC;our_cerD[iter] = our_ev$cerD;our_cerE[iter] = our_ev$cerE
    cp_result = cp_kmeans(try$x,k=k,r=r,l=l,multiplicative = NULL,max.s=10)
    cp_ev = sparse.evaluate(cp_result,try,CER=TRUE,show=FALSE)
    cp_cerC[iter] = cp_ev$cerC;cp_cerD[iter] = cp_ev$cerD;cp_cerE[iter] = cp_ev$cerE
  }
  our_cerC_mean[i] = mean(our_cerC); our_cerD_mean[i] = mean(our_cerD); our_cerE_mean[i] = mean(our_cerE)
  cp_cerC_mean[i] = mean(cp_cerC); cp_cerD_mean[i] = mean(cp_cerD); cp_cerE_mean[i] = mean(cp_cerE)
  our_cerC_sd[i] = sd(our_cerC); our_cerD_sd[i] = sd(our_cerD); our_cerE_sd[i] = sd(our_cerE)
  cp_cerC_sd[i] = sd(cp_cerC); cp_cerD_sd[i] = sd(cp_cerD); cp_cerE_sd[i] = sd(cp_cerE)
  cat("The cerC of our method is:", round(mean(our_cerC),4),"(", round(sd(our_cerC),4),");\n")
  cat("The cerD of our method is:", round(mean(our_cerD),4),"(", round(sd(our_cerD),4),");\n")
  cat("The cerE of our method is:", round(mean(our_cerE),4),"(", round(sd(our_cerE),4),");\n")
  cat("the cerC of cp k-means is:", round(mean(cp_cerC), 4),"(", round(sd(cp_cerC),4), ");\n")
  cat("the cerD of cp k-means is:", round(mean(cp_cerD), 4),"(", round(sd(cp_cerD),4), ");\n")
  cat("the cerE of cp k-means is:", round(mean(cp_cerE), 4),"(", round(sd(cp_cerE),4), ").\n")
}
figure4_2 = data.frame(cerC_mean = c(our_cerC_mean,cp_cerC_mean), 
                     method = c("us","us","us","CPD k-means","CPD k-means","CPD k-means"),
                     noise = c(0,10,20,0,10,20))

pdf("constantdata.pdf",width=5,height=3)
ggplot(data=figure4_2, aes(x=noise,y=cerC_mean))+geom_line(aes(color=method))
dev.off()






