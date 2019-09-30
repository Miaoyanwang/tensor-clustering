source('functions_synthesis_all_miaoyan.r')

seed = 24;whole_shape = c(68,68,136); core_shape = c(2,2,3);p=c(-1,2,-1);signal=10;dup=2;dist="binary";

data = gene_data_all(seed,whole_shape , core_shape, p, dist,dup,signal)

cons="non"; lambda=0.1; alpha=7; solver="GC";Nsim=3; solver="CG";dist="binary";
tsr=data$tsr[[1]];X_covar1=data$X_covar1;X_covar2=data$X_covar2;X_covar3=data$X_covar3;

result=update_binary_all(tsr,X_covar1,X_covar2,X_covar3,core_shape,Nsim,cons,lambda,alpha,solver,dist)
plot(result$lglk)

plot(result$C_ts,data$C_ts)
plot(result$U,data$U)

BIC=sele_rank(tsr,X_covar1 , X_covar2 ,X_covar3 ,rank = 2:4, Nsim,cons = 'non',dist)
####-------------------------------  data analysis ####-------------------------------  
data=load("../../../data/binary_tucker/HCP.RData")
tsr=tensor
table(attr[,5])
#####
##22-25 26-30 31+
##35    58    43 

X=attr[,4:5]
levels(X[,2])=c("22-25","26-30","31+","31+") ## three groups

X_covar3=model.matrix(~-1+as.factor(X[,1])+as.factor(X[,2])) ## baseline female and age 22-25

X_covar1=X_covar2=NULL


BIC=sele_rank(tsr,NULL, NULL ,X_covar3 ,rank1 = 8:10,rank2 = 8:10,rank3 = 3:4, Nsim,cons = 'non',dist)


result=update_binary_all(tsr,NULL,NULL,X_covar3,c(10,10,4),Nsim,cons,lambda,alpha,solver,dist)
## baseline 22-25 and 


####-------------------------------  convergence rate
p = seq(5,19,2)
re = c()

for(i in seq(length(p))){
  con = conv_rate(seed = 24,d = seq(20,70,10),r = rep(3,7), p1 = rep(p[i],7), p2 = rep(p[i],7), 
                  dis = 'gaussian',gs_mean = 0,gs_sd = 10, 
                   dup=2, Nsim =50, linear = FALSE, cons = 'penalty' ,lambda = 1,
                   solver = 'CG')
  
  con = as.data.frame(con)
  con$d = seq(20,70,10); con$rank = rep(3,7); con$p1 = rep(p[i],7); con$p2 = rep(p[i],7)
  re = rbind(re,con)
}


write.csv(re,"rank3.csv",row.names = FALSE)

## visual
re = read.csv('rank3.csv',header = TRUE)

# re$RMSE = re$RMSE*sqrt(re$d^3)
# re$rate = sqrt(re$rate)

library(ggplot2)

#re[re$d != 20,]
ggplot(re, aes(x = d, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 1.5)  +
  geom_point( size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggplot(re[re$d != 20 ,], aes(x = d, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 1.5)  +
  geom_point( size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggplot(re, aes(x = p1, y = RMSE)) + geom_line(aes(color = as.factor(d)),size = 1.5)  +
  geom_point(aes(shape =  as.factor(d)), size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(re[re$d != 20,], aes(x = p1, y = RMSE)) + geom_line(aes(color = as.factor(d)),size = 1.5)  +
  geom_point(aes(shape =  as.factor(d)), size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(re, aes(x = rate, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 2)  +
  geom_point(size = 3)  + 
  geom_smooth(method = 'lm',formula = y~x) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(re[re$d != 20 ,], aes(x = rate, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 2)  +
  geom_point( size = 3)  + 
  geom_smooth(method = 'loess',formula = y~x) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))




