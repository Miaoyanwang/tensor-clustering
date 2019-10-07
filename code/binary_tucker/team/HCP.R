.libPaths( c( .libPaths(), "/workspace/miaoyan/x86_64-pc-linux-gnu-library/3.4") )
source('functions_synthesis_all_miaoyan.R')

##BIC=sele_rank(tsr,X_covar1 , X_covar2 ,X_covar3 ,rank = 2:4, Nsim,cons = 'non',dist)


data=load("../../../data/binary_tucker/HCP.RData")
load("output_HCP.RData")
tsr=tensor
table(attr[,5])
#####
##22-25 26-30 31+
##35    58    43 

X=attr[,4:5]
levels(X[,2])=c("22-25","26-30","31+","31+") ## three groups

contrasts(X[,1]) <- contr.sum 
contrasts(X[,2]) <- contr.sum 

X_covar3=model.matrix(~as.factor(X[,1])+as.factor(X[,2])) ## baseline female and age 22-25

### intercept, F 1, M,-1; 22-25 [1,0]; 26-30 [0 1]
X_covar1=X_covar2=NULL

#intercept=sqrt(sum(X_covar3[,1]^2))
#gender=sqrt(sum(X_covar3[,2]^2))
#young=sqrt(sum(X_covar3[,3]^2))
#mid=sqrt(sum(X_covar3[,4]^2))
##BIC=sele_rank(tsr,NULL, NULL ,X_covar3 ,rank1 = 4:11,rank2 = 4:11,rank3 = 2:4, Nsim,cons = 'non',dist)

massive=massive_glm(tsr,X_covar3,"binary")
coef=result$C_ts
old=-coef[,,3]-coef[,,4]
#####

core_shape=c(10,10,4)
Nsim=50
cons="non"
lambda=0.1
alpha=10
solver="GC"
dist="binary"

result=update_binary_all(tsr,NULL,NULL,X_covar3,core_shape,Nsim,cons,lambda,alpha,solver,dist)

save(result,file="output.RData")
coef_old=coef

coef=mass=array(0,dim=c(68,68,5))

coef[,,1]=coef_old[,,1]
coef[,,2]=gender
coef[,,3]=age1
coef[,,4]=age2
coef[,,5]=age3

mass[,,1]=massive[,,1]
mass[,,2]=massive[,,2]
mass[,,3]=massive[,,3]
mass[,,4]=massive[,,4]
mass[,,5]=-massive[,,4]-massive[,,3]

m=tensor[,,1]
for(i in 2:136){
m=m+tensor[,,i]
}
m=m/136
m=1*((m==1)|(m==0))


pdf("compare.pdf",width=8,height=11)
par(mfrow=c(4,2))
hist(c(mass[,,1][m!=1]),nclass=40,xlab="Effect size",main="Intercept (Classical GLM)")
hist(c(coef[,,1][m!=1]),nclass=40,xlab="Effect size",main="Intercept (Tensor regression)")

hist(c(mass[,,2][m!=1]),nclass=40,xlab="Effect size",main="Gender (Classical GLM)")
hist(c(coef[,,2][m!=1]),nclass=40,xlab="Effect size",main="Gender (Tensor regression)")

hist(c(mass[,,4][m!=1]),nclass=40,xlab="Effect size",main="Age 26-30 (Classical GLM)")
hist(c(coef[,,4][m!=1]),nclass=40,xlab="Effect size",main="Age 26-30 (Tensor regression)")

hist(c(mass[,,5][m!=1]),nclass=40,xlab="Effect size",main="Age 31+ (Classical GLM)")
hist(c(coef[,,5][m!=1]),xlab="Effect size",nclass=40,main="Age 31+ (Tensor regression)")

coef_new=array(0,dim=c(68,68,5))
for(i in 2:5){
    x=coef[,,i]
    x[m==1]=0
    coef_new[,,i]=x
}
list=NULL
s=5
for(i in 1:68){
    for(j in 1:68){    
list=rbind(list,c(i,j,coef_new[i,j,s]))
    }
}

list[sort(list[,3],index=T)$ix,]

