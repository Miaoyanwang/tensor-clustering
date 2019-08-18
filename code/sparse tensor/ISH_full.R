install.packages("tensorsparse_1.0.tar.gz",repos=NULL,type="source")
library("tensorsparse")

require("sparseBCnew")
require("tensorsparse")
if(!require("clues")){
  install.packages("clues")
  stopifnot(require("clues"))
}

ISH=read.csv("../../data/sparse_tensor/entrezlist_ISH_full.csv")
position = ISH$rowid
genefamily_truth = ISH$gene_family

load("../../data/sparse_tensor/brain_tensor.rdata")
rawdata1 = data[[1]]
sample = rawdata1[position, , ]
sample=sample-mean(sample)
save(sample,file="../../data/sparse_tensor/brain_final.rdata")

load("../../data/sparse_tensor/brain_final.rdata")

n = 362; p = 193; q = 13; #dimension
k = 5; #gene clusters, set as 6 as the top 6 gene types from ISH data
r = 2; #individual clusters
l = 4; #tissue cluster

lambda_0 =c(floor(n*p*q/k/r/l))

method = "L0"
if (method == "L0") lambda =  sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.05)
if (method == "L1") lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.05)

verlam = chooseLambda(sample,k,r,l,lambda=lambda,method=method)
cat("the lambda tensorsparse choose is", verlam$lambda,".\n")
cat("the nonzeromus calculated by tensorsparse is ", verlam$nonzeromus, ".\n")
cat("=================================================================\n")

sim = label2(sample,k,r,l,threshold=1e-10,lambda=verlam$lambda,sim.times=5,trace=FALSE,method=method)
judgeX = sim$judgeX

#load each of the R file in tensorsparse/R/
#then choose best krl by: res
res = choosekrl_bic(sample, 1:10,1:5,1:10)
res[1]

krange = 1:50
rrange = 1:10
lrange = 1:10
kres = choosekrl_bic(sample, krange, 5, 4)
kres1 = as.vector(kres[2]$BIC)
rres = choosekrl_bic(sample, 10, rrange, 5)
rres1 = as.vector(rres[2]$BIC)
lres = choosekrl_bic(sample, 10, 5, lrange)
lres1 = as.vector(lres[2]$BIC)

plot(krange, kres1, type = 'p', main = 'K range 1-50 VS score')
#best 9,26
plot(rrange, rres1, type = 'p', main = 'R range 1-10 VS score')
#best 5
plot(lrange, lres1, type = 'p', main = 'L range 1-10 VS score')
#best 4

k = 26; #gene clusters, set as 6 as the top 6 gene types from ISH data
r = 5; #individual clusters
l = 4; #tissue cluster

lambda=chooseLambda(sample,26,5,4,method="L0") $lambda=71.048
result=classify2(sample,k,r,l,lambda$lambda)


t=1
which(result$mus[,,t]>=max(result$mus[,,t])-4,arr.ind=T)
max(result$mus[1,,t])
gene[which(result$Cs==1),4]
names(tensor[1,1,])[which(result$Es==t)]


which(result$mus[,,t]<=(min(result$mus[,,t])+2),arr.ind=T)
min(result$mus[23,,t])
gene[which(result$Cs==23),4]



######
library(R.matlab)
data=readMat("../../data/sparse_tensor/dnations.mat")
sample=data$R
sample[is.na(sample)]=0

method = "L0"
lambda = seq(0,2^2,by=0.05)
k=4
r=4
l=12

verlam = chooseLambda(sample,k,r,l,lambda=lambda,method=method)
cat("the lambda tensorsparse choose is", verlam$lambda,".\n")
cat("the nonzeromus calculated by tensorsparse is ", verlam$nonzeromus, ".\n")
cat("=================================================================\n")

###
unlist(data$relnnames)
unlist(data$countrynames)

table(res$Cs)
unlist(data$countrynames)[sort(res$Cs,index=T,decreasing=F)$ix]
table(res$Ds)
unlist(data$countrynames)[sort(res$Ds,index=T,decreasing=F)$ix]
###

##sparse
res = label2(sample,k,r,l,threshold=1e-10,lambda=verlam$lambda,sim.times=5,trace=FALSE,method=method)
validate(res$Cs,res$Ds,res$Es,sample,res)

##non-sparse
res = label2(sample,k,r,l,threshold=1e-10,lambda=0,sim.times=5,trace=FALSE,method=method)
validate(res$Cs,res$Ds,res$Es,sample,res)


### combinatorial 
matT<-unfold(as.tensor(sample),row_idx=1,col_idx=c(2,3))
Cs=kmeans(attr(matT,"data"),k)$cluster
matT<-unfold(as.tensor(sample),row_idx=2,col_idx=c(1,3))
Ds=kmeans(attr(matT,"data"),r)$cluster
matT<-unfold(as.tensor(sample),row_idx=3,col_idx=c(1,2))
Es=kmeans(attr(matT,"data"),l)$cluster
validate(Cs,Ds,Es,sample)

### CP
decomp=cp(as.tensor(sample),max(k,r,l))
Cs=kmeans(decomp$U[[1]],k)$cluster
Ds=kmeans(decomp$U[[2]],r)$cluster
Es=kmeans(decomp$U[[3]],l)$cluster
validate(Cs,Ds,Es,sample)

### Tucker
decomp=tucker(as.tensor(sample),c(k,r,l))
Cs=kmeans(decomp$U[[1]],k)$cluster
Ds=kmeans(decomp$U[[2]],r)$cluster
Es=kmeans(decomp$U[[3]],l)$cluster
validate(Cs,Ds,Es,sample)

#### brain 

validate=function(Cs,Ds,Es,sample,res=NULL){
    index=array(0,dim=dim(sample))
    k=max(Cs)
    r=max(Ds)
    l=max(Es)
    
    for(i in 1:k){
        for(j in 1:r){
            for(m in 1:l){
                if(length(res)!=0){
                    if(res$mus[i,j,m]==0) index[which(Cs==i),which(Ds==j),which(Es==m)]=0
                    else
                    index[which(Cs==i),which(Ds==j),which(Es==m)]=paste(i,j,m)
                }
                else
                index[which(Cs==i),which(Ds==j),which(Es==m)]=paste(i,j,m)
            }
        }
    }
    block_id=as.factor(c(index))
    levels(block_id)=1:(k*r*l)
    
    fit=lm(c(sample)~block_id)
    return(summary(fit)$r.squared)
    #dist=as.matrix(dist(c(sample)))
    #return(summary(silhouette(as.numeric(block_id),dist))) ## 0.598035 
}

