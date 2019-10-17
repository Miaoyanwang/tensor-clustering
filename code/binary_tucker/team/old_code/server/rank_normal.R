source('functions_synthesis_all_miaoyan.R')
#library(ggplot2)

### all covariates
seed=24
d=c(20,40)
dist="normal"


rank_range=array(0,dim=c(3,5,3))
rank_range[1,,]=rbind(c(3,3,3),c(3,3,2),c(3,2,2),c(2,2,2),c(3,2,3))
rank_range[2,,]=rbind(c(4,4,6),c(4,4,5),c(3,3,6),c(3,3,5),c(3,4,6))
rank_range[3,,]=rbind(c(6,8,8),c(6,7,7),c(6,7,8),c(5,7,7),c(5,8,8))
core_shape=rbind(c(3,3,3),c(4,4,6),c(6,8,8))

#############

#############
dup=10
final=array(0,dim=c(3,2,dup,3))
for(s in 1:3){
    for(i in 1:2){
        whole_shape = c(d[i],d[i],d[i])
        p=0.4*whole_shape
        data=gene_data_all(seed, whole_shape = whole_shape, core_shape=core_shape[s,],p=c(12,12,12),dist=dist, dup=dup, signal=4)
        #rank = expand.grid(rank1,rank2,rank3)
        #rank=rank[which(apply(rank,1,function(x){max(x)<= sqrt(prod(x))})==1),]
        #rank_range=rank[which(apply(rank,1,function(x){sum(x<=p)==3})),]
        table=NULL
        
        for(j in 1:dup){
            res=sele_rank(data$tsr[[j]],X_covar1 = data$X_covar1, X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,rank_range=rank_range[s,,],Nsim=10,cons = 'non',dist)
            table=rbind(table,res$rank)
        }
        
        final[s,i,,]=table
    }
}


save(final,file="normal_rank.RData")
