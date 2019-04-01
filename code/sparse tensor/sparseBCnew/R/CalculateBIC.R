CalculateBIC <- function(x,bires,method="L0"){
    RSS=log(sum((x-bires$mus)^2))*nrow(x)*ncol(x)
    reducedCs=ifelse(dim(bires$Mus)[1]>0,dim(bires$Mus)[1],0)
    reducedDs=ifelse(dim(bires$Mus)[2]>0,dim(bires$Mus)[2],0)
    
    df_clustering=min(reducedCs*log(nrow(x))+reducedDs*log(ncol(x)),log(nrow(x)*ncol(x))*(sum(bires$Mus!=0)))
    
    df = log(nrow(x) * ncol(x)) * (sum(bires$Mus != 0) + df_clustering)
    
    df=df_clustering
    return(RSS+df)
    
    
    #if(method=="L1"){
    #RSS=log(sum((x-bires$mus)^2))*nrow(x)*ncol(x)
    #df=log(nrow(x)*ncol(x))*sum(bires$Mus!=0)
    #return(RSS+df)
    #} 
    #if(method=="L0") {
    #mat <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    #Cs <- bires$Cs
    #Ds <- bires$Ds
    #for(i in unique(Cs)){
    #  for(j in unique(Ds)){
    #    if(bires$Mus[i,j]!=0){
    #      mat[Cs==i,Ds==j] <- mean(x[Cs==i,Ds==j])
    #    }
    #  }
    #}
    #mat[abs(bires$mus)<=1e-8] <- mean(x[abs(bires$mus)<=1e-8])
    #return(log(sum((x-mat)^2))*nrow(x)*ncol(x) + log(nrow(x)*ncol(x))*sum(bires$Mus!=0))
    # } else {stop("No such kind of method:", method,".\n")}
}
