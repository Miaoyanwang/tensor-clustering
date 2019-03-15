matrixBC.BIC <-
function(x,k,r,lambda,alpha=0.2,beta=0.2,nstart=20,Sigma.init=NULL,Delta.init=NULL){
    x<-x-mean(x)
    BIC<-rep(NA,length(lambda))
    nonzero<-rep(NA,length(lambda))
    
    for(i in 1:length(lambda)){
        mclustering<-matrixBC(x,k,r,lambda=lambda[i],nstart=nstart,alpha=alpha,beta=beta,Sigma.init=Sigma.init,Delta.init=Delta.init)
        BIC[i]<-CalculateBICMatrix(x,mclustering)
        nonzero[i]<-sum(mclustering$Mus!=0)
        
    }
    return(list(lambda=lambda[which(BIC==min(BIC))[1]],BIC=BIC,nonzeromus=nonzero))
}


############################################
# Function to calculate BIC as in Section 5.2
############################################
CalculateBIC <- function(x,biclustobj){
    mat <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    Cs <- biclustobj$Cs
    Ds <- biclustobj$Ds
    for(i in unique(Cs)){
        for(j in unique(Ds)){
            if(biclustobj$Mus[i,j]!=0){
                mat[Cs==i,Ds==j] <- mean(x[Cs==i,Ds==j])
            }
        }
    }
    mat[biclustobj$mus==0] <- mean(x[biclustobj$mus==0])
    return(log(sum((x-mat)^2))*nrow(x)*ncol(x) + log(nrow(x)*ncol(x))*sum(biclustobj$Mus!=0))
}
