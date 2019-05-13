#' Calculate BIC
#' 
#' Calculate BIC
#' @param x a three-dimensional array
#' @param clusterobj the return object of label2()
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicates Lasso penalty.
#' @return a vector [1]BIC, [2]nonzeromus
#' 
tensor.calculateBIC = function (x, clusterobj,method="L0") 
{
    
    ## add the degree of freedom due to clustering. 
    npq = dim(x)
    n = npq[1]; p = npq[2]; q = npq[3]
    RSS=log(sum((x-clusterobj$judgeX)^2))*(n*p*q)
    reducedCs=ifelse(dim(clusterobj$mus)[1]>1,dim(clusterobj$mus)[1],0)
    reducedDs=ifelse(dim(clusterobj$mus)[2]>1,dim(clusterobj$mus)[2],0)
    reducedEs=ifelse(dim(clusterobj$mus)[3]>1,dim(clusterobj$mus)[3],0)
    
    df=log(n*p*q)*(sum(clusterobj$mus != 0)+df_clustering)
    return(RSS+df)
    
    #if (method == "L1"){
    #npq = dim(x)
    #n = npq[1]; p = npq[2]; q = npq[3]
    #RSS=log(sum((x-clusterobj$judgeX)^2))*(n*p*q)
    #df=log(n*p*q)*sum(clusterobj$mus != 0)
    #return(RSS+df)
    #}
    #if (method == "L0"){
    #mat <- array(rep(0,length(c(x))), dim=dim(x))
    #Cs <- clusterobj$Cs
    #Ds <- clusterobj$Ds
    #Es <- clusterobj$Es
    #for (i in unique(Cs)) {
    # for (j in unique(Ds)) {
    #    for (m in unique(Es)) {
    #      if (clusterobj$mus[i, j, m] != 0) {
    #        mat[Cs==i, Ds==j, Es==m] <- mean(x[Cs==i, Ds==j, Es==m])
    #      }
    #    }
    #  }
    #}
    #mat[abs(clusterobj$judgeX) <= 1e-6] <- mean(x[abs(clusterobj$judgeX) <= 1e-6])
    # return(log(sum((x - mat)^2))*dim(x)[1]*dim(x)[2]*dim(x)[3]+log(dim(x)[1]* 
    #                                                                   dim(x)[2]*dim(x)[3]) * sum(clusterobj$mus!=0))
    # }
    # stop("No such kind of method:", method, ".\n")
}
