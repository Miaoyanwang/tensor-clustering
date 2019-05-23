#' Calculate BIC
#' 
#' Calculate BIC
#' @param x a three-dimensional array
#' @param clusterobj the return object of label2()
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @param apply "main": apply in the main formula; "cp": apply in the CPD k-means.
#' @return a vector [1]BIC, [2]nonzeromus
#' 
#' @export
tensor.calculateBIC = function (x, clusterobj,method="L0",apply="main") 
{
    if(apply!="main" & apply!="cp") stop("parameter apply does not take such a value named", apply, ".\n")
    
    ## Modified; add the degree of freedom due to clustering. 
    npq = dim(x)
    n = npq[1]; p = npq[2]; q = npq[3]
    RSS=log(sum((x-clusterobj$judgeX)^2))*(n*p*q)
    if (apply == "main"){
      reducedCs=ifelse(dim(clusterobj$mus)[1]>1,dim(clusterobj$mus)[1],0)
      reducedDs=ifelse(dim(clusterobj$mus)[2]>1,dim(clusterobj$mus)[2],0)
      reducedEs=ifelse(dim(clusterobj$mus)[3]>1,dim(clusterobj$mus)[3],0)
      
      #df_clustering=min(reducedCs*log(n)+reducedDs*log(p)+reducedEs*log(q),sum(clusterobj$mus != 0)*log(n*p*q))
      df_clustering=reducedCs*log(n)+reducedDs*log(p)+reducedEs*log(q)
      
      df=log(n*p*q)*(sum(clusterobj$mus != 0)+df_clustering)
      #df=df_clustering
      return(RSS+df)
    }
   if (apply == "cp"){
     df = log(n*p*q)/n/p/q*(n+p+q-2)*clusterobj$s
     return(RSS+df)
   }
    
   
}
