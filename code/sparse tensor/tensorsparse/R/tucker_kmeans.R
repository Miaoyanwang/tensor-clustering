#' tensor clustering by performing k-means on tucker decomposition
#' 
#' tensor clustering by performing k-means on tucker decomposition
#' @param x a three-dimensional array
#' @param k ...
#' @param r ...
#' @param l ...
#' 
#' @return a list with only one element: judgeX
#' 



tucker_means = function(x,k,r,l){
    tucker_result = tucker(as.tensor(data$x),c(k,r,l))
    fitted=attributes(tucker_result$est)$data
    
    
    n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
    mus = array(rep(0,n*p*q),c(n,p,q))
    Cs = kmeans(cp_result$U[[1]],k)$cluster
    Ds = kmeans(cp_result$U[[2]],r)$cluster
    Es = kmeans(cp_result$U[[3]],l)$cluster
    
    for(i in 1:k){
        for(j in 1:r){
            for(k in 1:l){
                mu[Cs==i,Ds==j,Es==l]=mean(x[Cs==1,Ds==j,Es==l])
            }
        }
    }
    # for (s in 1:multiplicative){
    # m1 = 1:n; m2 = 1:p; m3 = 1:q
    # for (i in 1:k) m1[Cs == i] = mean(cp_result$U[[1]][,s][Cs==i])
    # for (i in 1:r) m2[Ds == i] = mean(cp_result$U[[2]][,s][Ds==i])
    # for (i in 1:l) m3[Es == i] = mean(cp_result$U[[3]][,s][Es==i])
    # m1 = new("Tensor",3L,c(dim(x)[1],1L,1L),data=m1)
    # m2 = matrix(m2,ncol=1)
    #  m3 = matrix(m3,ncol=1)
    #  mus = mus + lambda[s]*ttm(ttm(m1,m2,2),m3,3)@data
    #}
    return(list(judgeX=fitted,Cs=Cs,Ds=Ds,Es=Es,blockmean=mus))
}
