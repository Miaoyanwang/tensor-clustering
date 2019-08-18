#' Perform tensor clustering via cp decomposition
#' 
#' Perform tensor clustering via cp decomposition.
#' @param x a three-dimensional array; data tensor
#' @param k \eqn{d_1}: the clusters number of mode 1
#' @param r \eqn{d_2}: the clusters number of mode 2
#' @param l \eqn{d_3}: the clusters number of mode 3
#' @param multiplicative the number of components
#' @param max.s the max value when selecting s
#' 
#' @return a list with only one element: judgeX
#' 
#' @export


cp_kmeans= function(x,k,r,l,multiplicative=NULL,max.s=NULL){

  #krl = matrix(c(rep(1:length(range.k),each=length(range.r)*length(range.l)),
  #               rep(1:length(range.r),times=length(range.k)*length(range.l)),
  #               rep(rep(1:length(range.l),each=length(range.r)),times=length(range.k))),byrow=TRUE,
  #             nrow=3)
  #krl_list = as.list(as.data.frame(krl))
  
  #if (.Platform$OS.type == "windows") {
  #bires = apply(krl,MARGIN=2,label_for_cp,x=x,multiplicative=multiplicative)
  #CBIC = lapply(bires,tensor.calculateBIC, x=x)
  #BIC = unlist(CBIC)
  #} else {
  #  bires = mclapply(krl_list, label_for_cp, x=x, multiplicative=multiplicative, mc.cores=n.cores)
  #  CBIC = mclapply(bires,tensor.calculateBIC, x=x, mc.cores = n.cores)
  #  BIC = unlist(CBIC)
  #}
  #names(BIC) = apply(krl,MARGIN=2,FUN=function(x) paste(range.k[x[1]],range.r[x[2]],range.l[x[3]]))
  #best = krl_list[[which(BIC == min(BIC))[1]]]
  
  if (is.null(multiplicative)){
    multiplicative = 1:max.s
    bires = lapply(multiplicative,label_for_cp,x=x,k=k,r=r,l=l)
    CBIC = lapply(bires,tensor.calculateBIC,x=x,apply="cp")
    BIC = unlist(CBIC)
    best = which(BIC == max(BIC))
  } else {best = multiplicative}

  return(label_for_cp(best,x,k,r,l))
}

