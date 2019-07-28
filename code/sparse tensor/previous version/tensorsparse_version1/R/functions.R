n.cores = detectCores()
ReNumber2 = function(Cs){
  #Cs=c(2,2)
  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  num = c()
  for (i in 1:length(uniq)) {
    num[i] = sum(Cs==i)
  }
  newCs[order(Cs)] = rep(order(num),num)
  return(newCs)
}

reorderClusters = function(x,Cs,Ds,Es){
  #x=truthX;Cs=truthCs;Ds=truthDs;Es=truthEs
  Cs = ReNumber2(Cs);Ds = ReNumber2(Ds);Es = ReNumber2(Es)
  return(list("x"=x[order(Cs),order(Ds),order(Es), drop=FALSE],"Cs"=sort(Cs),
              "Ds"=sort(Ds),"Es"=sort(Es)))
}

tensor.unfold = function(tensor,dim=1){
  if (dim == 1) unfold = aperm(tensor,c(3,2,1))
  if (dim == 2) unfold = aperm(tensor,c(1,3,2))
  if (dim == 3) unfold = tensor
  unfold = apply(unfold,3,c)
  if (is.vector(unfold)) return(as.matrix(unfold)) else return(t(unfold))
}

design.row = function(sp,C1,D1,E1,k,r,l){
  labels = array(rep(0,k*r*l),c(k,r,l))
  #print(C1[sp[1]])
  labels[C1[sp[1]],D1[sp[2]],E1[sp[3]]] = 1
  return(c(labels))
}

Objective = function (x, mu.array, Cs, Ds, Es, lambda = 0, method="L0") {
  if(method=="L0") return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2)+2*lambda*sum(mu.array !=0 ))
  if(method=="L1") return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2)+2*lambda*sum(abs(mu.array)))
  stop("No such kind of method:", method,".\n")
}

#dim should be a vector
tensor.index = function(index,dims){
  index = index-1
  Cs = (index %% dims[1])+1
  Ds = (index %% (dims[1]*dims[2]))%/%dims[1] +1
  Es = (index %/% (dims[1]*dims[2]))+1
  return(c(Cs,Ds,Es))
}


Soft = function (a, b, method="L0"){
  if (b < 0) 
    stop("Can soft-threshold by a nonnegative quantity only.")
  if(method == "L0") return(sign(a) * ifelse(abs(a)>b, abs(a), 0))
  if(method == "L1") return(sign(a) * pmax(0, abs(a) - b))
  stop("No such kind of method:", method,".\n")
}


#give x as an array
UpdateMus.tensor = function (x, Cs, Ds, Es, lambda=0, method="L0") {
  uniqCs = sort(unique(Cs))
  uniqDs = sort(unique(Ds))
  uniqEs = sort(unique(Es))
  mus = array(NA, c(length(uniqCs), length(uniqDs), length(uniqEs)))
  for (k in uniqCs){
    for (r in uniqDs){
      for (l in uniqEs){
        if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
        if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l)),method=method)
        if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
      }
    }
  }
  return(mus)
}

UpdateClusters.tensor = function (x, mus, curCs, curDs) {
  Cs.new <- rep(NA, length(curCs))
  uniq <- 1:max(curCs)
  mus.expandcolumns <- mus[, curDs, drop = FALSE]
  #if (dim(mus.expandcolumns)[1]==1) mus
  for (i in 1:nrow(x)) {
    dist2.clust <- NULL
    for (k in 1:length(uniq)) {
      #see which cluster is closest to one sample
      dist2.clust <- c(dist2.clust, sum((x[i, , drop = FALSE] - 
                                           mus.expandcolumns[k, , drop = FALSE])^2))
    }
    wh <- which(dist2.clust == min(dist2.clust))
    Cs.new[i] <- wh[1]
  }
  return(Cs.new)
}

ReNumber = function (Cs) 
{
  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  for (i in 1:length(uniq)) {
    newCs[Cs == uniq[i]] <- i
  }
  return(newCs)
}

label_for_krl = function(abc,k,r,l,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,sim.times,lambda,xmiss,method="L0",crossvalidation = TRUE){
  a = abc[1]; b = abc[2]; c = abc[3]
  if(is.null(Cs.init) == FALSE) Cs.init = Cs.init[,a]
  if(is.null(Ds.init) == FALSE) Ds.init = Ds.init[,b]
  if(is.null(Es.init) == FALSE) Es.init = Es.init[,c]
  if (crossvalidation == TRUE){
    return(list(label2(xmiss, k[a], r[b], l[c], lambda = lambda, 
                       Cs.init = Cs.init, Ds.init = Ds.init, 
                       Es.init = Es.init, sim.times=sim.times, method=method)$judgeX))} else {
                         return(label2(xmiss, k[a], r[b], l[c], lambda = lambda, 
                                       Cs.init = Cs.init, Ds.init = Ds.init, 
                                       Es.init = Es.init, sim.times=sim.times, method=method))
                       }
}

positionfun=function(d){
  x=rep(1,d[2])%x%rep(1,d[3])%x%c(1:d[1])
  y=rep(1,d[3])%x%(c(1:d[2])%x%rep(1,d[1]))
  z=c(1:d[3])%x%rep(1,d[1])%x%rep(1,d[2])
  position=cbind(x,y,z)
  return(list("position"=position))
}







