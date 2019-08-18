#library("parallel")
#library("rTensor")
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
  if(method=="L1"){
  for (k in uniqCs){
    for (r in uniqDs){
      for (l in uniqEs){
        if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
        if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/(2*sum(Cs==k)*sum(Ds==r)*sum(Es==l)),method=method)
        if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
      }
    }
  }
  }## added 
  if(method=="L0"){
      for (k in uniqCs){
          for (r in uniqDs){
              for (l in uniqEs){
                  if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
                  if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),sqrt(lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l))),method=method)
                  ### modified
                  ## Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/sum(Cs==k)*sum(Ds==r)*sum(Es==l),method=method)
                  if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
              }
          }
      }
  }## added
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


mse = function(bires, data){
    npq = dim(bires$judgeX)
    n = npq[1]; p = npq[2]; q = npq[3]
    return(sum((bires$judgeX-data$truthX)^2)/n/p/q)
} 



label_for_cp = function(multiplicative=1,x,k,r,l){
  cp_result = cp(as.tensor(x),num_components = multiplicative)
  lambda = cp_result$lambdas
  fitted=attributes(cp_result$est)$data
  

  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
  mus = array(rep(0,n*p*q),c(n,p,q))
  Cs = kmeans(cp_result$U[[1]],k)$cluster
  Ds = kmeans(cp_result$U[[2]],r)$cluster
  Es = kmeans(cp_result$U[[3]],l)$cluster
  
  for(i in 1:k){
      for(j in 1:r){
          for(k in 1:l){
              mus[Cs==i,Ds==j,Es==k]=mean(x[Cs==1,Ds==j,Es==k])
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
  return(list(judgeX=fitted,s=multiplicative,Cs=Cs,Ds=Ds,Es=Es,blockmean=mus,mus=fitted))
}


cluster2block=function(mu,Cs,Ds,Es){
    d=c(length(Cs),length(Ds),length(Es))
    r=c(length(unique(Cs)),length(unique(Ds)),length(unique(Es)))
    block=array(0,dim=d)
    for(i in 1:r[1]){
        for(j in 1:r[2]){
            for(k in 1:r[3]){
                if(mu[i,j,k]==0){
                    block[which(Cs==i),which(Ds==j),which(Es==k)]="0 0 0"
                }
                else{
               block[which(Cs==i),which(Ds==j),which(Es==k)]=paste(unique(Cs)[i],unique(Ds)[j],unique(Es)[k]) 
                }
            }
        }
    }
    return(block)
    }


tensor.unfold4 = function(tensor,dim=1){
  return(t(apply(tensor,dim,c)))
}


Objective4 = function (x, mu.array, Cs, Ds, Es, Fs, lambda = 0, method="L0") {
  if(method=="L0") return(sum((x - mu.array[Cs, Ds, Es, Fs, drop=FALSE])^2)+2*lambda*sum(mu.array !=0 ))
  if(method=="L1") return(sum((x - mu.array[Cs, Ds, Es, Fs, drop=FALSE])^2)+2*lambda*sum(abs(mu.array)))
  stop("No such kind of method:", method,".\n")
}

UpdateMus.tensor4 = function (x, Cs, Ds, Es, Fs, lambda=0, method="L0") {
  uniqCs = sort(unique(Cs))
  uniqDs = sort(unique(Ds))
  uniqEs = sort(unique(Es))
  uniqFs = sort(unique(Fs))
  mus = array(NA, c(length(uniqCs), length(uniqDs), length(uniqEs), length(uniqFs)))
  if(method=="L1"){
    for (k in uniqCs){
      for (r in uniqDs){
        for (l in uniqEs){
          for(m in uniqFs){
            if (lambda == 0) mus[k,r,l,m] = mean(x[Cs==k,Ds==r,Es==l,Fs=m])
            if (lambda > 0) mus[k,r,l,m] = Soft(mean(x[Cs==k,Ds==r,Es==l,Fs=m]),lambda/(2*sum(Cs==k)*sum(Ds==r)*sum(Es==l)*sum(Fs==m)),method=method)
            if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
          }
        }
      }
    }
  }## added 
  if(method=="L0"){
    for (k in uniqCs){
      for (r in uniqDs){
        for (l in uniqEs){
          for (m in uniqFs) {
            if (lambda == 0) mus[k,r,l,m] = mean(x[Cs==k,Ds==r,Es==l,Fs==m])
            if (lambda > 0) mus[k,r,l,m] = Soft(mean(x[Cs==k,Ds==r,Es==l,Fs==m]),sqrt(lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l)*sum(Fs==m))),method=method)
            ### modified
            ## Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/sum(Cs==k)*sum(Ds==r)*sum(Es==l),method=method)
            if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
          }
        }
      }
    }
  }## added
  return(mus)
}
