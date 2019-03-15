if(!require("glmnet")){
  install.packages("glmnet")
  stopifnot(require("glmnet"))
}
if(!require("rgl")){
  install.packages("rgl")
  stopifnot(require("rgl"))
}
if(!require("reshape")){
  install.packages("reshape")
  stopifnot(require("reshape"))
}
if(!require("mvtnorm")){
  install.packages("mvtnorm")
  stopifnot(require("mvtnorm"))
}
if(!require("HDCI")){
  install.packages("HDCI")
  stopifnot(require("HDCI"))
}
if(!require("clues")){
  install.packages("clues")
  stopifnot(require("clues"))
}



ReNumber2 = function(Cs){
  #Cs=truthDs
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
  #Cs=truthCs;Ds=truthDs;Es=truthEs
  Cs = ReNumber2(Cs);Ds = ReNumber2(Ds);Es = ReNumber2(Es)
  return(list("x"=x[order(Cs),order(Ds),order(Es)],"Cs"=sort(Cs),
              "Ds"=sort(Ds),"Es"=sort(Es)))
}

#sort=TRUE: we put the blocks with the same mean together
get.data = function(n,p,q,k,r,l,error=3,sort=TRUE,sparse.percent=0){
  #n=20;p=20;q=20;k=2;r=2;l=2;error=3
  mus = runif(k*r*l, -3,3)#take the mean of k*r*l biclusters/cubes
  if(sparse.percent!=0) mus[1:floor(k*r*l*sparse.percent)] = 0
  mus = array(mus,c(k,r,l))
  truthCs = sample(1:k,n,rep=TRUE)
  truthDs = sample(1:r,p,rep=TRUE)
  truthEs = sample(1:l,q,rep=TRUE)
  x = array(rnorm(n*p*q,mean=0,sd=error),c(n,p,q))
  truthX = array(rep(0,n*p*q),c(n,p,q))
  for(i in 1:max(truthCs)){
    for(j in 1:max(truthDs)){
      for(m in 1:max(truthEs)){
        x[truthCs==i, truthDs==j, truthEs==m] = x[truthCs==i, truthDs==j, truthEs==m] + mus[i,j,m]
        truthX[truthCs==i, truthDs==j, truthEs==m] =  mus[i,j,m]
      }
    }
  }
  #暂时去掉了中心化
  #x = x - mean(x)#minus the overall mean
  if(sort==TRUE){
    truthX = reorderClusters(truthX,truthCs,truthDs,truthEs)$x
    neworder = reorderClusters(x,truthCs,truthDs,truthEs)
  }
  result = list("x"=neworder$x,"truthX"=truthX,"truthCs"=neworder$Cs,"truthDs"=neworder$Ds,"truthEs"=neworder$Es,"mus"=mus)
  return(result)
}



tensor.unfold = function(tensor,dim=1){
  #dim=1: unfold by row; 2: by column; 3: by the 3rd dimension
  #tensor=mu.array;dim=3
  #tensor=guess;dim=3
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

#x is array, mus is array
Objective = function (x, mu.array, Cs, Ds, Es, lambda = 0) {
  #print(dim(x))
  #print(dim(mu.array[Cs, Ds, Es]))
  #print(Cs)
  #sum((x - mu.array[Cs, Ds, Es])^2)
  return(sum((x - mu.array[Cs, Ds, Es])^2)+2*lambda*sum(abs(mu.array)))
}

#dim should be a vector
tensor.index = function(index,dims){
  index = index-1
  Cs = (index %% dims[1])+1
  Ds = (index %% (dims[1]*dims[2]))%/%dims[1] +1
  Es = (index %/% (dims[1]*dims[2]))+1
  return(c(Cs,Ds,Es))
}


Soft = function (a, b){
  if (b < 0) 
    stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a) * pmax(0, abs(a) - b))
}


#give x as an array
UpdateMus.tensor = function (x, Cs, Ds, Es, lambda=0) {
  uniqCs = sort(unique(Cs))
  uniqDs = sort(unique(Ds))
  uniqEs = sort(unique(Es))
  mus = array(NA, c(length(uniqCs), length(uniqDs), length(uniqEs)))
  for (k in uniqCs){
    for (r in uniqDs){
      for (l in uniqEs){
        if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
        if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l)))
        if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
      }
    }
  }
  return(mus)
}

UpdateClusters.tensor = function (x, mus, curCs, curDs) {
  #x = tensor.unfold(x,1); mus=tensor.unfold(mu.array,1); curCs=Cs; curDs=(rep(Ds,each=q)-1)*l+rep(Es,times=p)
  #x = tensor.unfold(x,2); mus=tensor.unfold(mu.array,2); curCs=Ds; curDs=(rep(Es,each=n)-1)*k+rep(Cs,times=q)
  #x = tensor.unfold(x,3); mus=tensor.unfold(mu.array,3); curCs=Es; curDs=(rep(Ds,each=n)-1)*k+rep(Cs,times=p)
  Cs.new <- rep(NA, length(curCs))
  #uniq is useless
  uniq <- 1:max(curCs)
  #curDs=(D1-1)*l+E1;curCs=C1;mus=tensor.unfold(mu.array);x=tensor.unfold(x)
  mus.expandcolumns <- mus[, curDs, drop = FALSE]
  if (dim(mus.expandcolumns)[1]==1) mus
  for (i in 1:nrow(x)) {
    dist2.clust <- NULL
    for (k in 1:length(uniq)) {
      #k=1;i=2
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


classify2 = function(x,k,r,l,lambda=1e-3,max.iter=30,threshold = 5e-3,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL){
  #x=test;k=4;r=3;l=2;step.max=500;threshold=5e-3;max.iter=40;Cs.init=NULL;Ds.init=NULL;Es.init=NULL;lambda=10000
  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
  if(is.null(Cs.init)){
    Cs  = kmeans(tensor.unfold(x,1),k)$cluster
  } else {
    Cs = Cs.init
  }
  if(is.null(Ds.init)){
    Ds  = kmeans(tensor.unfold(x,2),r)$cluster
  } else {
    Ds = Ds.init
  }
  if(is.null(Es.init)){
    Es  = kmeans(tensor.unfold(x,3),l)$cluster
  } else {
    Es = Es.init
  }
  design.sort = cbind(matrix(rep(1:n,each=p*q),nrow=n*p*q),
                      matrix(rep(rep(1:p,each=q),times=n),nrow=n*p*q),
                      matrix(rep(rep(1:q,times=p),times=n),ncol=1),
                      matrix(rep(0,n*p*q*k*r*l),nrow=n*p*q))
  objs <- 1e+15
  improvement <- 1e+10
  i <- 1
  mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda)
  while((improvement > threshold) & (i <= max.iter)){
    #print(1)
    design = t(apply(design.sort,MARGIN=1,FUN=design.row,Cs,Ds,Es,k,r,l))
    #print(1)
    #calculate the mu for each group
    #design.mu = array(apply(design.sort,MARGIN=1,FUN=function(sp,mu)return(mu[Cs[sp[1]],Ds[sp[2]],Es[sp[3]]]),mu=mu.array),c(n,p,q))
    #hold mus and change assignment of row clusters
    Cs = UpdateClusters.tensor(tensor.unfold(x),tensor.unfold(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda))
    Cs <- ReNumber(Cs)
    k = length(unique(Cs))
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda))
    Ds = UpdateClusters.tensor(tensor.unfold(x,2),tensor.unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda))
    Es = UpdateClusters.tensor(tensor.unfold(x,3),tensor.unfold(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda))
    Es <- ReNumber(Es)
    l = length(unique(Es))
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda))
    improvement <- abs(objs[length(objs)] - objs[length(objs) - 
                                                   6])/abs(objs[length(objs) - 6])
    i <- i + 1
    if(trace) cat("step",i,",improvement=",improvement,".\n")
    if (is.na(improvement)) break
  }
  if (i > max.iter) {
    warning("The algorithm has not converged by the specified maximum number of iteration.\n")
  }
  return(list("judgeX"=mu.array[Cs,Ds,Es],"Cs"=Cs,"Ds"=Ds,"Es"=Es,"objs"=objs[length(objs)], "mus"=mu.array))
}


label2 = function(x,k,r,l,lambda=1e-3,max.iter=30,threshold = 5e-3,sim.times=5,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL){
  #x=test;lambda=1e-3;max.iter=200;threshold = 5e-3;sim.times=10
  result = list(); objs = c()
  #objs_sparse = c(); result_sparse = list()
  #for (iteration in 1:sim.times){
  #  if (trace) cat("Iteration:", iteration, ".\n")
  #  result_sparse[[iteration]] = classify2(x,k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init)
  #  objs_sparse = c(objs_sparse, result_sparse[[iteration]]$objs)
  #}
  #result_sparse = result_sparse[[which(objs_sparse == min(objs_sparse))[1]]]
  
  for (iteration in 1:sim.times){
    if (trace) cat("Iteration:", iteration, ".\n")
    result[[iteration]] = classify2(x,k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init)
    objs = c(objs, result[[iteration]]$objs)
  }
  result = result[[which(objs == min(objs))[1]]]
  return(result)
}


simulation  = function(n,p,q,k,r,l,error,lambda,iteration=1){
  cer = c()
  for (i in 1:iteration){
    data = get.data(n,p,q,k,r,l,error=error)
    test = data$x
    truthCs = data$truthCs
    truthDs = data$truthDs
    truthEs = data$truthEs
    sim = label2(test,k,r,l,threshold=5e-2,lambda=0,sim.times=5,trace=FALSE)
    cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
    cer = c(cer,cerC,cerD,cerE)
  }
  cer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,mean)
  print(cer)
}


label_for_circle = function(abc,k,r,l,Cs.init,Ds.init,Es.init,sim.time,lambda,xmiss){
  #label_for_circle(abc=c(1,1,1),k,r,l,Cs.init,Ds.init,Es.init,20,lambda=0,xmiss)
  a = abc[1]; b = abc[2]; c = abc[3]
  return(list(label2(xmiss, k[a], r[b], l[c], lambda = lambda, 
         Cs.init = Cs.init[,a], Ds.init = Ds.init[,b], 
         Es.init = Es.init[,c], sim.time=sim.time)$judgeX))
}
  


#all k,r,l are vectors which is the selection range of k,r,l
sparse.choosekrl = function (x,k,r,l,lambda=0,percent=0.2,trace=FALSE) {
  #x=test;l=range.l;lambda=0;percent=0.2;trace=TRUE
  #k=2:4;r=2:4;l=2:4
  if ((1%%percent) != 0) 
    stop("1 must be divisible by the specified percentage")
  if (percent <= 0) 
    stop("percentage cannot be less than or equal to 0")
  if (percent >= 1) 
    stop("percentage cannot be larger or equal to 1")
  
  #Returns suitably lagged and iterated differences.
  if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0 || sum(diff(l) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
  n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
  numberoftimes <- 1/percent
  allresults <- array(NA, dim = c(numberoftimes, length(k), 
                                  length(r), length(l)))
  Cs.init <- matrix(NA, nrow = dim(x)[1], ncol = length(k))
  #put the kmeans results into columns
  for (i in 1:length(k)) {
    Cs.init[, i] <- kmeans(tensor.unfold(x), k[i], nstart = 20)$cluster
  }
  
  Ds.init <- matrix(NA, nrow = dim(x)[2], ncol = length(r))
  for (j in 1:length(r)) {
    Ds.init[, j] <- kmeans(tensor.unfold(x,2), r[j], nstart = 20)$cluster
  }
  
  Es.init <- matrix(NA, nrow = dim(x)[3], ncol = length(l))
  for (j in 1:length(l)) {
    Es.init[, j] <- kmeans(tensor.unfold(x,3), l[j], nstart = 20)$cluster
  }
  
  for (i in 1:numberoftimes) {#i = 1 
    if (trace == TRUE) 
      cat("Iteration", i, fill = TRUE)
    xmiss <- x
    missing <- sample(1:(n*p*q), round(n*p*q*percent), replace = FALSE)
    xmiss[missing] <- NA
    xmiss[missing] <- mean(xmiss, na.rm = TRUE)
    
    #for (a in 1:length(k)) {
    #  for (b in 1:length(r)) {
    #    for (c in 1:length(l)){#a=1;b=2;c=3
    #      res <- label2(xmiss, k[a], r[b], l[c], lambda = lambda, 
    #                    Cs.init = Cs.init[, a], Ds.init = Ds.init[, b], 
    #                    Es.init = Es.init[,c], sim.time=20)$judgeX
    #      allresults[i, a, b, c] <- sum((x[missing] - res[missing])^2)
    #    }
    #  }
    #}
    
    
    krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
                   rep(1:length(r),times=length(k)*length(l)),
                       rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
                 nrow=3)
    res = apply(krl,MARGIN=2,label_for_circle,k,r,l,Cs.init,Ds.init,Es.init,sim.time=5,lambda=lambda,xmiss=xmiss)
    for (a in 1:dim(krl)[2]){
      #print(krl[,a])
      #print(sum((x[missing]-res[[a]][[1]][missing])^2))
      allresults[i,krl[,a][1],krl[,a][2],krl[,a][3]] = sum((x[missing]-res[[a]][[1]][missing])^2)
    }
    
  }
  #allresults = use.allresults
  #allresults = array((rnorm(5*3*3*3))^2,dim=c(5,3,3,3));k=2:4;r=2:4;l=2:4;numberoftimes=5
  #calculate the standard deviation and mean for different selection of k,r,l
  results.se <- apply(allresults, MARGIN=c(2, 3, 4), sd)
  results.mean <- apply(allresults, c(2, 3, 4), mean)
  #comparing every mean with the mean with higher k,r,l
  IndicatorArray <- 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
                                                              1), 1:(length(l) - 1)] <= results.mean[2:length(k), 2:length(r), 2:length(l)] + results.se[2:length(k), 2:length(r), 2:length(l)])
  if (max(IndicatorArray) == 0) 
    return(list(estimated_krl = c(max(k), max(r), max(l))))
  
  #outer(outer(matrix(1:4,2),matrix(5:8,2)),matrix(1:4,2))
  ModeIndex <- outer(outer(k[-length(k)], r[-length(r)],  "*"), l[-length(l)], "*")
  smallestIndicatorTrue <- min(ModeIndex[IndicatorArray == TRUE])
  out <- which(IndicatorArray == TRUE & ModeIndex == 
                 smallestIndicatorTrue, arr.ind = TRUE)
  out <- array(c(k[out[,1]], r[out[,2]], l[out[,3]]), dim=dim(out))
  tempmode1 <- NULL
  tempmode2 <- NULL
  tempmode3 <- NULL
  for (i in 1:length(k)) {
    tempmode1 <- c(tempmode1, paste("K = ", k[i], sep = ""))
  }
  for (i in 1:length(r)) {
    tempmode2 <- c(tempmode2, paste("R = ", r[i], sep = ""))
  }
  for (i in 1:length(l)) {
    tempmode3 <- c(tempmode3, paste("L = ", l[i], sep = ""))
  }
  
  
  dimnames(results.se)[[1]] = tempmode1
  dimnames(results.se)[[2]] = tempmode2
  dimnames(results.se)[[3]] = tempmode3
  dimnames(results.mean)[[1]] = tempmode1
  dimnames(results.mean)[[2]] = tempmode2
  dimnames(results.mean)[[3]] = tempmode3
  return(list(estimated_krl = out, results.se = results.se, 
              results.mean = results.mean))
}







sim.choosekrl <- function(n,p,q,k,r,l,error=3){
  #n=20;p=40;q=30;k=2;r=4;l=3
  classification<-list()
  for(a in 1:5){
    cat("Starting", a, fill=TRUE)
    #k=3;r=2;l=4;n=30;p=30;q=50
    x = get.data(n,p,q,k,r,l,error=error)$x
    classification[[a]]<-sparse.choosekrl(x,k=2:5,r=2:5,l=2:5)$estimated_kr#estimate clusters
    #if(is.null(classification[[a]])) 
  }
  return(classification)
}




Calculate<-function(true,results){
  real<-matrix(true,ncol=3)
  percent<-0
  for(i in 1:length(results)){
    if(nrow(results[[i]])>1){
      for(a in 1:nrow(results[[i]])){
        #  	  cat("iteration",a,fill=TRUE)
        if(sum(results[[i]][a,]==real)==3){
          percent<-percent+(1/nrow(results[[i]]))
        }
      }
    }
    else if(nrow(results[[i]])<2){
      if(sum(results[[i]]==real)==3){
        percent<-percent+1
      }	
    }
  }
  return(percent/length(results))
}



Calculatekrl<-function(results){
  k<-rep(NA,length(results))
  r<-rep(NA,length(results))
  l<-rep(NA,length(results))
  #length(result): the number of samples
  for(i in 1:length(results)){
    #i = 1
    if(nrow(results[[i]]>1)){
      tempk<-0
      tempr<-0
      templ<-0
      for(a in 1:nrow(results[[i]])){
        #because the return value of Do function sometimes there are not only one classification result in one iteration
        tempk<-tempk+(results[[i]][a,1]/nrow(results[[i]]))
        tempr<-tempr+(results[[i]][a,2]/nrow(results[[i]]))
        templ<-templ+(results[[i]][a,3]/nrow(results[[i]]))
      }
      k[i]<-tempk
      r[i]<-tempr
      l[i]<-templ
    } else if(nrow(results[[i]]<2)){
      k[i]<-results[[i]][1,1]
      r[i]<-results[[i]][1,2]
      l[i]<-results[[i]][1,3]
    }
  }
  return(list(meank=mean(k),meanr=mean(r),meanl=mean(l),sdek=sd(k)/sqrt(length(k)),sder=sd(r)/sqrt(length(r)),sdel=sd(l)/sqrt(length(l)) ))
  #return(list(averagek=k/length(results),averager=r/(length(results))))
}



#clusterobj is the return value of label2()
CalculateBIC = function (x, clusterobj,trace=FALSE) 
{
  mat <- array(rep(0,length(c(x))), dim=dim(x))
  Cs <- clusterobj$Cs
  Ds <- clusterobj$Ds
  Es <- clusterobj$Es
  for (i in unique(Cs)) {
    for (j in unique(Ds)) {
      for (m in unique(Es)) {
        if (clusterobj$mus[i, j, m] != 0) {
          mat[Cs==i, Ds==j, Es==m] <- mean(x[Cs==i, Ds==j, Es==m])
        }
      }
    }
  }
  mat[clusterobj$mus == 0] <- mean(x[clusterobj$mus == 0])
  if(trace==TRUE) cat(log(sum((x - mat)^2)), sum(clusterobj$mus != 0), "\n")
  return(log(sum((x - mat)^2))*dim(x)[1]*dim(x)[2]*dim(x)[3]+log(dim(x)[1]* 
                  dim(x)[2]*dim(x)[3]) * (sum(clusterobj$mus != 0)+1))
}





#lambda: A range of values of tuning parameters to be considered. All values must be non-negative.
chooseLambda = function (x, k, r, l, lambda=NULL) {
  if (is.null(lambda)){
    lambda_0 =c(floor(dim(x)[1]*dim(x)[2]*dim(x)[3]/k/r/l))
    lambda = c(floor(lambda_0/50), floor(lambda_0/10), floor(lambda_0/4), floor(lambda_0/3), floor(lambda_0/2), floor(lambda_0/3*2), lambda_0, floor(lambda_0*3/2), lambda_0*2)
  } 
  x <- x - mean(x)
  BIC <- rep(NA, length(lambda))
  nonzero <- rep(NA, length(lambda))
  for (i in 1:length(lambda)) {
    bires <- label2(x, k, r, l, lambda = lambda[i])
    BIC[i] <- CalculateBIC(x, bires)
    nonzero[i] <- sum(bires$mus != 0)
  }
  return(list(lambda = lambda[which(BIC == min(BIC))[1]], BIC = BIC, 
              nonzeromus = nonzero))
}



base_chooseLambda2 = function(lambda,x,k,r,l){
  bires <- label2(x, k, r, l, lambda)
  return(CalculateBIC(x, bires))
}


chooseLambda2 = function (x, k, r, l) {
  #x = test; k = 4; r = 3; l = 2; lambda=1000; base_chooseLambda2(lambda,x,r,k,l)
  result = optimize(f=base_chooseLambda2, interval=c(0,dim(x)[1]*dim(x)[2]*dim(x)[3]/k/r/l), x=x,k=k,r=r,l=l)
  return(result)
}

label3 = function(x,k,r,l,lambda=1e-3,max.iter=30,threshold = 5e-3,sim.times=5,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL){
  result = label2(x,k,r,l,0,max.iter,threshold,sim.times,trace,Cs.init,Ds.init,Es.init)
  result$judgeX[result$judgeX<=lambda & result$judgeX>=-lambda] = 0
  result$mus[result$mus<=lambda & result$mus>=-lambda] = 0
  return(result)
}

base_chooseLambda3 = function(lambda,x,k,r,l){
  bires <- label3(x, k, r, l, lambda)
  return(CalculateBIC(x, bires))
}

chooseLambda3 = function (x, k, r, l) {
  #x = test; k = 4; r = 3; l = 2; lambda=1000; base_chooseLambda2(lambda,x,r,k,l)
  result = optimize(f=base_chooseLambda3, interval=c(0,1), x=x,k=k,r=r,l=l)
  return(result)
}


error_lambda = function(x,truth){
  return(sum((x-truth)^2)/sum(truth^2))
}


chooseLambda4 = function (x, k, r, l, lambda=NULL,trace=TRUE) {
  if (is.null(lambda)){
    lambda_0 =c(floor(dim(x)[1]*dim(x)[2]*dim(x)[3]/k/r/l)/50)
    lambda = c(0, floor(lambda_0/50), floor(lambda_0/10), floor(lambda_0/4), floor(lambda_0/3), floor(lambda_0/2), floor(lambda_0/3*2), lambda_0, floor(lambda_0*3/2), lambda_0*2)
  } 
  x <- x - mean(x)
  BIC <- rep(NA, length(lambda))
  nonzero <- rep(NA, length(lambda))
  for (i in 1:length(lambda)) {
    bires <- label3(x, k, r, l, lambda = lambda[i])
    BIC[i] <- CalculateBIC(x, bires,trace=trace)
    nonzero[i] <- sum(bires$mus != 0)
  }
  return(list(lambda = lambda[which(BIC == min(BIC))[1]], BIC = BIC, 
              nonzeromus = nonzero))
}


sim.lambda_lasso = function(lambda,n,p,q,k,r,l,sparse.percent,sim.times=20){
  error = c()
  for (i in 1:sim.times){
    data = get.data(n,p,q,k,r,l,sparse.percent = sparse.percent)
    x = data$x
    truth = data$truthX
    x = label3(x,k,r,l,lambda=lambda)$judgeX
    correct_zeros = c(correct_zeros, sum(c(x==0) & c(truth == 0))/sum(c(truth == 0)))
    wrong_zeros = c(wrong_zeros, sum(c(x==0) & c(truth != 0))/sum(truth!=0))
    error = c(error, error_lambda(x,truth))
  }
  return(list(error = c(mean(error),sd(error)), correct_zeros = c(mean(correct_zeros),sd(correct_zeros)), wrong_zeros = c(mean(wrong_zeros), sd(wrong_zeros))))
}



