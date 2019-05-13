#' calculate the correct rate of the result of sparse.choosekrl()
#' 
#' @param true true k, r, l
#' @param results result returned by sparse.choosekrl()
#' 
#' 
Calculate<-function(true,results){
    real<-matrix(true,ncol=3)
    percent<-0
    for(i in 1:length(results)){
        if(nrow(results[[i]])>1){
            for(a in 1:nrow(results[[i]])){
                #        cat("iteration",a,fill=TRUE)
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


#' summary the result of sparse.choosekrl
#' 
#' summary the result of sparse.choosekrl
#' @param results result returned by sparse.choosekrl()
#' 
#' @return \code{meank}
#'         \code{meanr}
#'         \code{meanl}
#'         \code{sdk}
#'         \code{sdr}
#'         \code{sdl}
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
}
#' calculate clustering error rate
#' 
#' calculate clustering error rate
#' @param bires the return list of label2()
#' @param data the return list of get.data()
#' 
#' @return the cer in all modes.
#' 
cer = function(bires,data){
    cerC<-1-adjustedRand(data$truthCs,bires$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(data$truthDs,bires$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(data$truthEs,bires$Es,randMethod=c("Rand"))
    cat("The CER(clustering error rate) is ",cerC,",",cerD,",",cerE,".\n")
    return(list(cerC = cerC, cerD = cerD, cerE = cerE))
}

#' choose the best k,r,l to do clustering when the ranges are given
#' 
#' choose the best k,r,l to do clustering when the ranges are given by minimizing the BIC.
#' @param x a three-dimensional array
#' @param k a vector, the possible clusters number of mode 1
#' @param r a vector, the possible clusters number of mode 2
#' @param l a vector, the possible clusters number of mode 3
#' @param lambda a numeric value
#' @param sim.times the same as label2()
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @return a list \code{estimated_krl} ...\\
#'                \code{BIC} ...
#' 
choosekrl_bic = function (x,k,r,l,lambda=0,sim.times=1,method="L0"){
    ## x = x - mean(x) ## commented out by Miaoyan
    if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0 || sum(diff(l) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
    n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
    krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
    rep(1:length(r),times=length(k)*length(l)),
    rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
    nrow=3)
    krl_list = as.list(as.data.frame(krl))
    if (.Platform$OS.type == "windows") {
        bires = apply(krl,MARGIN=2,label_for_krl,k,r,l,sim.times=sim.times,lambda=lambda,xmiss=x,method=method,crossvalidation=FALSE)
        CBIC = lapply(bires,tensor.calculateBIC, x=x,method=method)
        BIC = unlist(CBIC)
    } else {
        bires = mclapply(krl_list, label_for_krl,k,r,l,sim.times=sim.times,lambda=lambda,xmiss=x,method=method,crossvalidation=FALSE,mc.cores=n.cores)
        CBIC = mclapply(bires,tensor.calculateBIC, x=x,method=method,mc.cores = n.cores)
        BIC = unlist(CBIC)
    }
    names(BIC) = apply(krl,MARGIN=2,FUN=function(x)paste(k[x[1]],r[x[2]],l[x[3]]))
    best = krl_list[[which(BIC == min(BIC))[1]]]
    return(list(estimated_krl = t(as.matrix(c(k[best[1]],r[best[2]],l[best[3]]))), BIC = BIC))
}

#' choose the best lambda when the range is given
#' 
#' choose the best lambda for a certain range.
#' @param x a three-dimensional array
#' @param k the clusters number of mode 1
#' @param r the clusters number of mode 2
#' @param l the clusters number of mode 3
#' @param lambda a vector of possible lambda, for example: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @return a list \code{lambda} the lambda with lowest BIC;
#'                \code{BIC} the corresponding BIC for each lambda in the given range;
#'                \code{nonzeromus} the number clusters with non-zero mean.
#' 
chooseLambda = function (x, k, r, l, lambda=NULL,method="L0") {
    ##x = x - mean(x) commented out by Miaoyan
    if (is.null(lambda)){
        n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
        if (method == "L0") lambda = (n*p*q)/(k*r*l)*seq(0,0.8,by=0.04)
        if (method == "L1") lambda = (n*p*q)/(k*r*l)*seq(0,1,by=0.05)
        if (is.null(lambda)) stop("No such kind of method:", method, ".\n")
    } 
    if (.Platform$OS.type == "windows") {
        bires = lapply(lambda,FUN=classify2,x=x,k=k,r=r,l=l,method=method)
        CBIC = lapply(bires,tensor.calculateBIC, x=x,method=method)
        BIC = unlist(CBIC)
        nonzero = unlist(lapply(bires, FUN=function(bires){return(sum(bires$mus!=0))}))
    } else {
        bires = mclapply(lambda,FUN=classify2,x=x,k=k,r=r,l=l,method=method,mc.cores = n.cores)
        CBIC = mclapply(bires,tensor.calculateBIC, x=x,method=method,mc.cores = n.cores)
        BIC = unlist(CBIC)
        nonzero = unlist(mclapply(bires, FUN=function(bires){return(sum(bires$mus!=0))},mc.cores = n.cores))
    }
    return(list(lambda = lambda[which(BIC == min(BIC))[1]], BIC = BIC, 
    nonzeromus = nonzero))
}

#' perform tensor clustering
#' 
#' perform tensor clustering
#' @param x a three-dimensional array
#' @param k the number of clusters in mode 1
#' @param r the number of clusters in mode 2
#' @param l the number of clusters in mode 3
#' @param lambda the sparsity parameter in the sparse clustering
#' @param max.iter maximum times of iteration
#' @param threshold ...
#' @param trace ...
#' @param Cs.init ...
#' @param Ds.init ...
#' @param Es.init ...
#' @param nstart ...
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicates Lasso penalty, and "L0" indicates sparse subset penalty
#' @param center ...
#' @return a list \code{judgeX}
#'                \code{Cs} clustering result in mode 1 
#'                \code{Ds} clustering result in mode 2
#'                \code{Es} clustering result in mode 3
#'                \code{mus}
#' 
#' 
classify2 = function(x,k,r,l,lambda=0,max.iter=1000,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,nstart=20,method="L0",center=FALSE){
    n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
    if(center == TRUE) {
        mustemp <- mean(x)
        x <- x-mustemp
    }
    if(is.null(Cs.init)){
        if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor.unfold(x,1),k,nstart = nstart)$cluster}
    } else {
        Cs = Cs.init
    }
    if(is.null(Ds.init)){
        if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor.unfold(x,2),r,nstart = nstart)$cluster}
    } else {
        Ds = Ds.init
    }
    if(is.null(Es.init)){
        if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor.unfold(x,3),l,nstart = nstart)$cluster}
    } else {
        Es = Es.init
    }
    objs <- 1e+15
    improvement <- 1e+10
    i <- 1
    mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
    while((improvement > threshold) & (i <= max.iter)){
        Cs = UpdateClusters.tensor(tensor.unfold(x),tensor.unfold(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
        Cs <- ReNumber(Cs)
        k = length(unique(Cs))
        ##mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
        Ds = UpdateClusters.tensor(tensor.unfold(x,2),tensor.unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
        Ds <- ReNumber(Ds)
        r = length(unique(Ds))
        ##mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda,method=method)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
        Es = UpdateClusters.tensor(tensor.unfold(x,3),tensor.unfold(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
        Es <- ReNumber(Es)
        l = length(unique(Es))
        mu.array = UpdateMus.tensor(x,Cs,Ds,Es,lambda, method=method)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
        improvement <- abs(objs[length(objs)] - objs[length(objs) - 
        6])/abs(objs[length(objs) - 6])
        i <- i + 1
        if(trace) cat("step",i,",improvement=",improvement,".\n")
        if (is.na(improvement)) break
    }
    if (i > max.iter) {
        warning("The algorithm has not converged by the specified maximum number of iteration.\n")
    }
    if(center==TRUE){
        mu.array <- mu.array + mustemp
    }
    mu.array[abs(mu.array)<=1e-6] = 0
    return(list("judgeX"=mu.array[Cs,Ds,Es, drop=FALSE],"Cs"=Cs,"Ds"=Ds,"Es"=Es,"objs"=objs[length(objs)], "mus"=mu.array))
}

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
                    if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l)),method=method)
                    if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
                }
            }
        }
    }## added by Miaoyan  
    if(method=="L0"){
        for (k in uniqCs){
            for (r in uniqDs){
                for (l in uniqEs){
                    if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
                    if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/sqrt(sum(Cs==k)*sum(Ds==r)*sum(Es==l)),method=method)
                    ### modified by Miaoyan
                    ## Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/sum(Cs==k)*sum(Ds==r)*sum(Es==l),method=method)
                    if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
                }
            }
        }
    }## added by Miaoyan  
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

#' generate a random tensor
#' 
#' generate a random tensor
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k the clusters number of mode 1 
#' @param r the clusters number of mode 2
#' @param l the clusters number of mode 3 
#' @param error noise in rnorm
#' @param sort if TRUE, the data belongs to the same clusters would assemble together
#' @param sparse.percent the proportion of 0 mean
#' @param center if center sum(x) = 0
#' @param mumin ...
#' @param mumax ...
#' @return a list \code{x}
#'                \code{truthX}
#'                \code{truthCs} true distribution in mode 1
#'                \code{truthDs} true distribution in mode 2
#'                \code{truthEs} true distribution in mode 3
#'                \code{mus} 
#'                \code{binaryX}
#' 
get.data = function(n,p,q,k,r,l,error=3,sort=TRUE,sparse.percent=0,center=FALSE,mumin = -3, mumax = 3){
    mus = runif(k*r*l,mumin,mumax)#take the mean of k*r*l biclusters/cubes
    if(sparse.percent!=0) mus[1:floor(k*r*l*sparse.percent)] = 0
    mus = array(mus,c(k,r,l))
    if (k!=1) truthCs = ReNumber(sample(1:k,n,replace=TRUE)) else truthCs = rep(1,n)
    if (r!=1) truthDs = ReNumber(sample(1:r,p,replace=TRUE)) else truthDs = rep(1,p)
    if (l!=1) truthEs = ReNumber(sample(1:l,q,replace=TRUE)) else truthEs = rep(1,q)
    
    ##### added by Miaoyan
    if(sort==TRUE){
        truthCs=sort(truthCs)
        truthDs=sort(truthDs)
        truthEs=sort(truthEs)
    }
    ######
    
    x = array(rnorm(n*p*q,mean=0,sd=error),dim = c(n,p,q))
    truthX = array(rep(0,n*p*q),c(n,p,q))
    for(i in 1:max(truthCs)){
        for(j in 1:max(truthDs)){
            for(m in 1:max(truthEs)){
                x[truthCs==i, truthDs==j, truthEs==m] = x[truthCs==i, truthDs==j, truthEs==m] + mus[i,j,m]
                truthX[truthCs==i, truthDs==j, truthEs==m] =  mus[i,j,m]
            }
        }
    }
    if (center == TRUE) x = x - mean(x)
    
    #### removed by Miaoyan
    #if(sort==TRUE){
    # truthX = reorderClusters(truthX,truthCs,truthDs,truthEs)$x
    # neworder = reorderClusters(x,truthCs,truthDs,truthEs)
    # binaryX = (truthX!=0)*1
    # result = list("x"=neworder$x,"truthX"=truthX,"truthCs"=neworder$Cs,"truthDs"=neworder$Ds,"truthEs"=neworder$Es,"mus"=mus,"binaryX"=binaryX)
    #} else {
    binaryX = (truthX!=0)*1
    result = list("x"=x,"truthX"=truthX,"truthCs"=truthCs,"truthDs"=truthDs,"truthEs"=truthEs,"mus"=mus,"binaryX"=binaryX)
    #}
    
    return(result)
}

#' perform tensor clustering
#' 
#' perform tensor clustering in a more stable way
#' @param x a three-dimensional array
#' @param k the number of clusters in mode 1
#' @param r the number of clusters in mode 2
#' @param l the number of clusters in mode 3
#' @param lambda the sparsity parameter in the sparse clustering
#' @param max.iter maximum times of iteration
#' @param threshold ...
#' @param sim.times the times of calling classify2() with different seeds.
#' @param trace ...
#' @param Cs.init ...
#' @param Ds.init ...
#' @param Es.init ...
#' @param nstart ...
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicates Lasso penalty, and "L0" indicates sparse subset penalty
#' @param center ...
#' @return a list \code{judgeX}
#'                \code{Cs} clustering result in mode 1 
#'                \code{Ds} clustering result in mode 2
#'                \code{Es} clustering result in mode 3
#'                \code{mus}
#' 
#' 
label2 = function(x,k,r,l,lambda=0,max.iter=1000,threshold = 1e-10,sim.times=1,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,method="L0"){
    #x=test;lambda=1e-3;max.iter=200;threshold = 5e-3;sim.times=10
    if (sim.times == 1) return(classify2(x,k,r,l,lambda=lambda,max.iter = max.iter,threshold = threshold,Cs.init = Cs.init,Ds.init = Ds.init,Es.init = Es.init,method=method))
    if (.Platform$OS.type == "windows") {
        result = lapply(rep(list(x),sim.times), classify2, k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,method=method)
        objs = unlist(lapply(result, function(result){result$objs}))
    } else {
        result = mclapply(rep(list(x),sim.times), classify2, k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,nstart = sample(1:1000,1),method=method,mc.cores = n.cores)
        objs = unlist(lapply(result, function(result){result$objs}))
    }
    result = result[[which(objs == min(objs))[1]]]
    return(result)
}

#' draw 3d plot of tensor
#' 
#' draw 3d plot of tensor
#' @param tensor a three-dimensional array
#' 
plot_tensor=function(tensor){
    
    n=prod(dim(tensor))   
    color_choice=min(round(prod(dim(tensor))/(6*6*6)),100)+1
    marker=viridis_pal(option = "B")(color_choice)
    
    position=positionfun(dim(tensor))$position
    quan=c(quantile(tensor,(0:color_choice)/color_choice))
    col=tensor
    for(i in 1:color_choice){
        col[(tensor>=quan[i])&(tensor<quan[i+1])]=marker[i]
    }
    col[tensor==quan[i+1]]=marker[i]
    
    plot3d(position[,1],position[,2],position[,3],col=col,alpha=0.3,size=5,xlab="",ylab="",zlab="")
}

#' simulation: choose the best k,r,l to do clustering when the ranges are given
#' 
#' simulation: choose the best k,r,l to do clustering when the ranges are given
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k true clusters number of mode 1
#' @param r true clusters number of mode 2
#' @param l true clusters number of mode 3
#' @param error noise when producing the data
#' @param sim.times simulation times
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @param mode two options: "bic", "crossvalidation".
#' @return A list of the krl result in different iteration.
#' 


sim.choosekrl <- function(n,p,q,k,r,l,error=1,sim.times=5,method="L0",mode="bic"){
    #n=20;p=40;q=30;k=2;r=4;l=3
    classification<-list()
    if(mode == "bic"){
        for(a in 1:sim.times){
            cat("Starting", a, fill=TRUE)
            x = get.data(n,p,q,k,r,l,error=error)$x
            classification[[a]]<-choosekrl_bic(x,k=2:5,r=2:5,l=2:5,method=method)$estimated_krl#estimate clusters
        }
    }
    if(mode == "crossvalidation"){
        for(a in 1:sim.times){
            cat("Starting", a, fill=TRUE)
            x = get.data(n,p,q,k,r,l,error=error)$x
            classification[[a]]<-sparse.choosekrl(x,k=2:5,r=2:5,l=2:5,method=method)$estimated_krl#estimate clusters
        }
    }
    if(mode!="bic" & mode!="crossvalidation") stop("No such kind of mode:", mode, ".\n")
    return(classification)
}

#' simulation: choose the best lambda when the range is given
#' 
#' choose the best lambda for a certain range.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k the clusters number of mode 1
#' @param r the clusters number of mode 2
#' @param l the clusters number of mode 3
#' @param iteration iteration times
#' @param lambda a vector of possible lambda, for example: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)
#' @param standarddeviation the standard deviation when producing data
#' @param center if True, then x = x- mean(x)
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @return a list \code{selectedlambda} a vector of lambdas with lowest BIC in each iteration
#' 
sim.chooseLambda = function(n,p,q,k,r,l,iteration,lambda,standarddeviation=4,center = FALSE,method="L0"){
    selectedLambda = 1:iteration
    for(iter in 1:iteration){
        cat("Iteration",iter,fill=TRUE)
        set.seed(iter)
        smp = get.data(n,p,q,k,r,l,error=standarddeviation,sort=FALSE)$x
        if(center == TRUE)smp = smp - mean(smp)
        selectedLambda[iter] = chooseLambda(smp,k,r,l,lambda=lambda,method=method)$lambda
    }
    return(selectedLambda)
}

#' simulation: choose the best lambda when the range is given
#' 
#' choose the best lambda for a certain range.
#' @param n the dimension of mode 1
#' @param p the dimension of mode 2
#' @param q the dimension of mode 3
#' @param k the clusters number of mode 1
#' @param r the clusters number of mode 2
#' @param l the clusters number of mode 3
#' @param error the standard deviation when producing data
#' @param lambda a vector of possible lambda, for example: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)
#' @param iteration iteration times
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
simulation  = function(n,p,q,k,r,l,error,lambda,iteration=1,method="L0"){
    cer = c()
    for (i in 1:iteration){
        data = get.data(n,p,q,k,r,l,error=error)
        test = data$x
        truthCs = data$truthCs
        truthDs = data$truthDs
        truthEs = data$truthEs
        sim = label2(test,k,r,l,threshold=5e-2,lambda=0,trace=FALSE,method=method)
        cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
        cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
        cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
        cer = c(cer,cerC,cerD,cerE)
    }
    cer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,mean)
    print(cer)
}

#' choose the best k,r,l to do clustering when the ranges are given
#' 
#' choose the best k,r,l to do clustering when the ranges are given using cross validation.
#' @param x a three-dimensional array
#' @param k a vector, the possible clusters number of mode 1
#' @param r a vector, the possible clusters number of mode 2
#' @param l a vector, the possible clusters number of mode 3
#' @param lambda a numeric value
#' @param percent a numeric value between 0 and 1
#' @param trace ...
#' @param nstart ...
#' @param sim.times the same as label2()
#' @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
#' @return a list \code{estimated_krl} ...\\
#'                \code{results.se} ...\\
#'                \code{results.maen} ...\\
#' 
sparse.choosekrl = function (x,k,r,l,lambda=0,percent=0.2,trace=FALSE,nstart=20,sim.times=1,method="L0") {
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
    miss <- sample(1:(n*p*q), n*p*q, replace = FALSE)
    numberoftimes <- 1/percent
    allresults <- array(NA, dim = c(numberoftimes, length(k), 
    length(r), length(l)))
    Cs.init <- matrix(NA, nrow = dim(x)[1], ncol = length(k))
    #put the kmeans results into columns
    for (i in 1:length(k)) {
        Cs.init[, i] <- kmeans(tensor.unfold(x), k[i], nstart = nstart)$cluster
    }
    
    Ds.init <- matrix(NA, nrow = dim(x)[2], ncol = length(r))
    for (j in 1:length(r)) {
        Ds.init[, j] <- kmeans(tensor.unfold(x,2), r[j], nstart = nstart)$cluster
    }
    
    Es.init <- matrix(NA, nrow = dim(x)[3], ncol = length(l))
    for (j in 1:length(l)) {
        Es.init[, j] <- kmeans(tensor.unfold(x,3), l[j], nstart = nstart)$cluster
    }
    
    for (i in 1:numberoftimes) {#i = 1 
        if (trace == TRUE) 
        cat("Iteration", i, fill = TRUE)
        xmiss <- x
        missing <- miss[1:round(n*p*q*percent)]
        xmiss[missing] <- NA
        xmiss[missing] <- mean(xmiss, na.rm = TRUE)
        
        
        krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
        rep(1:length(r),times=length(k)*length(l)),
        rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
        nrow=3)
        
        
        if (.Platform$OS.type == "windows") {
            res = apply(krl,MARGIN=2,label_for_krl,k,r,l,Cs.init,Ds.init,Es.init,sim.times=sim.times,lambda=lambda,xmiss=xmiss,method=method)
        } else {
            krl_list = as.list(as.data.frame(krl))
            res = mclapply(krl_list, label_for_krl,k,r,l,Cs.init,Ds.init,Es.init,sim.times=sim.times,lambda=lambda,xmiss=xmiss,method=method,mc.cores=n.cores)
        }
        
        
        for (a in 1:dim(krl)[2]){
            #print(krl[,a])
            #print(sum((x[missing]-res[[a]][[1]][missing])^2))
            allresults[i,krl[,a][1],krl[,a][2],krl[,a][3]] = sum((x[missing]-res[[a]][[1]][missing])^2)
        }
        miss <- miss[-1:-(dim(x)[1] * dim(x)[2] * dim(x)[3]/numberoftimes)]
    }
    results.se <- apply(allresults, MARGIN=c(2, 3, 4), sd)/sqrt(numberoftimes)
    results.mean <- apply(allresults, c(2, 3, 4), mean)
    #comparing every mean with the mean with higher k,r,l
    IndicatorArray <- 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
    1), 1:(length(l) - 1)] <= results.mean[2:length(k), 2:length(r), 2:length(l)] + results.se[2:length(k), 2:length(r), 2:length(l)])
    if (max(IndicatorArray) == 0) {
        warning("The krl has reached the upper boundary. Please enlarger the range.")
        return(list(estimated_krl = c(max(k),max(r),max(l))))}
    
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

#' evaluate the accuracy of the model when inputting a sparse tensor
#' 
#' evaluate the accuracy of the model when inputting a sparse tensor
#' @param bires the return list of label2()
#' @param data the return list of get.data()

sparse.evaluate = function(bires, data){
    npq=dim(data$x);n=npq[1];p=npq[2];q=npq[3]
    restensor<-(abs(bires$judgeX)>1e-10)*1
    totalzero<-sum(restensor==0)/(n*p*q)
    correctzero<-sum(restensor[which(data$binaryX==0)]==0)/sum(data$binaryX==0)
    correctone<-sum(restensor[which(data$binaryX==1)]==1)/(n*p*q-sum(data$binaryX==0))
    totalincorrect<-sum(abs(restensor-data$binaryX))/(n*p*q)
    cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect,".\n")
}

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






