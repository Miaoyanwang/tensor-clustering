UpdateMus = function (x, Cs, Ds, lambda = 0) {
  uniqCs <- sort(unique(Cs))
  uniqDs <- sort(unique(Ds))
  mus <- matrix(NA, nrow = length(uniqCs), ncol = length(uniqDs))
  for (k in uniqCs) {
    for (r in uniqDs) {
      if (lambda == 0) 
        mus[k, r] <- mean(x[Cs == k, Ds == r])
      if (lambda > 0) 
        mus[k, r] <- Soft(mean(x[Cs == k, Ds == r]), 
                          lambda/(sum(Cs == k) * sum(Ds == r)))
      if (lambda < 0) 
        stop("Cannot have a negative tuning parameter value.")
    }
  }
  return(mus)
}

Soft = function (a, b){
  if (b < 0) 
    stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a) * pmax(0, abs(a) - b))
}

UpdateClusters = function (x, mus, curCs, curDs) {
  Cs.new <- rep(NA, length(curCs))
  #uniq is useless
  uniq <- 1:max(curCs)
  #curDs=Ds;curCs=Cs
  mus.expandcolumns <- mus[, curDs, drop = FALSE]
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

Objective = function (x, mus, Cs, Ds, lambda = 0) {
  return(sum((x - mus[Cs, Ds])^2) + 2 * lambda * sum(abs(mus)))
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

#here x is the data matrix
sparseBC = function (x, k, r, lambda, nstart = 20, Cs.init = NULL, Ds.init = NULL, 
                     max.iter = 1000, threshold = 1e-10, center = TRUE) {
  #x = matrix(1:250,ncol=25);k=2;r=2;lambda=1;nstart = 20;Cs.init = NULL;Ds.init = NULL;max.iter = 1000;threshold = 1e-10;center = TRUE
  if (is.null(Cs.init)) {
    Cs <- (kmeans(x, k, nstart = nstart)$cluster)
  } else {
    Cs <- Cs.init
  }
  if (is.null(Ds.init)) {
    Ds <- (kmeans(t(x), r, nstart = nstart)$cluster)
  } else {
    Ds <- Ds.init
  }
  if (center == TRUE) {
    mustemp <- mean(x)
    x <- x - mustemp
  }
  #match.call returns a call in which all of the specified arguments are specified by their full names.
  cl <- match.call()
  mus <- UpdateMus(x, Cs, Ds, lambda = lambda)
  objs <- 1e+15
  improvement <- 1e+10
  i <- 1
  while (improvement > (threshold) && i <= max.iter) {
    #update Cs
    
    Cs <- UpdateClusters(x, mus, Cs, Ds)
    objs <- c(objs, Objective(x, mus, Cs, Ds, lambda = lambda))
    Cs <- ReNumber(Cs)
    mus <- UpdateMus(x, Cs, Ds, lambda = lambda)
    objs <- c(objs, Objective(x, mus, Cs, Ds, lambda = lambda))
    Ds <- UpdateClusters(t(x), t(mus), Ds, Cs)
    objs <- c(objs, Objective(x, mus, Cs, Ds, lambda = lambda))
    Ds <- ReNumber(Ds)
    mus <- UpdateMus(x, Cs, Ds, lambda = lambda)
    objs <- c(objs, Objective(x, mus, Cs, Ds, lambda = lambda))
    improvement <- abs(objs[length(objs)] - objs[length(objs) - 
                                                   4])/abs(objs[length(objs) - 4])
    i <- i + 1
  }
  if (i > max.iter) {
    warning("The algorithm has not converged by the specified maximum number of iteration")
  }
  if (center == TRUE) {
    mus <- mus + mustemp
  }
  out <- list()
  class(out) <- "sparseBC"
  out$Cs <- Cs
  out$Ds <- Ds
  out$objs <- objs
  out$mus <- mus[Cs, Ds]
  out$Mus <- mus
  out$iteration <- i
  out$cl <- cl
  return(out)
}

sparseBC.choosekr = function (x, k, r, lambda, percent = 0.1, trace = FALSE) {
  if ((1%%percent) != 0) 
    stop("1 must be divisible by the specified percentage")
  if (percent <= 0) 
    stop("percentage cannot be less than or equal to 0")
  if (percent >= 1) 
    stop("percentage cannot be larger or equal to 1")
  if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
  miss <- sample(1:(nrow(x) * ncol(x)), nrow(x) * ncol(x), 
                 replace = FALSE)
  numberoftimes <- 1/percent
  allresults <- array(NA, dim = c(numberoftimes, length(k), 
                                  length(r)))
  Cs.init <- matrix(NA, nrow = nrow(x), ncol = length(k))
  for (i in 1:length(k)) {
    Cs.init[, i] <- kmeans(x, k[i], nstart = 20)$cluster
  }
  Ds.init <- matrix(NA, nrow = ncol(x), ncol = length(r))
  for (j in 1:length(r)) {
    Ds.init[, j] <- kmeans(t(x), r[j], nstart = 20)$cluster
  }
  for (i in 1:numberoftimes) {
    if (trace == TRUE) 
      cat("Iteration", i, fill = TRUE)
    xmiss <- x
    missing <- miss[1:(nrow(x) * ncol(x))/numberoftimes]
    xmiss[missing] <- NA
    xmiss[missing] <- mean(xmiss, na.rm = TRUE)
    for (a in 1:length(k)) {
      for (b in 1:length(r)) {
        res <- sparseBC(xmiss, k[a], r[b], lambda = lambda, 
                        Cs.init = Cs.init[, a], Ds.init = Ds.init[, 
                                                                  b])$mus
        allresults[i, a, b] <- sum((x[missing] - res[missing])^2)
      }
    }
    miss <- miss[-1:-(dim(x)[1] * dim(x)[2]/numberoftimes)]
  }
  results.se <- apply(allresults, c(2, 3), sd)/sqrt(numberoftimes)
  results.mean <- apply(allresults, c(2, 3), mean)
  IndicatorMatrix <- 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
                                                               1)] <= results.mean[2:length(k), 2:length(r)] + results.se[2:length(k), 
                                                                                                                          2:length(r)])
  if (max(IndicatorMatrix) == 0) 
    return(list(bestK = max(k), bestR = max(r)))
  RowIndexPlusColIndex <- outer(k[-length(k)], r[-length(r)], 
                                "*")
  smallestIndicatorTrue <- min(RowIndexPlusColIndex[IndicatorMatrix == 
                                                      TRUE])
  out <- which(IndicatorMatrix == TRUE & RowIndexPlusColIndex == 
                 smallestIndicatorTrue, arr.ind = TRUE)
  out <- matrix(c(k[out[, 1]], r[out[, 2]]), nrow(out), ncol(out))
  temprow <- NULL
  tempcol <- NULL
  for (i in 1:length(k)) {
    temprow <- c(temprow, paste("K = ", k[i], sep = ""))
  }
  for (i in 1:length(r)) {
    tempcol <- c(tempcol, paste("R = ", r[i], sep = ""))
  }
  rownames(results.se) <- temprow
  colnames(results.se) <- tempcol
  rownames(results.mean) <- temprow
  colnames(results.mean) <- tempcol
  return(list(estimated_kr = out, results.se = results.se, 
              results.mean = results.mean))
}