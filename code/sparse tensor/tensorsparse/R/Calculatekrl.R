#' Summarize the result of estimating the true c(d_1,d_2,d_3) in a simulation.
#' 
#' Summarize the return objct of sparse.choosekrl().
#' @param results List. A list consists of 3-dimensional vectors and each vector is estimated c(d_1,d_2,d_3). (The result returned by sparse.choosekrl().)
#' 
#' @return \code{meank} mean estimated d_1
#'         \code{meanr} mean estimated d_2
#'         \code{meanl} mean estimated d_3
#'         \code{sdk} the standard deviation of estimated d_1
#'         \code{sdr} the standard deviation of estimated d_2
#'         \code{sdl} the standard deviation of estimated d_3
#' @export
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
