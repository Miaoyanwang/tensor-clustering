
### simulation
### This script set simulation and using constrained simulation

gene_data = function(whole_shape = c(20,20,20), core_shape = c(3,3,3)){
  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  ####-------- generate data
  set.seed(24)  # 24 # 37  #  347
  A = randortho(d1)[,1:r1]     ## factor matrix
  B = randortho(d2)[,1:r2]     ## factor matrix
  C = randortho(d3)[,1:r3]     ## factor matrix
  
  ### G: core tensor
  G = as.tensor(array(data = rnorm(r1*r2*r3,mean = 0,sd = 10),dim = core_shape))
  ###
  #G = as.tensor(array(data = runif(r1*r2*r3,0,1),dim = core_shape))
  
  ### U: ground truth
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  
  #hist(U); quantile(U,c(0,0.01,0.25,0.75,0.99,1))
  #hist(inv.logit(U)); quantile(inv.logit(U),c(0,0.01,0.25,0.75,0.99,1))
  
  ### ts:binary tensor
  ts = list()
  for (i in 1:5) {
    binary = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
    ts[[i]] = as.tensor(array(binary,dim = c(d1,d2,d3)))@data
  }
  
  return(list(U = U,ts = ts))
}


whole_shape = c(20,20,20); core_shape = c(3,3,3)
data = gene_data()
U = data$U
ts = data$ts


#####---------------------------  selecting lambda
d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]

lambda = c(0.01,0.1,0.5,1,5, 500 )
re = list()
for(i in 1:6){
  upp2 = update_binary_cons(ts,c(r1,r2,r3),Nsim = 100, lambda = lambda[i], alpha = 10* max(U))
  re[[i]] = upp2
}

MSE = rep(0,6)
for (i  in 1:6) {
  upp2 = re[[1]]
  U_est2 = ttl(upp2$G,list(upp2$A,upp2$B,upp2$C),ms = c(1,2,3))@data
  
  MSE[i] = sum((U_est2 - U)^2)/prod(whole_shape)
}

plot(MSE)


plot(as.vector(U),as.vector(U_est2));abline(0,1)
plot(as.vector(inv.logit(U)),as.vector(inv.logit(U_est2)));abline(0,1)


######------------------------------   convergence rate 
## --------------------   using penalty
##-------  lambda = 1
conv_rate = function(d, r){
  rate = rep(0,length(d))
  RMSE = rep(0,length(d))
  for (i in 1:length(d)) {
    data = gene_data(rep(d[i],3), rep(r[i],3))
    U = data$U
    ts = data$ts
    RMSEi = rep(0,5)
    for (j in 1:5) {
      upp = update_binary_cons(ts[[j]],rep(r[i],3),50, lambda = 1, alpha = 10*max(U))
      U_est = ttl(upp$G,list(upp$A,upp$B,upp$C),ms = c(1,2,3))@data
      RMSEi[j] = sqrt(sum((U_est - U)^2)/(d[i]^3))
      print(paste(j,"-th observation ---- when dimension is ",d[i],"-- rank is ",r[i]," ---------"))
    }
    RMSE[i] = mean(RMSEi)
    rate[i] = r[i]^2/d[i]^2
  }
  return(list(RMSE = RMSE, rate = rate))
}

re3 = conv_rate(c(20,50,80),c(3,3,3))
re5 = conv_rate(c(20,50,80),c(5,5,5))

re3$RMSE*sqrt((20^3))/sqrt(c(20^3,50^3,80^3))
re3$rate

plot(re3$RMSE*sqrt((20^3))/sqrt(c(20^3,50^3,80^3)),
re3$rate)

re = data.frame(rate = 1:10, RMSE = rnorm(10))
library(ggplot2)
ggplot(re, aes(x = rate, y = RMSE)) + geom_line(color = 'navyblue') +
  geom_point(color = 'red', size = 4)


## --------------------   using momentum method
##-------
conv_rate = function(d, r){
  rate = rep(0,length(d))
  RMSE = rep(0,length(d))
  for (i in 1:length(d)) {
    data = gene_data(rep(d[i],3), rep(r[i],3))
    U = data$U
    ts = data$ts
    RMSEi = rep(0,5)
    for(j in 1:5){
      upp = update_binary_momentum(ts[[j]],rep(r[i],3),40, alpha = 10*max(U))
      U_est = ttl(upp$G,list(upp$A,upp$B,upp$C),ms = c(1,2,3))@data
      RMSEi[j] = sqrt(sum((U_est - U)^2)/d[i]^3)
      print(paste(j,"-th observation ---- when dimension is ",d[i],"-- rank is ",r[i]," ---------"))
    }
    
    RMSE[i] = mean(RMSEi)
    rate[i] = r[i]^2/d[i]^2
  }
  return(list(RMSE = RMSE, rate = rate))
}

re3 = conv_rate(c(20,50,80),c(3,3,3))
re5 = conv_rate(c(20,50,80),c(5,5,5))

re3

