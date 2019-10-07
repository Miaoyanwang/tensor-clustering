
###  select rank using BIC criterion

gene_data = function(whole_shape = c(30,30,30), core_shape = c(3,3,3)){
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
  ts = rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)) ) )
  ts = as.tensor(array(ts,dim = c(d1,d2,d3)))@data
  return(list(U = U,ts = ts))
}


whole_shape = c(30,30,30); core_shape = c(3,3,3)
data = gene_data()
U = data$U
ts = data$ts


sele_rank = function(rank){
  BIC = rep(0,length(rank))
  for(i in 1:length(rank)){
    upp = update_binary_cons(ts,rep(rank[i],3),100, lambda = 1, alpha = 10*max(U))
    log_Lik = max(upp$lglk)
    BIC[i] = -2log_Lik + (rank[i]^3 + sum(whole_shape*rank[i]))*log(prod(whole_shape))
  }
  return(rank(which(BIC = min(BIC))))
}








