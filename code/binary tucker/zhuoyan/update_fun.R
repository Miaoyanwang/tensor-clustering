#################  update

update_binary = function(ts, Y_covar , core_shape, penalty_par, Nsim){
  ## get initialization
  ts1 = 10*(2*ts - 1)
  ts1 = as.tensor(ts1)
  ts = as.tensor(ts)
  tckr = tucker(ts1, ranks = core_shape)
  A = tckr$U[[1]] ; B = tckr$U[[2]] ; C = tckr$U[[3]]
  G = tckr$Z
  d1 = dim(ts)[1] ; d2 = dim(ts)[2] ; d3 = dim(ts)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3] 
  Y_1 = unfold(ts, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(ts, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(ts, row_idx = 3, col_idx = c(1,2))@data
  
  alpha = penalty_par[1]
  beta = penalty_par[2]
  lambda = penalty_par[3]
  gamma = penalty_par[4]
  
  num_feature = dim(Y_covar)[2]
  W = matrix(rnorm(r1*num_feature),r1,num_feature)
  V = matrix(rnorm(r1*num_feature),r1,num_feature)
  
  cross_entropy = rep(0,4*Nsim)
  for(n in 1:Nsim){
    
    ###### update A
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    
    re = glm_mat_regu(t(A),t(Y_1),t(Y_covar),t(G_BC1),alpha,t(W))
    
    #Y_1_and_Start = cbind(Y_1,A)
    #re = glm_mat(t(Y_1_and_Start),t(G_BC1))
    #re = glm_mat(t(Y_1),start = t(A),t(G_BC1))
    
    A = t(re[[1]])
    cross_entropy[4*n - 3] = re[[2]] + beta*norm(Y_covar - B%*%V,type = 'F')
    ## orthogonal A*
    U = ttl(G,list(A,B,C),ms = c(1,2,3))
    tuk = tucker(U, ranks = core_shape)
    G = tuk$Z
    A = tuk$U[[1]]
    B = tuk$U[[2]]
    C = tuk$U[[3]]
    print("A Done------------------")
    
    ##### update W
    #W = glmnet_mat_regu(Y_covar,A,lambda)
    W = solve(t(A)%*%A)%*%t(A)%*%Y_covar
    
    ##### update B
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
    
    re = glm_mat_regu(t(B),t(Y_2),t(Y_covar),t(G_AC2),beta,t(V))
    
    #Y_2_and_Start = cbind(Y_2,B)
    #re = glm_mat(t(Y_2_and_Start),t(G_AC2))
    #re = glm_mat(t(Y_2),start = t(B),t(G_AC2))
    
    B = t(re[[1]])
    cross_entropy[4*n - 2] = re[[2]]  + alpha*norm(Y_covar - A%*%W,type = 'F')
    ## orthogonal B*
    U = ttl(G,list(A,B,C),ms = c(1,2,3))
    tuk = tucker(U, ranks = core_shape)
    G = tuk$Z
    A = tuk$U[[1]]
    B = tuk$U[[2]]
    C = tuk$U[[3]]
    print("B Done------------------")
    
    ##### update V
    #V = glmnet_mat_regu(Y_covar,B,gamma)
    V = solve(t(B)%*%B)%*%t(B)%*%Y_covar
    
    ###### update C
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
    
    #Y_3_and_Start = cbind(Y_3,C)
    #re = glm_mat(t(Y_3_and_Start),t(G_AB3))
    re = glm_mat(t(Y_3),start = t(C),t(G_AB3))
    
    C = t(re[[1]])
    cross_entropy[4*n - 1] = -re[[2]] + alpha*norm(Y_covar - A%*%W,type = 'F') + beta*norm(Y_covar - B%*%V,type = 'F')
    ## orthogonal C*
    U = ttl(G,list(A,B,C),ms = c(1,2,3))
    tuk = tucker(U, ranks = core_shape)
    G = tuk$Z
    A = tuk$U[[1]]
    B = tuk$U[[2]]
    C = tuk$U[[3]]
    print("C Done------------------")
    
    M_long = matrix(0,nrow = d1*d2*d3, ncol = r1*r2*r3)
    m=1
    ## update G
    for(k in 1:d3){
      for(j in 1:d2){
        for(i in 1:d1){
          M_long[m,] = kronecker_list(list(C[k,],B[j,],A[i,]))
          m = m + 1
        }
      }
    }
    
    mod_re = glm_modify(as.vector(ts@data), M_long, as.vector(G@data))
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    cross_entropy[4*n] = - mod_re[[2]] + alpha*norm(Y_covar - A%*%W,type = 'F') + beta*norm(Y_covar - B%*%V,type = 'F')
    print("G Done------------------")
    print(n)
    #if(cross_entropy[4*n] - cross_entropy[4*(n-1)+1] <= 0.005) {print('break');break}
  }
  return(list(A = A,B = B,C = C,G = G,cross_entropy = cross_entropy))
}

up = update_binary(tensor, features , c(5,5,6), rep(0.005,4), 3)

plot(up$cross_entropy)


dim(features)[2]









