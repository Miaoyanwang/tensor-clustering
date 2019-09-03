#' Solve unsupervised binary tucker decomposition with no covariate or covariates are identity matrix
#'
#' This is main updating function for the update of factor matrces and  core tensor
#'             \eqn{U = G \times_1 A \times_2 B \times_3 C}
#' @param tsr         the data vector
#' @param core_shape  the tucker rank of the given tensor
#' @param Nsim        simulation times
#' @param cons        the constrain method, "non" for without constrain, "vanilla" for global scale down once at each iteration,
#'
#'                    "penalty" for add log-barrier penalty to object function.
#' @param lambda      penalty coefficient when use "penalty" constrain.
#' @param alpha       max norm consrain on ground truth
#' @param solver      solver for solving object function when using "penalty" constrain, see "details"
#' @return     a list containing following components:
#'
#'                    \code{W1} coefficient matrix on mode-1
#'
#'                    \code{W2} coefficient matrix on mode-2
#'
#'                    \code{C} factor matrix on mode-3
#'
#'                    \code{G} core tensor
#'
#'                    \code{lglk} the lig-likelihood
#'
#' @details    When you choose constrain to be "penalty". It actually add log-barrier as regularizer to
#'            general object function (negative log-likelihood). It use solver in function "optim" to
#'            solve the objective function. You can choose any solver contained in function "optim".
#'
#' @export


#################  update
###--------   unsupervised
update_binary_un = function(tsr, core_shape, Nsim, cons, lambda = 1, alpha = 1, solver = NULL){
    ## get initialization
    tsr1 = 10*(2*tsr - 1)
    tsr1 = as.tensor(tsr1)
    tsr = as.tensor(tsr)
    tckr = tucker(tsr1, ranks = core_shape)
    A = tckr$U[[1]] ; B = tckr$U[[2]] ; C = tckr$U[[3]]
    G = tckr$Z
    d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
    r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
    Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
    Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
    Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
    
    violate = c() ## which iteration violate constrain
    
    lglk = rep(0,4*Nsim)
    for(n in 1:Nsim){
        
        ###### update A
        G_BC = ttl(G, list(B,C), ms = c(2,3))
        G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
        
        re = glm_mat(t(Y_1),start = t(A),t(G_BC1))
        
        if(dim(A)[2]==1) A=as.matrix(re[[1]])
        else A = t(re[[1]])
        
        lglk[4*n - 3] = re[[2]]
        ## orthogonal A*
        U = ttl(G,list(A,B,C),ms = c(1,2,3))
        tuk = tucker(U, ranks = core_shape)
        G = tuk$Z
        A = tuk$U[[1]]
        B = tuk$U[[2]]
        C = tuk$U[[3]]
        print("A Done------------------")
        
        ##### update B
        G_AC = ttl(G, list(A,C), ms = c(1,3))
        G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
        
        re = glm_mat(t(Y_2),start = t(B),t(G_AC2))
        
        
        if(dim(B)[2]==1) B=as.matrix(re[[1]])
        else B = t(re[[1]])
        
        lglk[4*n - 2] = re[[2]]
        ## orthogonal B*
        U = ttl(G,list(A,B,C),ms = c(1,2,3))
        tuk = tucker(U, ranks = core_shape)
        G = tuk$Z
        A = tuk$U[[1]]
        B = tuk$U[[2]]
        C = tuk$U[[3]]
        print("B Done------------------")
        
        ###### update C
        G_AB = ttl(G, list(A,B), ms = c(1,2))
        G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
        
        re = glm_mat(t(Y_3),start = t(C),t(G_AB3))
        
        if(dim(C)[2]==1) C=as.matrix(re[[1]])
        else C = t(re[[1]])
        
        lglk[4*n - 1] = re[[2]]
        
        #########-----------------------------------------------
        ###  then we apply out constrain
        ######---- differnent version of contrains
        U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
        
        if(cons == 'non'){U = U}
        else if(max(abs(U)) <= alpha){U = U}
        else if(cons == 'vanilla'){
            U = U/max(abs(U))*alpha
            print("Violate constrain ------------------")
            violate = c(violate,n)
        }
        else{
            U = U/max(abs(U))*(alpha-0.01)
            print("Violate constrain ------------------")
            violate = c(violate,n)
        }
        
        ## orthogonal C*
        
        U = as.tensor(U)
        tuk = tucker(U, ranks = core_shape)
        G = tuk$Z
        A = tuk$U[[1]]
        B = tuk$U[[2]]
        C = tuk$U[[3]]
        print("C Done------------------")
        
        ##### update G
        M_long = kronecker_list(list(C,B,A)) ## form M_long
        
        if(cons == 'penalty'){
            mod_re = optim(par = as.vector(G@data),loss,loss_gr,y = as.vector(tsr@data),
            X = M_long, lambda = lambda, alpha = alpha, method = solver)
            coe = mod_re$par
            G = as.tensor(array(data = coe,dim = core_shape))
            lglk[4*n] = -mod_re$value
        }
        else {
            mod_re = glm_modify(as.vector(tsr@data), M_long, as.vector(G@data))
            coe = mod_re[[1]]
            G = as.tensor(array(data = coe,dim = core_shape))
            lglk[4*n] = mod_re[[2]]
        }
        
        print("G Done------------------")
        
        print(paste(n,"-th  iteration ---- when dimension is ",d1,"-- rank is ",r1," -----------------"))
        #if(abs(lglk[4*n-1] - lglk[4*n-2]) <= 0.0005) break
        #if use vanilla constrain, it may happen than on the first iteration, the likelihood would drop.
        if(cons == "vanilla") {
            if(lglk[4*n] - lglk[4*n-1] <= 0.00005 & n >= 20) break
            
        }
        else{
            if(lglk[4*n-1] - lglk[4*n-2] <= 0.00005) break
        }
        
    }
    return(list(A = A,B = B,C = C,G = G,lglk = lglk, violate = violate))
}


