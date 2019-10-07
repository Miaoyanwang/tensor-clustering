#' Simulation function
#'
#' This is the function for pack of simulations using the updating function to validate the rates
#' @param seed        seed for generated data
#' @param d           the vector contains different tensor size, please see "details"
#' @param r           the vector contains different tensor rank, please see "details"
#' @param p1          the vector contains  different number of features on mode 1
#' @param p2          the vector contains  different number of features on mode 2
#' @param dis         distribution of generated core tensor
#' @param gs_mean     please see "details"
#' @param gs_sd       please see "details"
#' @param unf_a       please see "details"
#' @param unf_b       please see "details"
#'
#' @param dup         the duplicates of simulated tensor for each ground truth
#' @param Nsim        please see  help documents of "update_binary"
#'
#' @param linear      please see  help documents of "update_binary"
#' @param cons        please see  help documents of "update_binary"
#' @param lambda      please see  help documents of "update_binary"
#' @param solver      please see  help documents of "update_binary"
#'
#' @return     a list containing following components:
#'
#'                    \code{RMSE} RMSE for different \eqn{d}
#'
#'                    \code{rate} RMSE for different \eqn{d}
#'
#'
#' @details    In this scenario, for each \eqn{d_i} in vector d, we consider
#'               \deqn{d_1 = d_2 = d_3 = d}
#'              r is the same
#'
#' @export



conv_rate = function(seed,d,r, p1, p2, dis,gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 1,
                     dup, Nsim, linear = TRUE, cons = 'vanilla' ,lambda = 1,
                     solver = NULL){
  #cons can be "non","vanilla","penalty"
  rate = rep(0,length(d))
  RMSE = rep(0,length(d))
  for (i in 1:length(d)) {
    data = gene_data(seed,rep(d[i],3), rep(r[i],3), p1[i], p2[i], dis, gs_mean, gs_sd, unf_a, unf_b, dup)
    X_covar1 = data$X_covar1
    X_covar2 = data$X_covar2
    C_ts = data$C_ts
    U = data$U
    tsr = data$tsr
    RMSEi = rep(0,dup)
    for (j in 1:dup) {
      upp = update_binary(tsr = tsr[[j]], X_covar1 = X_covar1, X_covar2 = X_covar2,
                          core_shape =  rep(r[i],3), Nsim, linear, cons, lambda = lambda,
                          alpha = 10*max(abs(U)), solver = NULL)

      C_ts_est = ttl(upp$G,list(upp$W1,upp$W2,upp$C),ms = c(1,2,3))@data
      #U_est = ttl(upp$G,list(X_covar1%*%upp$W1,X_covar2%*%upp$W2,upp$C),ms = c(1,2,3))@data
      RMSEi[j] = sqrt(sum((C_ts_est - C_ts)^2)/(d[i]^3))
      print(paste(j,"-th observation ---- when dimension is ",d[i],"-- rank is ",r[i]," ---------"))
    }
    RMSE[i] = mean(RMSEi)
    rate[i] = r[i]^2*(d[i] + p1[i] + p2[i])/d[i]^3
  }
  return(list(RMSE = RMSE, rate = rate))
}
