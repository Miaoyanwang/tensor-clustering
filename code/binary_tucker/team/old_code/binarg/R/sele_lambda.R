#' Select penalty coefficient in constrained optimization in penalised algorithm
#'
#' This is the function for selecting lambda in penalised optimization in supervised setting
#' @param seed        seed for generated data
#' @param lambda      a vector containing different penalty coefficient
#' @param ...         other passed argument(see "details")

#' @return     a list containing different RMSE when choosing different lambda, if we set dimension
#'             to be a vector, the element in list is a vector containing different RMSE correspond to
#'             different dimension(so is rank, p1, p2,etc.)
#'
#' @details    other passed argument please see "conv_rate"
#'
#' @export


sele_lambda = function(seed, lambda, ...){
  #lambda = as.list(lambda)
  re = lapply(lambda, FUN = conv_rate, seed = seed, ...)
  re = lapply(seq(length(re)), function(x) re[[x]]$RMSE)
  return(re)
}



