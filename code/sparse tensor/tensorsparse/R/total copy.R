## simulations
## block model
data = get.data(50,50,50,4,4,4,error=1)
krl=choosekrl_bic(data$x,2:6,2:6,2:6)$estimated_krl
