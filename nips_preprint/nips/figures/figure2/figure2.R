rm(list=ls())
require(tensorsparse)

########################################
###   Figure 2    ######################
########################################

n=60;p=60;q=60;k=6;r=6;l=7;error=1;sparse.percent=0.8
set.seed(1)
data = get.data(n,p,q,k,r,l,error,TRUE,sparse.percent)
truth = data$truthX
input = data$x

selected_krl = choosekrl_bic(input,k=5:8,r=5:8,l=5:8)$estimated_krl
lambda = chooseLambda(input,selected_krl[1,1],selected_krl[1,2],selected_krl[1,3])$lambda

output = classify2(input,selected_krl[1,1],selected_krl[1,2],selected_krl[1,3],lambda = lambda)
plot_tensor(input)
Sys.sleep(1)
plot_tensor(truth)
Sys.sleep(1)
plot_tensor(output$judgeX)


cat("The evaluation of our method:\n")
sparse.evaluate(output,data)
