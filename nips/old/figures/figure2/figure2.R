rm(list=ls())
require(tensorsparse)

########################################
###   Figure 2    ######################
########################################

n=20;p=20;q=20;k=5;r=6;l=7;error=3;sparse.percent=0.8
set.seed(1)
data = get.data(n,p,q,k,r,l,error,sparse.percent=0)
truth = data$truthX
input = data$x
plot_tensor(truth)

selected_krl = choosekrl_bic(input,k=1:8,r=1:8,l=1:8)$estimated_krl
lambda = chooseLambda(input,selected_krl[1,1],selected_krl[1,2],selected_krl[1,3])$lambda

output = classify2(input,selected_krl[1,1],selected_krl[1,2],selected_krl[1,3],lambda = lambda)
plot_tensor(input)
Sys.sleep(1)
plot_tensor(truth)
Sys.sleep(1)
plot_tensor(output$judgeX)


cat("The evaluation of our method:\n")
sparse.evaluate(output,data)
