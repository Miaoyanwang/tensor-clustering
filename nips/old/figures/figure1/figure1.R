rm(list=ls())
require(tensorsparse)

########################################
###   Figure 1    ######################
########################################

n=60;p=60;q=60;k=5;r=5;l=5;error=10
set.seed(1)
data = get.data(n,p,q,k,r,l,error)
truth = data$truthX
input = data$x
output_k_means = classify2(input,k,r,l,max.iter = 0)
output = classify2(input,k,r,l)
plot_tensor(input)
Sys.sleep(1)
plot_tensor(truth)
Sys.sleep(1)
plot_tensor(output$judgeX)
Sys.sleep(1)
plot_tensor(output_k_means$judgeX)
cat("The evaluation of k means:\n")
sparse.evaluate(output_k_means,data)
cat("The evaluation of our method:\n")
sparse.evaluate(output,data)
