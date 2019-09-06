#######---------  codes for visualization
re_non = read.csv('~/GitHub/tensor-clustering/code/binary tucker/simulation_result/Unsupervised/seed24_non_constrain.csv',header = TRUE)

re_va = read.csv('~/GitHub/tensor-clustering/code/binary tucker/simulation_result/Unsupervised/seed24_vanilla.csv',header = TRUE)

re_mom = read.csv('~/GitHub/tensor-clustering/code/binary tucker/simulation_result/Unsupervised/seed24_momentum.csv',header = T)

library(ggplot2)

ggplot(re_non, aes(x = d, y = RMSE)) + geom_line(aes(color = as.factor(rank)),size = 1.5)  +
  geom_point(aes(shape =  as.factor(rank)), size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggplot(re_non, aes(x = rate, y = RMSE)) + geom_line(aes(color = as.factor(rank)),size = 2)  +
  geom_point(aes(shape =  as.factor(rank)), size = 3)  + 
  geom_smooth(method = 'lm',formula = y~x) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


####----  reproduce figure 2 and 6 result
##--============ figure 2
data = gene_data_un(37, whole_shape = c(10,10,10), core_shape = c(1,1,1),dis = 'gaussian', gs_mean = 0,gs_sd = 1,unf_a = 0,unf_b = 1, 1,signal=10)
U = data$U
tsr = data$ts[[1]]
upp = update_binary_un(tsr, core_shape = c(1,1,1), Nsim = 5, cons = 'vanilla', lambda = 0.1, alpha = max(abs(U)), solver = "CG")

### U_est: estimated ground truth
U_est = ttl(upp$G,list(upp$A,upp$B,upp$C),ms = c(1,2,3))@data

plot(as.vector(U),as.vector(U_est));abline(0,1)   ## plot U vs U_est


plot(as.vector(inv.logit(U)),as.vector(inv.logit(U_est)));abline(0,1)  ## plot U vs U_est in logical scale

##--============ figure 6

data = gene_data_un(37, whole_shape = c(20,20,20), core_shape = c(3,3,3),dis = 'uniform',
                    gs_mean = 0,gs_sd = 10,unf_a = 0,unf_b = 100, 1)
U = data$U
tsr = data$ts[[1]]
upp = update_binary_un(tsr, core_shape = c(3,3,3), Nsim = 20, cons = 'non', lambda = 1, alpha = 1, solver = NULL)

### U_est: estimated ground truth
U_est = ttl(upp$G,list(upp$A,upp$B,upp$C),ms = c(1,2,3))@data


plot(as.vector(U),as.vector(U_est));abline(0,1)   ## plot U vs U_est
plot(as.vector(inv.logit(U)),as.vector(inv.logit(U_est)));abline(0,1)  ## plot U vs U_est in logical scale

########################--------------------------------------------------------------------
####################################################################################

###---  newest simulation code, still waiting for result
rm(list = ls())
source('functions_synthesis_all.r')

####-------------------------------  convergence rate
p = seq(5,19,2)
re = c()

for(i in seq(length(p))){
  con = conv_rate(seed = 24,d = seq(20,70,10),r = rep(3,7), p1 = rep(p[i],7), p2 = rep(p[i],7), p3 = NULL,
                  dis = 'gaussian',gs_mean = 0,gs_sd = 10, 
                  dup=2, signal = 10,Nsim =50, linear = FALSE, cons = 'penalty' ,lambda = 1,
                  solver = 'CG')
  
  con = as.data.frame(con)
  con$d = seq(20,70,10); con$rank = rep(3,7); con$p1 = rep(p[i],7); con$p2 = rep(p[i],7)
  re = rbind(re,con)
}




####  supervised simulation result for most recent note
###---  "GitHub\tensor-clustering\note\zhuoyan\supervised_sim2"


source('functions_synthesis.r')

####-------------------------------  convergence rate
p = seq(5,19,2)
re = c()

for(i in seq(length(p))){
  con = conv_rate(seed = 24,d = seq(20,70,10),r = rep(3,7), p1 = rep(p[i],7), p2 = rep(p[i],7), 
                  dis = 'gaussian',gs_mean = 0,gs_sd = 10, 
                   dup=2, Nsim =50, linear = FALSE, cons = 'penalty' ,lambda = 1,
                   solver = 'CG')
  
  con = as.data.frame(con)
  con$d = seq(20,70,10); con$rank = rep(3,7); con$p1 = rep(p[i],7); con$p2 = rep(p[i],7)
  re = rbind(re,con)
}


write.csv(re,"rank3.csv",row.names = FALSE)

## visual
re = read.csv('rank3.csv',header = TRUE)

re$RMSE = re$RMSE*sqrt(re$d^3)
re$rate = sqrt(re$rate)

library(ggplot2)

#re[re$d != 20,]
ggplot(re, aes(x = d, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 1.5)  +
  geom_point( size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggplot(re[re$d != 20 ,], aes(x = d, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 1.5)  +
  geom_point( size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggplot(re, aes(x = p1, y = RMSE)) + geom_line(aes(color = as.factor(d)),size = 1.5)  +
  geom_point(aes(shape =  as.factor(d)), size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(re[re$d != 20,], aes(x = p1, y = RMSE)) + geom_line(aes(color = as.factor(d)),size = 1.5)  +
  geom_point(aes(shape =  as.factor(d)), size = 3) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(re, aes(x = rate, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 2)  +
  geom_point(size = 3)  + 
  geom_smooth(method = 'lm',formula = y~x) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(re[re$d != 20 ,], aes(x = rate, y = RMSE)) + geom_line(aes(color = as.factor(p1)),size = 2)  +
  geom_point( size = 3)  + 
  geom_smooth(method = 'loess',formula = y~x) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))




