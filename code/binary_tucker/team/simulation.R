source('functions_synthesis.r')


seed = 24;d = seq(20,70,10);r = rep(3,7); p1 = rep(p[i],7); p2 = rep(p[i],7); 
p3 = rep(p[i],7); 
dis = 'gaussian';gs_mean = 0;gs_sd = 10;dup=2;Nsim =50; linear = FALSE; cons = 'penalty';lambda = 1; solver = 'CG'

 data = gene_data_all(seed,c(20,20,20),c(3,3,3), p1 = p_1,  p2 = p_2,p3 = p_3, dis, gs_mean, gs_sd,unf_a, unf_b, dup,signal)


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

# re$RMSE = re$RMSE*sqrt(re$d^3)
# re$rate = sqrt(re$rate)

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




