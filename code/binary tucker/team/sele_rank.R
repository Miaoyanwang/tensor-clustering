rm(list = ls())
source('functions_synthesis_all.r')

data = gene_data_all(24, whole_shape = c(20,20,20), core_shape = c(4,4,5),
                    p1 = NULL,p2 = NULL, p3 = NULL,dis = 'gaussian', gs_mean = 0,gs_sd = 10,
                    unf_a = 0,unf_b = 1,dup = 5, signal = 10)
# data$tsr[[1]] - data$tsr[[2]]
# data$U
re = lapply(X = data$tsr, FUN = sele_rank, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,
         rank = seq(4,6), Nsim = 30,linear = FALSE, cons = 'non')


-----   lapply(re,toString)
re = unlist(lapply(re,toString))

write.csv(re,'rank445.csv',row.names = F)




####   ------    import 

# c = read.csv('a.csv')
# 
# lapply(strsplit(as.character(c$x),split = ','), as.numeric)
# 


data2 = gene_data_all(24, whole_shape = c(20,20,20), core_shape = c(5,5,8),
                     p1 = NULL,p2 = NULL, p3 = NULL,dis = 'gaussian', gs_mean = 0,gs_sd = 10,
                     unf_a = 0,unf_b = 1,dup = 5, signal = 10)


re2 = lapply(X = data2$tsr, FUN = sele_rank, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,
            rank = seq(4,9), Nsim = 30,linear = FALSE, cons = 'non')


re2 = unlist(lapply(re2,toString))

write.csv(re2,'rank558.csv',row.names = F)


# 
# lapply(strsplit(a,split = ','), as.numeric)
# 
# 
# 
# lapply(strsplit(b$`12`,split = ','), as.numeric)

data3 = gene_data_all(24, whole_shape = c(30,30,30), core_shape = c(8,8,10),
                      p1 = NULL,p2 = NULL, p3 = NULL,dis = 'gaussian', gs_mean = 0,gs_sd = 10,
                      unf_a = 0,unf_b = 1,dup = 5, signal = 10)


re3 = lapply(X = data3$tsr, FUN = sele_rank, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,
             rank = seq(7,11), Nsim = 25,linear = FALSE, cons = 'non')


re3 = unlist(lapply(re3,toString))

write.csv(re3,'rank8810.csv',row.names = F)

#####----   read
re = read.csv('rank445.csv')
re2 = read.csv('rank558.csv')
re3 = read.csv('rank8810.csv')

rank = cbind(re,re2,re3)
colnames(rank) = c('4,4,5','5,5,8','8,8,10')
rownames(rank) = c('est1','est2','est3','est4','est5')



write.csv(rank,'rank_sele.csv',row.names = TRUE)

