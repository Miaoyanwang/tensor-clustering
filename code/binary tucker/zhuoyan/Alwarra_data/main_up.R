### apply on real data
# 1. Kinship (Nickel et al., 2011, alyawarradata.mat):
# This is a 104 ¡Á 104 ¡Á 26 binary tensor consisting of 26 types of relations among a
# set of 104 individuals in Australian Alyawarra tribe. The data was first collected by
# Denham and White (2005) to study the kinship system in the Alyawarra language. 
# The tensor entry Y(i,j,k) is 1 if individual i used the kinship term k to refer to
# individual j, and 0 otherwise.
# 
# data$Rs: 104 ¡Á 104 ¡Á 26
# data$features: 14 features for the 104 individuals

library(R.matlab)
data = readMat('alyawarradata.mat')
atts = data$atts
tensor = data$Rs
features = data$features
ss = data$datass

?isSymmetric
?isSymmetric.matrix

### To check whether it is symmetric
for(i in 1:dim(tensor)[3]){
  print(isSymmetric(tensor[,,i]))
}


#up = update_binary(tensor,c(5,5,6),5)
up = update_binary(tensor, features , c(5,5,6), rep(0.005,4), 3)

plot(up$cross_entropy,type = 'b')
up$cross_entropy

              




















