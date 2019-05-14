
# implement the simulation
set.seed(1)
library(xtable)
source("Table6Function.R")
biclust0<-Do(100,method="biclustering",40)
biclust1<-Do(100,method="biclustering",60)
biclust2<-Do(100,method="biclustering",80)
biclust3<-Do(100,method="biclustering",100)
biclust4<-Do(100,method="biclustering",120)

ssvd<-Do(100,method="ssvd",0)
plaid<-Do(100,method="plaid1",0)
plaid2<-Do(100,method="plaid2",0)

table<-rbind(c("Method","Totalzero","CorrectZero","CorrectOne","Incorrect"),c("biclust lambda=40",biclust0$totalzero,biclust0$correctzero,biclust0$correctone,biclust0$totalincorrect),c("biclust lambda=60",biclust1$totalzero,biclust1$correctzero,biclust1$correctone,biclust1$totalincorrect),c("biclust lambda=80",biclust2$totalzero,biclust2$correctzero,biclust2$correctone,biclust2$totalincorrect),c("biclust lambda=100",biclust3$totalzero,biclust3$correctzero,biclust3$correctone,biclust3$totalincorrect),c("biclust lambda=120",biclust4$totalzero,biclust4$correctzero,biclust4$correctone,biclust4$totalincorrect),c("ssvd",ssvd$totalzero,ssvd$correctzero,ssvd$correctone,ssvd$totalincorrect),c("plaid",plaid$totalzero,plaid$correctzero,plaid$correctone,plaid$totalincorrect),c("plaid constant",plaid2$totalzero,plaid2$correctzero,plaid2$correctone,plaid2$totalincorrect))

table2<-rbind(c("Method","Totalzero","CorrectZero","CorrectOne","Incorrect"),c("biclust lambda=40",biclust0$sdtotalzero,biclust0$sdcorrectzero,biclust0$sdcorrectone,biclust0$sdtotalincorrect),c("biclust lambda=60",biclust1$sdtotalzero,biclust1$sdcorrectzero,biclust1$sdcorrectone,biclust1$sdtotalincorrect),c("biclust lambda=80",biclust2$sdtotalzero,biclust2$sdcorrectzero,biclust2$sdcorrectone,biclust2$sdtotalincorrect),c("biclust lambda=100",biclust3$sdtotalzero,biclust3$sdcorrectzero,biclust3$sdcorrectone,biclust3$sdtotalincorrect),c("biclust lambda=120",biclust4$sdtotalzero,biclust4$sdcorrectzero,biclust4$sdcorrectone,biclust4$sdtotalincorrect),c("ssvd",ssvd$sdtotalzero,ssvd$sdcorrectzero,ssvd$sdcorrectone,ssvd$sdtotalincorrect),c("plaid",plaid$sdtotalzero,plaid$sdcorrectzero,plaid$sdcorrectone,plaid$sdtotalincorrect),c("plaid constant",plaid2$sdtotalzero,plaid2$sdcorrectzero,plaid2$sdcorrectone,plaid2$sdtotalincorrect))

print(xtable(table),file="a.txt")
print(xtable(table2),file="b.txt")
