# May 7 2012
# implement the simulation
set.seed(1)
library(xtable)
source("Table6FunctionLambda.R")
biclustauto<-Do(100,lambda=c(0,10,20,30,40,50,60,70,80))

table<-rbind(c("Method","Totalzero","CorrectZero","CorrectOne","Incorrect"),c(paste("biclust auto = ",round(biclustauto$lambda,digits=3),sep=","),biclustauto$totalzero,biclustauto$correctzero,biclustauto$correctone,biclustauto$totalincorrect))

table2<-rbind(c("Method","Totalzero","CorrectZero","CorrectOne","Incorrect"),c(paste("biclust auto = ",round(biclustauto$lambda,digits=3),sep=","),biclustauto$sdtotalzero,biclustauto$sdcorrectzero,biclustauto$sdcorrectone,biclustauto$sdtotalincorrect))


print(xtable(table),file="a2.txt")
print(xtable(table2),file="b2.txt")
