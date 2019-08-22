#**Setting for the workpath and clean the environment**

#rm(list = ls())
#setwd("/Users/March/Desktop/Tensor/code/Binarytensor/")


#**Source the functions**
# library(knitr)
# source_rmd <- function(file, local = FALSE, ...){
#   options(knitr.duplicate.label = 'allow')
#   
#   tempR <- tempfile(tmpdir = ".", fileext = ".R")
#   on.exit(unlink(tempR))
#   knitr::purl(file, output=tempR, quiet = TRUE)
#   
#   envir <- globalenv()
#   source(tempR, local = envir, ...)
# }
# 
# source_rmd("upgrade_new_new.Rmd")

source("upgrade_new_new.R")


#**The function of the whole pipline**

binary_tensor_fctr_km_pip = function(tensor,k,maxiter){
  #tensor can be an array
  #k is the rank of tucker factorization
  
  #replace NaN to 0
  tensor[is.na(tensor)] = 0
  #set tensor as a tensor class
  tensor = as.tensor(tensor)
  
  #here use tucker factorization to find the initial set
  inital_tucker = tucker(as.tensor(10*(2*tensor@data-1)),ranks = k)
  C_0 = inital_tucker$Z
  M1_0 = as.matrix(inital_tucker$U[[1]])
  M2_0 = as.matrix(inital_tucker$U[[2]])
  M3_0 = as.matrix(inital_tucker$U[[3]])
  
  likelihood = cal_likelihood(tensor,C_0,M1_0,M2_0,M3_0)#the first likelihood
  modes = list("U"=C_0,"M1"=M1_0,"M2" =M2_0,"M3"=M3_0)#the initial set
  
  for (iter in 2:maxiter) {
    modes_new = upgrade_glm_tucker(tensor,modes$U,modes$M1,modes$M2,modes$M3)
    likelihood = c(likelihood,modes_new$likelihood)
    print(iter)
    if((likelihood[(1+4*(iter-1))]-likelihood[(4*(iter-2)+1)])<= 0.005){
      #compare the full upgrade likelihood
      modes = modes_new
      break
    }
    modes = modes_new
  }
  
  M1_km = kmeans(modes$M1,centers = k[1])
  M2_km = kmeans(modes$M2,centers = k[2])
  M3_km = kmeans(modes$M3,centers = k[3])
  
  return(list("modes"=modes,"M1_CL" = M1_km,"M2_CL" = M2_km
              ,"M3_CL" = M3_km,"likelihood"=likelihood ))
}

#**test**
dnoations = readMat("dnations.mat")
dnoations_tr = dnoations$R
test_2220  =binary_tensor_fctr_km_pip(dnoations_tr,c(2,2,2),20)
