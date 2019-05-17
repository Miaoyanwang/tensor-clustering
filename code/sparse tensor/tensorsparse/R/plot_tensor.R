#' draw 3d plot of tensor
#' 
#' draw 3d plot of tensor
#' @param tensor a three-dimensional array
#' 
#' @export
plot_tensor=function(tensor){
  
  n=prod(dim(tensor))   
  color_choice=min(round(prod(dim(tensor))/(6*6*6)),100)+1
  marker=viridis_pal(option = "B")(color_choice)
  
  position=positionfun(dim(tensor))$position
  quan=c(quantile(tensor,(0:color_choice)/color_choice))
  col=tensor
  for(i in 1:color_choice){
    col[(tensor>=quan[i])&(tensor<quan[i+1])]=marker[i]
  }
  col[tensor==quan[i+1]]=marker[i]
  
  plot3d(position[,1],position[,2],position[,3],col=col,alpha=0.3,size=5,xlab="",ylab="",zlab="")
}
