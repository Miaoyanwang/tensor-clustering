###Class Definition

#'S4 Class for a Tensor
#'
#'An S4 class for a tensor with arbitrary number of modes. The Tensor class extends the base "array" class to include additional tensor manipulation (folding, unfolding, reshaping, subsetting) as well as a formal class definition that enables more explicit tensor algebra.
#'
#'@section Slots:
#' \describe{
#'	\item{num_modes}{number of modes (integer)}
#'  \item{modes}{vector of modes (integer), aka sizes/extents/dimensions}
#'  \item{data}{actual data of the tensor, which can be 'array' or 'vector'}
#' }
#'@name Tensor-class
#'@rdname Tensor-class
#'@aliases Tensor Tensor-class
#'@docType class
#'@exportClass Tensor
#'@author James Li \email{jamesyili@@gmail.com}
#'@note All of the decompositions and regression models in this package require a Tensor input.
#'@references James Li, Jacob Bien, Martin T. Wells (2018). rTensor: An R Package for Multidimensional Array (Tensor) Unfolding, Multiplication, and Decomposition. Journal of Statistical Software, Vol. 87, No. 10, 1-31. URL: http://www.jstatsoft.org/v087/i10/.
#'@seealso \code{\link{as.tensor}}
setClass("Tensor",
         representation(num_modes = "integer", modes = "integer", data="array"),
         validity = function(object){
           num_modes <- object@num_modes
           modes <- object@modes
           errors <- character()
           if (any(modes <= 0)){
             msg <- "'modes' must contain strictly positive values; if any mode is 1, consider a smaller num_modes"
             errors <- c(errors, msg)
           }
           if(length(errors)==0) TRUE else errors
         })


setMethod(f="initialize",
          signature="Tensor",
          definition = function(.Object, num_modes=NULL, modes=NULL, data=NULL){
            if(is.null(num_modes)){
              if (is.vector(data)) num_modes <- 1L
              else{num_modes <- length(dim(data))}
            }
            if(is.null(modes)){
              if (is.vector(data)) modes <- length(data)
              else{modes <- dim(data)}
            }
            .Object@num_modes <- num_modes
            .Object@modes <- modes
            if (is.vector(data)){
              .Object@data <- array(data,dim=modes)
            }else{
              .Object@data <- data
            }
            validObject(.Object)
            .Object
          })


setMethod(f="show",
          signature="Tensor",
          definition=function(object){
            cat("Numeric Tensor of", object@num_modes, "Modes\n", sep=" ")
            cat("Modes: ", object@modes, "\n", sep=" ")
            cat("Head of entries: \n")
            print(head(object@data))
          })

###Generic Definitions
#'Mode Getter for Tensor
#'
#'Return the vector of modes from a tensor
#'
#'@name dim-methods
#'@details \code{dim(x)}
#'@export
#'@aliases dim,Tensor-method
#'@docType methods
#'@rdname dim-methods
#'@param x the Tensor instance
#'@return an integer vector of the modes associated with \code{x}
#'@examples
#'tnsr <- rand_tensor()
#'dim(tnsr)
setMethod(f = 'dim', signature = 'Tensor', function(x) x@modes)



#'Tensor Unfolding
#'
#'Unfolds the tensor into a matrix, with the modes in \code{rs} onto the rows and modes in \code{cs} onto the columns. Note that \code{c(rs,cs)} must have the same elements (order doesn't matter) as \code{x@@modes}. Within the rows and columns, the order of the unfolding is determined by the order of the modes. This convention is consistent with Kolda and Bader (2009).
#'
#For Row Space Unfolding or m-mode Unfolding, see \code{\link{rs_unfold-methods}}. For Column Space Unfolding or matvec, see \code{\link{cs_unfold-methods}}.

#'
#\code{\link{vec-methods}} returns the vectorization of the tensor.

#'
#'@details \code{unfold(tnsr,row_idx=NULL,col_idx=NULL)}
#'@export
#'@docType methods
#'@name unfold-methods
#'@rdname unfold-methods
#'@aliases unfold unfold,Tensor-method
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009, Vol. 51, No. 3 (September 2009), pp. 455-500. URL: https://www.jstor.org/stable/25662308.
#'@param tnsr the Tensor instance
#'@param row_idx the indices of the modes to map onto the row space
#'@param col_idx the indices of the modes to map onto the column space
#'@return matrix with \code{prod(row_idx)} rows and \code{prod(col_idx)} columns
#@seealso \code{\link{k_unfold-methods}} and \code{\link{matvec-methods}}

#'@examples
#'tnsr <- rand_tensor()
#'matT3<-unfold(tnsr,row_idx=2,col_idx=c(3,1))
setGeneric(name="unfold",
           def=function(tnsr,row_idx,col_idx){standardGeneric("unfold")})

setMethod("unfold", signature="Tensor",
          definition=function(tnsr,row_idx=NULL,col_idx=NULL){
            #checks
            rs <- row_idx
            cs <- col_idx
            if(is.null(rs)||is.null(cs)) stop("row and column indices must be specified")
            num_modes <- tnsr@num_modes
            if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
            if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
            perm <- c(rs,cs)
            if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
            modes <- tnsr@modes
            mat <- tnsr@data
            new_modes <- c(prod(modes[rs]),prod(modes[cs]))
            #rearranges into a matrix
            mat <- aperm(mat,perm)
            dim(mat) <- new_modes
            as.tensor(mat)
          })


#'@docType methods
#'@name rs_unfold-methods
#'@details \code{rs_unfold(tnsr,m=NULL)}
#'@param tnsr Tensor instance
#'@param m mode to be unfolded on
#'@export
#'@rdname rs_unfold-methods
#'@aliases rs_unfold rs_unfold,Tensor-method
####aliases rs_unfold,ANY-method
setGeneric(name="rs_unfold",
           def=function(tnsr,m){standardGeneric("rs_unfold")})

#'@rdname rs_unfold-methods
#'@aliases rs_unfold,Tensor-method
setMethod("rs_unfold", signature="Tensor",
          definition=function(tnsr,m=NULL){
            if(is.null(m)) stop("mode m must be specified")
            num_modes <- tnsr@num_modes
            rs <- m
            cs <- (1:num_modes)[-m]
            unfold(tnsr,row_idx=rs,col_idx=cs)
          })



#'
#'@docType methods
#'@name cs_unfold-methods
#'@details \code{cs_unfold(tnsr,m=NULL)}
#'@param tnsr Tensor instance
#'@param m mode to be unfolded on
#'@export
#'@rdname cs_unfold-methods
#'@aliases cs_unfold cs_unfold,Tensor-method
setGeneric(name="cs_unfold",
           def=function(tnsr,m){standardGeneric("cs_unfold")})

#'@rdname cs_unfold-methods
#'@aliases cs_unfold,Tensor-method
setMethod("cs_unfold", signature="Tensor",
          definition=function(tnsr,m=NULL){
            if(is.null(m)) stop("mode m must be specified")
            num_modes <- tnsr@num_modes
            rs <- (1:num_modes)[-m]
            cs <- m
            unfold(tnsr,row_idx=rs,col_idx=cs)
          })
options(warn=1)

###Matrix Foldings

#'General Folding of Matrix
#'
#'General folding of a matrix into a Tensor. This is designed to be the inverse function to \code{\link{unfold-methods}}, with the same ordering of the indices. This amounts to following: if we were to unfold a Tensor using a set of \code{row_idx} and \code{col_idx}, then we can fold the resulting matrix back into the original Tensor using the same \code{row_idx} and \code{col_idx}.
#'@export
#'@details This function uses \code{aperm} as the primary workhorse.
#'@name fold
#'@rdname fold
#'@aliases fold
#'@param mat matrix to be folded into a Tensor
#'@param row_idx the indices of the modes that are mapped onto the row space
#'@param col_idx the indices of the modes that are mapped onto the column space
#'@param modes the modes of the output Tensor
#'@return Tensor object with modes given by \code{modes}
#seealso \code{\link{unfold-methods}}, \code{\link{k_fold}}, \code{\link{unmatvec}}
#'@seealso \code{\link{unfold-methods}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009, Vol. 51, No. 3 (September 2009), pp. 455-500. URL: https://www.jstor.org/stable/25662308.
#'@examples
#'tnsr <- new('Tensor',3L,c(3L,4L,5L),data=runif(60))
#'matT3<-unfold(tnsr,row_idx=2,col_idx=c(3,1))
#'identical(fold(matT3,row_idx=2,col_idx=c(3,1),modes=c(3,4,5)),tnsr)
fold <- function(mat, row_idx = NULL, col_idx = NULL, modes=NULL){
  #checks
  rs <- row_idx
  cs <- col_idx
  if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
  if(is.null(modes)) stop("Tensor modes must be specified")
  if(!is(mat,"Tensor")){
    if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
  }else{
    stopifnot(mat@num_modes==2)
    mat <- mat@data			
  }
  num_modes <- length(modes)
  stopifnot(num_modes==length(rs)+length(cs))
  mat_modes <- dim(mat)
  if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
  #rearranges into array
  iperm <- match(1:num_modes,c(rs,cs))
  as.tensor(aperm(array(mat,dim=c(modes[rs],modes[cs])),iperm))
}



#'@export
#'@name Ops-methods
#'@docType methods
#'@aliases Ops-methods Ops,Tensor,Tensor-method Ops,Tensor,array-method Ops,Tensor,numeric-method Ops,array,Tensor-method Ops,numeric,Tensor-method
#'@param e1 left-hand object
#'@param e2 right-hand object
#'@examples
#'tnsr <- rand_tensor(c(3,4,5))
#'tnsr2 <- rand_tensor(c(3,4,5))
#'tnsrsum <- tnsr + tnsr2
#'tnsrdiff <- tnsr - tnsr2
#'tnsrelemprod <- tnsr * tnsr2
#'tnsrelemquot <- tnsr / tnsr2
#'for (i in 1:3L){
#'	for (j in 1:4L){
#'		for (k in 1:5L){
#'			stopifnot(tnsrsum@@data[i,j,k]==tnsr@@data[i,j,k]+tnsr2@@data[i,j,k])
#'			stopifnot(tnsrdiff@@data[i,j,k]==(tnsr@@data[i,j,k]-tnsr2@@data[i,j,k]))
#'			stopifnot(tnsrelemprod@@data[i,j,k]==tnsr@@data[i,j,k]*tnsr2@@data[i,j,k])
#'			stopifnot(tnsrelemquot@@data[i,j,k]==tnsr@@data[i,j,k]/tnsr2@@data[i,j,k])
#'}
#'}
#'}
setMethod("Ops", signature(e1="Tensor", e2="Tensor"),
          definition=function(e1,e2){
            e1@data<-callGeneric(e1@data, e2@data)
            validObject(e1)
            e1
          })
setMethod("Ops", signature(e1="Tensor", e2="array"),
          definition=function(e1,e2){
            e1@data<-callGeneric(e1@data,e2)
            validObject(e1)
            e1
          })
setMethod("Ops", signature(e1="array", e2="Tensor"),
          definition=function(e1,e2){
            e2@data<-callGeneric(e1,e2@data)
            validObject(e2)
            e2
          })
setMethod("Ops", signature(e1="Tensor", e2="numeric"),
          definition=function(e1,e2){
            e1@data<-callGeneric(e1@data,e2)
            validObject(e1)
            e1
          })
setMethod("Ops", signature(e1="numeric", e2="Tensor"),
          definition=function(e1,e2){
            e2@data<-callGeneric(e1,e2@data)
            validObject(e2)
            e2
          })


#####
.is_zero_tensor <- function(tnsr){
  if (sum(tnsr@data==0)==prod(tnsr@modes)) return(TRUE)
  return(FALSE)
}

#####Special Tensors

#'Tensor with Random Entries
#'
#'Generate a Tensor with specified modes with iid normal(0,1) entries.
#'@export
#'@name rand_tensor
#'@rdname rand_tensor
#'@aliases rand_tensor
#'@param modes the modes of the output Tensor
#'@param drop whether or not modes equal to 1 should be dropped
#'@return a Tensor object with modes given by \code{modes}
#'@note Default \code{rand_tensor()} generates a 3-Tensor with modes \code{c(3,4,5)}.
#'@examples
#'rand_tensor()
#'rand_tensor(c(4,4,4))
#'rand_tensor(c(10,2,1),TRUE)
rand_tensor <- function(modes=c(3,4,5),drop=FALSE){
  as.tensor(array(rnorm(prod(modes)), dim=modes),drop=drop)
}



###Creation of Tensor from an array/matrix/vector

#'Tensor Conversion
#'
#'Create a \code{\link{Tensor-class}} object from an \code{array}, \code{matrix}, or \code{vector}.
#'@export
#'@name as.tensor
#'@rdname as.tensor
#'@aliases as.tensor
#'@param x an instance of \code{array}, \code{matrix}, or \code{vector}
#'@param drop whether or not modes of 1 should be dropped
#'@return a \code{\link{Tensor-class}} object
#'@examples
#'#From vector
#'vec <- runif(100); vecT <- as.tensor(vec); vecT
#'#From matrix
#'mat <- matrix(runif(1000),nrow=100,ncol=10)
#'matT <- as.tensor(mat); matT
#'#From array
#'indices <- c(10,20,30,40)
#'arr <- array(runif(prod(indices)), dim = indices)
#'arrT <- as.tensor(arr); arrT
as.tensor <- function(x,drop=FALSE){
  stopifnot(is.array(x)||is.vector(x))
  if (is.vector(x)){
    modes <- c(length(x))
    num_modes <- 1L
    new("Tensor", num_modes, modes, data = x)
  }
  else {
    modes <- dim(x)
    num_modes <- length(modes)
    dim1s <- which(modes==1)
    if (drop && (length(dim1s)>0)){
      modes <- modes[-dim1s]
      num_modes <- num_modes-length(dim1s)
      new("Tensor",num_modes,modes,data=array(x,dim=modes))
    }
    else{
      new("Tensor",num_modes,modes,data=x)
    }
  }
}

