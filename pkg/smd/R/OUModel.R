
setClass(
         Class = "OUModel",
         representation = representation(
           A = "matrix",
           B = "matrix",
           D = "matrix"),
         contains = "MultDiffModel",
         validity = function(object){
           p <- dim(getValue(object@data))[2]
           if ((dim(object@A)[1] != p) || (dim(object@A)[2] != 1))
             stop("Mismatch in dimensions of 'data' and 'A'")
           if ((dim(object@B)[1] != p) || (dim(object@B)[2] != p))
             stop("Mismatch in dimensions of 'data' and 'B'")
           if ((dim(object@D)[1] != p) || (dim(object@D)[2] != p))
             stop("Mismatch in dimensions of 'data' and 'D'")
           # Check whether 'D' positive semi-definite - how?
         }
         )




      
         
