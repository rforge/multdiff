
setClass(
         Class = "OUModel",
         representation = representation(
           parameters = "list"),
         contains = "MultDiffModel",
         validity = function(object){
           par <- object@parameters
           if ((length(par)!=3) | (sum(c("A","B","C") %in% names(par))!=3))
             stop("'parameters' must be a list with components 'A', 'B' and 'C'")
           if (!is.matrix(par$A))
             stop("'A' must be a matrix")
           if (!is.matrix(par$B))
             stop("'B' must be a matrix")
           if (!is.matrix(par$C))
             stop("'C' must be a matrix")
           dimdata <- dim(getValue(object@data))[2]
           if ((dim(par$A)[1] != dimdata) | (dim(par$A)[2] != 1))
             stop("Mismatch in dimensions of 'data' and 'A'")
           if ((dim(par$B)[1] != dimdata) | (dim(par$B)[2] != dimdata))
             stop("Mismatch in dimensions of 'data' and 'B'")
           if ((dim(par$C)[1] != dimdata) | (dim(par$C)[2] != dimdata))
             stop("Mismatch in dimensions of 'data' and 'C'")
           if (!identical(par$C,t(par$C)))
             stop("'C' must be symmetric")
           if (!isTRUE(all.equal(min(eigen(par$C,only.values = TRUE, symmetric = TRUE)$values,0),0)))
             stop("'C' must be positive semidefinite")
         }
         )

