
setClass(
         Class = "OUModel",
         contains = "MultDiffModel",
         validity = function(object){
           par <- object@parameters
           if(length(par)!=0) {
             if ( (length(par)!=3) || (sum(c("A","B","C") %in% names(par))!=3) )
               stop("'parameters' must be a list with three components 'A', 'B' and 'C'")
             if ( (!is.matrix(par$A)) || (!is.matrix(par$B)) || (!is.matrix(par$C)) )
               stop("'A', 'B' and 'C' must all be of class 'matrix'")
             if (dim(par$A)[2] != 1)
               stop("'A' must be a one-column matrix'")
             if (dim(par$B)[1] != dim(par$B)[2])
               stop("'B' must be a square matrix'")
             if (dim(par$C)[1] != dim(par$C)[2])
               stop("'C' must be a square matrix'")
             if ((dim(par$A)[1] != dim(par$B)[1]) || (dim(par$A)[1] != dim(par$C)[1]))
               stop("Mismatch in dimensions of 'A', 'B' and 'C'")
             if (!identical(par$C,t(par$C)))
               stop("'C' must be symmetric")
             if (!isTRUE(all.equal(min(eigen(par$C,only.values = TRUE, symmetric = TRUE)$values,0),0)))
               stop("'C' must be positive semidefinite")
             p <- dim(getValue(object@data))[2]
             if (p!=0) {
               if (dim(par$A)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'A'")
               if (dim(par$B)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'B'")
               if (dim(par$C)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'C'")
             }
           }
         }
       )



setMethod(
          "condMeanVar",
          "OUModel",
          function(object,x,t){
            par <- object@parameters
            p <- dim(par$A)[1]
            if ((!is.matrix(x)) || (dim(x)[2] != 1))
              stop("'x' must be a one-column matrix of class 'matrix'")
            if (dim(x)[1] != p)
              stop("Mismatch in dimensions of 'x' and parameters of 'object'")
            if ((!is.numeric(t)) || (length(t)!=1) || (t<0))
              stop("'t' must be a non-negative real value")
            X <- matrix(,nrow = 2*p+1, ncol=2*p+1)
            X[1:p,] <- cbind(-par$B,par$A,par$C)
            X[p+1,] <- rep(0,2*p+1)
            X[(p+2):(2*p+1),] <- cbind(matrix(rep(0,p*(p+1)),nrow=p,ncol=p+1),t(par$B))
            exptX <- expm(t*X)
            alpha <- t(exptX[(p+2):(2*p+1),(p+2):(2*p+1)])%*%exptX[1:p,(p+1),drop=F]
            beta <- t(exptX[(p+2):(2*p+1),(p+2):(2*p+1)])
            Sigma <- t(exptX[(p+2):(2*p+1),(p+2):(2*p+1)])%*%exptX[1:p,(p+2):(2*p+1)]
            tmp <- list(x = x,
                        t = t,
                        condMean = alpha + beta %*% x,
                        condVar = Sigma
                        )
            return(tmp)

          }
          ) 





