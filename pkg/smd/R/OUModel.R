
setClass(
         Class = "OUModel",
         contains = "MultDiffModel",
         validity = function(object){
           parameters <- object@parameters
           if(length(parameters)!=0) {
             if ( (length(parameters)!=3) || (sum(c("A","B","C") %in% names(parameters))!=3) )
               stop("'parameters' must be a list with three components 'A', 'B' and 'C'")
             if ( (!is.matrix(parameters$A)) || (!is.matrix(parameters$B)) || (!is.matrix(parameters$C)) )
               stop("'A', 'B' and 'C' must all be of class 'matrix'")
             if (dim(parameters$A)[2] != 1)
               stop("'A' must be a one-column matrix'")
             if (dim(parameters$B)[1] != dim(parameters$B)[2])
               stop("'B' must be a square matrix'")
             if (dim(parameters$C)[1] != dim(parameters$C)[2])
               stop("'C' must be a square matrix'")
             if ((dim(parameters$A)[1] != dim(parameters$B)[1]) || (dim(parameters$A)[1] != dim(parameters$C)[1]))
               stop("Mismatch in dimensions of 'A', 'B' and 'C'")
             if (!identical(parameters$C,t(parameters$C)))
               stop("'C' must be symmetric")
             if (!isTRUE(all.equal(min(eigen(parameters$C,only.values = TRUE, symmetric = TRUE)$values,0),0)))
               stop("'C' must be positive semidefinite")
             p <- dim(getValue(object@data))[2]
             if (p!=0) {
               if (dim(parameters$A)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'A'")
               if (dim(parameters$B)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'B'")
               if (dim(parameters$C)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'C'")
             }
           }
         }
       )





setMethod(
          "condMeanVarFct",
          "OUModel",
          function(object,t){
            parameters <- object@parameters
            p <- dim(parameters$A)[1]
            if ((!is.numeric(t)) || (length(t)!=1) || (t<0))
              stop("'t' must be a non-negative real value")
            X <- matrix(,nrow = 2*p+1, ncol=2*p+1)
            X[1:p,] <- cbind(-parameters$B,parameters$A,parameters$C)
            X[p+1,] <- rep(0,2*p+1)
            X[(p+2):(2*p+1),] <- cbind(matrix(rep(0,p*(p+1)),nrow=p,ncol=p+1),t(parameters$B))
            exptX <- expm(t*X)
            G1 <- exptX[1:p,(p+1),drop=FALSE]
            H1 <- exptX[1:p,(p+2):(2*p+1)]
            F3 <- exptX[(p+2):(2*p+1),(p+2):(2*p+1)]
            alpha <- t(F3)%*%G1
            beta <- t(F3)
            Sigma <- t(F3)%*%H1
            tmp <- list(
                        alpha = alpha,
                        beta = beta,
                        Sigma = Sigma
                        )
            return(tmp)
          }
          )

setMethod(
          "condMeanVar",
          "OUModel",
          function(object,x,t){
            parameters <- object@parameters
            p <- dim(parameters$A)[1]
            if (!is.matrix(x))
              stop("'x' must be a matrix of class 'matrix'")
            if (dim(x)[2] != p)
              stop("Mismatch in dimensions of 'x' and parameters of 'object'")
            if (dim(x)[1] != length(t))
              stop("Mismatch in dimensions of 'x' and 't'")
            n <- length(t)
            f <- factor(t)
            tmpData <- cbind(x,t,f)
            tUnique <- sort(unique(t))
            tmpFunctions <- lapply(tUnique,function(y){condMeanVarFct(object,y)})
            condMean <- lapply(1:n,function(y)
                                 {
                                   tmp <- tmpFunctions[[f[y]]]$alpha + tmpFunctions[[f[y]]]$beta %*% t(x[y,,drop=FALSE])
                                   return(tmp)
                                 }
                                 )
            condVar <- lapply(1:n,function(y)
                              {
                                tmp <- tmpFunctions[[f[y]]]$Sigma
                                return(tmp)
                              }
                              )
            tmp <- list(
                        condMean = condMean,
                        condVar = condVar
                        )
            return(tmp)
          }
          )


### Works if x is vector, t is non-negative real number. Returns one condMean, one condVar

#setMethod(
#          "condMeanVar",
#          "OUModel",
#          function(object,x,t){
#            parameters <- object@parameters
#            p <- dim(parameters$A)[1]
#            x <- as.matrix(x)
#            if (length(x) != p)
#              stop("Mismatch in dimensions of 'x' and parameters of 'object'")
#            tmp <- condMeanVarFct(object,t)
#            tmp <- list(
#                        condMean = tmp$alpha + tmp$beta %*% x,
#                        condVar = tmp$Sigma
#                        )
#            return(tmp)
#          }
#          )


### Works if x is vector, t is non-negative real number. Returns one condMean, one condVar. Independent of condMeanVarFct

#setMethod(
#          "condMeanVar",
#          "OUModel",
#          function(object,x,t){
#            parameters <- object@parameters
#            p <- dim(parameters$A)[1]
#            x <- as.matrix(x)
#            if (length(x) != p)
#              stop("Mismatch in dimensions of 'x' and parameters of 'object'")
#            if ((!is.numeric(t)) || (length(t)!=1) || (t<0))
#              stop("'t' must be a non-negative real value")
#            X <- matrix(,nrow = 2*p+1, ncol=2*p+1)
#            X[1:p,] <- cbind(-parameters$B,parameters$A,parameters$C)
#            X[p+1,] <- rep(0,2*p+1)
#            X[(p+2):(2*p+1),] <- cbind(matrix(rep(0,p*(p+1)),nrow=p,ncol=p+1),t(parameters$B))
#            exptX <- expm(t*X)
#            G1 <- exptX[1:p,(p+1),drop=FALSE]
#            H1 <- exptX[1:p,(p+2):(2*p+1)]
#            F3 <- exptX[(p+2):(2*p+1),(p+2):(2*p+1)]
#            alpha <- t(F3)%*%G1
#            beta <- t(F3)
#            Sigma <- t(F3)%*%H1
#            tmp <- list(
#                        condMean = alpha + beta %*% x,
#                        condVar = Sigma
#                        )
#            return(tmp)
#
#          }
#          )
