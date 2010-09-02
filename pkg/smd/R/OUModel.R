### Begin: OUModel-class

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

### End: OUModel-class


### Begin: condMeanVar function

setMethod(
          "condMeanVar",
          "OUModel",
          function(object,parameters,x,t,var=FALSE){
            if (missing(parameters))
              {parameters <- object@parameters}
            p <- dim(parameters$A)[1]
            if (!is.matrix(x))
              stop("'x' must be a matrix of class 'matrix'")
            if (dim(x)[2] != p)
              stop("Mismatch in dimensions of 'x' and parameters of 'object'")
            if (dim(x)[1] != length(t))
              stop("Mismatch in dimensions of 'x' and 't'")
            if (!is.logical(var))
              stop("'var' must be logical")
            if (var==FALSE)
              {
                tmpMean <- apply(cbind(x,t),1,function(y){
                  parameters$A + expm(y[p+1]*parameters$B)%*%(matrix(y[1:p])-parameters$A)
                }
                    )
                return(list(condMean = unname(split(tmpMean, col(tmpMean)))))
              }
            if (var==TRUE)
              {
                n <- length(t)
                f <- factor(t)
                tUnique <- sort(unique(t))
                X <- matrix(,nrow = 2*p+1, ncol = 2*p+1)
                X[1:p,] <- cbind(-parameters$B,-parameters$B%*%parameters$A,parameters$C)
                X[p+1,] <- rep(0,2*p+1)
                X[(p+2):(2*p+1),] <- cbind(matrix(rep(0,p*(p+1)),nrow=p,ncol=p+1),t(parameters$B))
                tmpFunctions <- lapply(tUnique,function(y)
                                       {exptX <- expm(y*X)
                                        G1 <- exptX[1:p,(p+1),drop=FALSE]
                                        H1 <- exptX[1:p,(p+2):(2*p+1)]
                                        F3 <- exptX[(p+2):(2*p+1),(p+2):(2*p+1)]
                                        tmp <- list(
                                                    alpha = t(F3)%*%G1,
                                                    beta = t(F3),
                                                    Sigma = t(F3)%*%H1
                                                    )
                                        return(tmp)
                                        }
                                       )
                tmpMean <- apply(cbind(x,f),1,function(y)
                                  {
                                    tmpFunctions[[y[p+1]]]$alpha + tmpFunctions[[y[p+1]]]$beta %*% matrix(y[1:p])
                                  }
                                  )
                condVar <- lapply(f,function(y)
                                  {
                                    tmpFunctions[[y]]$Sigma
                                  }
                                  )

                return(list(
                            condMean = unname(split(tmpMean, col(tmpMean))),
                            condVar = condVar
                            ))
              }
          }
          )      

### End: condMeanVar function



### Begin: Functions to do with parameters
           
            
setMethod("validateParameters",
          "OUModel",
          function(object,parameters){
            if (length(parameters)!=0){
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
              return(TRUE)
            }
          }
          )


setMethod("parToNumeric",
          "OUModel",
          function(object,parameters){
            if (missing(parameters))
              {parameters <- object@parameters}
            validateParameters(object,parameters)
            return(as.numeric(cbind(parameters$A,parameters$B,parameters$C)))
          }
          )

setMethod("parToList",
          "OUModel",
          function(object,x){
            L <- length(x)
            p <- (-1+sqrt(1+8*L))/4
            if (abs(p-round(p)!=0))
              stop("Error in length of 'x'")
            parameters <- list(
                        A = matrix(x[1:p]),
                        B = matrix(x[p+1:(p+p*p)],nrow=p,ncol=p),
                        C = matrix(x[(p+p*p+1):(2*p*p+p)], ncol=p,nrow=p)
                        )
            validateParameters(object,parameters)
            return(parameters)
          }
          )       
            

### End: functions to do with parameters

