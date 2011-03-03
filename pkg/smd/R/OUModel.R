### Class: OUModel
### Represents a multivariate Ornstein Uhlenbeck Diffusion Model

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


### Method: condMeanVar
### Computes the conditional mean and variance of
### X_{t_0+t} given X_{t_0} = x.

setMethod(
          "condMeanVar",
          "OUModel",
          function(object, parameters, x, t, equitmp, var = FALSE){
            if (missing(parameters)) 
              parameters <- object@parameters
            if (length(parameters)==0)
              stop("'object' must contain 'parameters' or 'parameters' must be specified")
            p <- dim(parameters$A)[1]
            
            if (missing(t) & !missing(x))
              stop("If 'x' is specified, 't' must be specified as well")
            if (!missing(t) & missing(x))
              stop("If 't' is specified, 'x' must be specified as well")
            
            if (missing(t) & missing(x)){
              t <- equitmp
              #t <- getDistance(object@data)
              if (t==0 || length(t)==0)
                t <- diff(getTime(object@data))
              m <- length(t)
              if (m==0)
                stop("Object must contain data or 'x' and 't' must be specified")
              tmpx <- as.matrix(getValue(object@data))
              N <- dim(tmpx)[1]-1
              x <- tmpx[1:N,]
            } else {
              m <- length(t)
              N <- dim(x)[1]
              if (!is.matrix(x))
                stop("'x' must be a matrix of class 'matrix'")
              if (!(m %in% c(1,N)))
                stop("Mismatch in dimensions of 'x' and 't'")
            }

            if (p != dim(x)[2])
              stop("Mismatch in dimensions of 'x' and parameters of 'object'")
            if (!is.logical(var))
              stop("'var' must be logical")
            
            if (!var) # The conditional variance is not computed
              {
                tmpMean <-  array(, dim=c(p,N))
                if (m==1){
                  tmpExp <- expm(t*parameters$B)
                  for (i in 1:N){
                    tmpMean[,i] <- parameters$A + tmpExp%*%(matrix(x[i,])-parameters$A)
                  }
                } else {
                  for (i in 1:N){
                    tmpMean[,i] <- parameters$A + expm(t[i]*parameters$B)%*%(matrix(x[i,])-parameters$A)
                  }
                }
                return(tmpMean)
              } 
            if (var) # The conditional variance is computed
              {
                X <- matrix(,nrow=2*p+1, ncol=2*p+1)
                X[1:p,1:p] <- -parameters$B
                X[1:p,p+1] <- -parameters$B%*%parameters$A
                X[1:p,(p+2):(2*p+1)] <- parameters$C
                X[p+1,] <- rep(0,2*p+1)
                X[(p+2):(2*p+1),1:(p+1)] <- matrix(rep(0,p*(p+1)),nrow=p,ncol=p+1)
                X[(p+2):(2*p+1),(p+2):(2*p+1)] <- t(parameters$B)

                tmpMean <- array(, dim=c(p,N))
                tmpVar <- array(, dim=c(p,p,N))
                
                if (m==1){
                  exptX <- expm(t*X)
                  G <- exptX[1:p,(p+1),drop=FALSE]
                  H <- exptX[1:p,(p+2):(2*p+1)]
                  tF <- t(exptX[(p+2):(2*p+1),(p+2):(2*p+1)])

                  for (i in 1:N){
                    tmpMean[,i] <- tF%*%G + tF%*%matrix(x[i,])
                    tmpVar[,,i] <- tF%*%H
                  }
                  return(list(condMean = tmpMean, condVar = tmpVar))
                  } else {
                    for (i in 1:N){
                      exptX <- expm(t[i]*X)
                      G <- exptX[1:p,(p+1),drop=FALSE]
                      H <- exptX[1:p,(p+2):(2*p+1)]
                      tF <- t(exptX[(p+2):(2*p+1),(p+2):(2*p+1)])
                      tmpMean[,i] <- tF%*%G + tF%*%matrix(x[i,])
                      tmpVar[,,i] <- tF%*%H
                    }
                    return(list(condMean = tmpMean, condVar = tmpVar))
                  }
              }
          }
          )


setMethod("gradient",
          "OUModel",
          function(object, parameters, equitmp, lossType){
            data <- as.matrix(getValue(object@data))
            if (length(data)==0)
              stop("'object' must contain 'data'")
            n <- dim(data)[1]
            if (missing(parameters))
              {parameters <- object@parameters}
            p <- dim(parameters$A)[1]
            if (p != dim(data)[2])
              stop("Mismatch in dimensions of object data and 'parameters'")
            delta <- equitmp
            #delta <- getDistance(object@data)
            if (delta==0 || length(delta)==0)
              delta <- diff(getTime(object@data))
            m <- length(delta)
            tmpMean <- condMeanVar(object, parameters, equitmp=equitmp)
            centeredObs <- t(data[2:n,]) - tmpMean
            if (lossType==1){
              if (m==1){
                tmp1 <- vector("list",n-1)
                for (i in 1:(n-1)){
                  tmp1[[i]] <- (t(data[i,,drop=F])-parameters$A)%*%t(centeredObs[,i,drop=F])
                }
                E <- unname(Reduce("+", tmp1))
                tmp2 <- Reduce("+", lapply(1:(n-1),function(x){centeredObs[,x,drop=F]}))
                B <- -t(expmFrechet(delta*parameters$B,delta*E)$Lexpm)
                A <- -(diag(1,p)-expm(delta*t(parameters$B)))%*%tmp2
                return(list(A=A,B=B))
              } else {
                tmp1 <- matrix(,nrow=p*p,ncol=n-1)
                tmp2 <- matrix(,nrow=p,ncol=n-1)
                for (i in 1:(n-1)){
                  E <- (matrix(data[i,])-parameters$A)%*%matrix(centeredObs[,i], nrow=1)
                  tmp1[,i] <- as.numeric(-(expmFrechet(delta[i]*parameters$B, delta[i]*E)$Lexpm))
                  tmp2[,i] <- -(diag(1,p)-expm(delta[i]*t(parameters$B)))%*%matrix(centeredObs[,i])
                }
                B <- t(matrix(tmp1%*%matrix(rep(1,n-1)),nrow=p,ncol=p))
                A <- matrix(tmp2%*%matrix(rep(1,n-1)))
                return(list(A=A,B=B))
              }
            }
            if (lossType==2){
              stop("Not implemented yet")
            }
          }
          )
                             

setMethod("validateParameters",
          "OUModel",
          function(object,parameters){
            if (missing(parameters))
              {parameters <- object@parameters}
            if (length(parameters)==0)
              stop("'object' must contain 'parameters' or 'parameters' must be specified")
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
          )


setMethod("parToNumeric",
          "OUModel",
          function(object,parameters){
            if (missing(parameters))
              {parameters <- object@parameters}
            #validateParameters(object,parameters)
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
            #validateParameters(object,parameters)
            return(parameters)
          }
          )       
            

