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
             if ( !is.vector(A) )
               stop("'A' must be a vector")
             if ( (!is.matrix(parameters$B)) || (!is.matrix(parameters$C)) )
               stop("'B' and 'C' must all be of class 'matrix'")
             if (dim(parameters$B)[1] != dim(parameters$B)[2])
               stop("'B' must be a square matrix'")
             if (dim(parameters$C)[1] != dim(parameters$C)[2])
               stop("'C' must be a square matrix'")
             if ((length(parameters$A) != dim(parameters$B)[1]) || (length(parameters$A) != dim(parameters$C)[1]))
               stop("Mismatch in dimensions of 'A', 'B' and 'C'")
             if (!identical(parameters$C, t(parameters$C)))
               stop("'C' must be symmetric")
             if (!isTRUE(all.equal(min(eigen(parameters$C,only.values = TRUE, symmetric = TRUE)$values,0),0)))
               stop("'C' must be positive semidefinite")
             p <- dim(getValue(object@data))[2]
             if (p!=0) {
               if (length(parameters$A) != p)
                 stop("Mismatch in dimensions of 'data' and 'A'")
               if (dim(parameters$B)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'B'")
               if (dim(parameters$C)[1] != p)
                 stop("Mismatch in dimensions of 'data' and 'C'")
             }
           }
         }
       )

### Method: constructor

setMethod(
          "ouModel",
          c("ContinuousProcess"),
          function(data, A, B, C, computeSufficient = TRUE, ...) {
            
            .object <- new("OUModel",
                           data = data,
                           parameters = list(
                             A = A,
                             B = B,
                             C = C),
                           hasSufficientStat = FALSE
                           )

            if(computeSufficient)
              .object <- computeSufficientStat(.object)

            return(.object)            
          }
          )

setMethod(
          "computeSufficientStat",
          "OUModel",
          function(object, ...) {
            
            if(getEquiDistance(object@data) > 0) {
              sufficientStat = new.env(parent = .GlobalEnv)
              X <- getNumerics(object@data)
              n <- dim(X)[1]
              assign("normObs", sum(X[-1, ]^2), envir = sufficientStat)
              Sm1 <- colSums(X[-1, ])
              assign("Sm1", Sm1, envir = sufficientStat)
              Smn <- Sm1 + X[1, ] - X[n, ]
              assign("Smn", Smn, envir = sufficientStat)              
              assign("SS", crossprod(X[-n, ]), envir = sufficientStat)
              assign("SsS", crossprod(X[-1, ], X[-n, ]), envir = sufficientStat)
              object@sufficientStat <- sufficientStat
              object@hasSufficientStat <- TRUE
            } else {
              cat("Currently, computations of sufficient statistics for non-equidistant observations are not supported.")
            }
            return(object)
          }
          )

setMethod(
          "hasSufficientStat",
          "OUModel",
          function(object, ...) {
            return(object@hasSufficientStat)
          }
          )

setMethod(
          "loss",
          "OUModel",
          function(object, parameters = NULL, lossType = 1, useSufficientStat = TRUE, ...) {
            Delta <- getEquiDistance(object@data)
            useSufficientStat <- hasSufficientStat(object) & useSufficientStat
            if(lossType == 2 || !useSufficientStat || Delta == 0)
              return(callNextMethod(object = as(object, "MultDiffModel"),
                                 parameters = parameters,
                                 lossType = lossType))

            if(is.null(parameters))
              parameters <- getParameters(object)

            normObs <- getSufficientStat(object, "normObs")
            Sm1 <- getSufficientStat(object, "Sm1")
            Smn <- getSufficientStat(object, "Smn")
            SS <- getSufficientStat(object, "SS")
            SsS <- getSufficientStat(object, "SsS")
            n <- dim(object@data)[1]
            tmpExp <- expm(Delta*parameters$B)
            tmpExpA <- parameters$A - tmpExp %*% parameters$A
            l <- normObs +
              sum(crossprod(tmpExp) * SS) +
              (n-1)*sum(tmpExpA^2) +
              2*((tcrossprod(Smn, tmpExp) - Sm1) %*% tmpExpA -
                 sum(tmpExp * SsS))
            
            return(l)
          }
          )

### Method: condMeanVar
### Computes the conditional mean and variance of
### X_{t_0+t} given X_{t_0} = x.

setMethod(
          "condMeanVar",
          "OUModel",
          function(object, parameters, x, t, var = FALSE, ...){
            if (missing(parameters)) 
              parameters <- getParameters(object)
            if (length(parameters)==0)
              stop("'object' must contain 'parameters' or 'parameters' must be specified")
            p <- length(parameters$A)
            
            if (missing(t) & !missing(x))
              stop("If 'x' is specified, 't' must be specified as well")
            if (!missing(t) & missing(x))
              stop("If 't' is specified, 'x' must be specified as well")
            
            if (missing(t) & missing(x)){
              t <- getEquiDistance(object@data)
              if (t==0 || length(t)==0)
                t <- diff(getTime(object@data))
              m <- length(t)
              if (m==0)
                stop("Object must contain data or 'x' and 't' must be specified.")
              tmpx <- as.matrix(getValue(object@data))
              N <- dim(tmpx)[1]-1
              x <- tmpx[1:N, ]
            } else {
              m <- length(t)
              N <- dim(x)[1]
              if (!is.matrix(x))
                stop("'x' must be a matrix of class 'matrix'.")
              if (!(m %in% c(1,N)))
                stop("Mismatch in dimensions of 'x' and 't'.")
            }

            if (p != dim(x)[2])
              stop("Mismatch in dimensions of 'x' and parameters of 'object'.")
            if (!is.logical(var))
              stop("'var' must be logical.")
            
            if (!var) # The conditional variance is not computed
              {
                tmpMean <-  array(, dim=c(p,N))
                if (m == 1){
                  tmpExp <- expm(t*parameters$B)
                  ## TODO: Reformulate as matrix computations
                  for (i in 1:N){
                    tmpMean[ , i] <- parameters$A + tmpExp %*% (x[i, ]-parameters$A)
                  }
                } else {
                  for (i in 1:N){
                    tmpMean[ , i] <- parameters$A + expm(t[i]*parameters$B) %*% (x[i, ]-parameters$A)
                  }
                }
                return(tmpMean)
              } 
            if (var) # The conditional variance is computed
              {
                X <- matrix(0, nrow = 2*p+1, ncol = 2*p+1)
                X[1:p, 1:p] <- -parameters$B
                X[1:p, p+1] <- -parameters$B %*% parameters$A
                X[1:p, (p+2):(2*p+1)] <- parameters$C
                X[p+1, ] <- rep(0,2*p+1)
                X[(p+2):(2*p+1), 1:(p+1)] <- matrix(rep(0,p*(p+1)), nrow = p, ncol = p+1)
                X[(p+2):(2*p+1), (p+2):(2*p+1)] <- t(parameters$B)

                tmpMean <- array( , dim=c(p,N))
                tmpVar <- array( , dim=c(p,p,N))
                
                if (m == 1){
                  exptX <- expm(t*X)
                  G <- exptX[1:p, (p+1), drop=FALSE]
                  H <- exptX[1:p, (p+2):(2*p+1)]
                  tF <- t(exptX[(p+2):(2*p+1), (p+2):(2*p+1)])

                  for (i in 1:N){
                    tmpMean[,i] <- tF %*% G + tF %*% x[i, ]
                    tmpVar[,,i] <- tF %*% H
                  }
                  return(list(condMean = tmpMean, condVar = tmpVar))
                  } else {
                    for (i in 1:N){
                      exptX <- expm(t[i]*X)
                      G <- exptX[1:p, (p+1), drop=FALSE]
                      H <- exptX[1:p, (p+2):(2*p+1)]
                      tF <- t(exptX[(p+2):(2*p+1), (p+2):(2*p+1)])
                      tmpMean[ , i] <- tF %*% G + tF %*% x[i, ]
                      tmpVar[ , , i] <- tF %*% H
                    }
                    return(list(condMean = tmpMean, condVar = tmpVar))
                  }
              }
          }
          )


setMethod("gradient",
          "OUModel",
          function(object, parameters = NULL, lossType = 1, useSufficientStat = TRUE, ...){

            if (is.null(parameters))
              parameters <- object@parameters
            Delta <- getEquiDistance(object@data)
            
            p <- length(parameters$A)
            if (p != dim(object@data)[2])
              stop("Mismatch in dimensions of object data and 'parameters'")
            if (lossType==1){
              if(useSufficientStat & hasSufficientStat(object) & Delta != 0) {
                Sm1 <- getSufficientStat(object, "Sm1")
                Smn <- getSufficientStat(object, "Smn")
                SS <- getSufficientStat(object, "SS")
                SsS <- getSufficientStat(object, "SsS")
                n <- dim(object@data)[1]
                Delta <- getEquiDistance(object@data)
                DeltaB <- Delta*parameters$B
                tmpExp <- expm(DeltaB)
                tmpExpA <- parameters$A - tmpExp %*% parameters$A
                E <- t(SsS + tmpExpA %*% ((n-1)*t(parameters$A) - Smn) +
                       tmpExp %*% t(parameters$A %*% t(Smn) - SS)) - parameters$A %*% t(Sm1)
                B <- -2*t(expmFrechet(DeltaB, Delta*E, expm = FALSE)$Lexpm)
                A <- Sm1 - (n-1)*tmpExpA - tmpExp %*% Smn
                A <- 2* t(t(A) %*% tmpExp - t(A))
              } else {  ## Not using sufficient stat
                data <- getNumerics(object@data)
                if (length(data)==0)
                  stop("'object' must contain 'data'.")
                n <- dim(data)[1]
                if (Delta == 0 || length(Delta) == 0)
                  Delta <- diff(getTime(object@data))
                m <- length(Delta)
                tmpMean <- condMeanVar(object, parameters)
                centeredObs <- t(data[2:n, ]) - tmpMean
                if (m == 1){
                  tmp1 <- vector("list", n-1)
                  for (i in 1:(n-1)){
                    tmp1[[i]] <- (t(data[i,,drop=F])-parameters$A) %*% t(centeredObs[,i,drop=F])
                  }
                  E <- unname(Reduce("+", tmp1))
                  tmp2 <- Reduce("+", lapply(1:(n-1),function(x){centeredObs[,x,drop=F]}))
                  B <- -2*t(expmFrechet(Delta*parameters$B, Delta*E)$Lexpm)
                  A <- -2*(diag(1,p)-expm(Delta*t(parameters$B)))%*%tmp2
                  return(list(A=A,B=B))
                } else {
                  tmp1 <- matrix(,nrow=p*p,ncol=n-1)
                  tmp2 <- matrix(,nrow=p,ncol=n-1)
                  for (i in 1:(n-1)){
                    E <- (matrix(data[i,])-parameters$A)%*%matrix(centeredObs[,i], nrow=1)
                    tmp1[,i] <- as.numeric(-(expmFrechet(Delta[i]*parameters$B, Delta[i]*E)$Lexpm))
                    tmp2[,i] <- -(diag(1,p)-expm(Delta[i]*t(parameters$B)))%*%matrix(centeredObs[,i])
                  }
                  B <- 2*t(matrix(tmp1%*%matrix(rep(1,n-1)),nrow=p,ncol=p))
                  A <- 2*matrix(tmp2%*%matrix(rep(1,n-1)))
                }
              }
            }
            if (lossType==2){
              stop("Not implemented yet")
            }
            return(list(A=A,B=B))
          }
          )
                             

setMethod("validateParameters",
          "OUModel",
          function(object, parameters){
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
          function(object, parameters, ...){
            if (missing(parameters))
              parameters <- object@parameters
            #validateParameters(object,parameters)
            return(unlist(parameters))
          }
          )

setMethod("parToList",
          "OUModel",
          function(object, x, ...){
            L <- length(x)
            p <- (-1+sqrt(1+8*L))/4
            if (abs(p-round(p)) != 0)
              stop("Error in length of 'x'")
            parameters <- list(
                        A = x[1:p],
                        B = matrix(x[p+1:(p+p*p)], nrow=p, ncol=p),
                        C = matrix(x[(p+p*p+1):(2*p*p+p)], ncol=p, nrow=p)
                        )
            #validateParameters(object,parameters)
            return(parameters)
          }
          )       
            

