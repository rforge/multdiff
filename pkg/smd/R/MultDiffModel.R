# (Virtual) class: MultDiffModel
# Represents a multivariate diffusion model

setClass(
         Class = "MultDiffModel",
         representation = representation(
           data = "ContinuousProcess",
           parameters = "list",
           sufficientStat = "environment",
           hasSufficientStat = "logical",   ## Flag to indicate if sufficient statistics are computed 
           "VIRTUAL"
           )
         )

# Method: getData
# Retrieves data from an object inheriting from class MultDiffModel

setMethod(
          "getData",
          "MultDiffModel",
          function(object, ...){
            return(object@data)
          }
          )

# Method: getParameters
# Retrieves parameters from an object inheriting from class MultDiffModel


setMethod(
          "getParameters",
          "MultDiffModel",
          function(object, ...){
            return(object@parameters)
          }
          )


# Method: loss
# Generic calculations of the value of the type 1 and 2 loss functions

setMethod(
          "loss",
          "MultDiffModel",
          function(object, parameters = NULL, lossType = 1, ...){
            data <- getNumerics(object@data)
            if (length(data)==0)
              stop("'object' must contain 'data'")
            n <- dim(data)[1]
            if (is.null(parameters))
              parameters <- object@parameters
            ## Not meaningfull to make general test of dimensions here
            ## p <- dim(parameters$A)[1]
            ## if (p != dim(data)[2])
            ##  stop("Mismatch in dimensions of object data and 'parameters'")
            if (lossType == 1)
              {
                tmpMean <- condMeanVar(object, parameters)
                return(sum((t(data[-1, ])-tmpMean)^2))
              }
            if (lossType == 2)
              {
                tmpMeanVar <- condMeanVar(object, parameters, var=TRUE)
                tmpMean <- tmpMeanVar$condMean
                tmpVar <- tmpMeanVar$condVar
                centeredObs <- t(data[2:n, ]) - tmpMean
                tmp1 <- rep(0, n-1)
                tmp2 <- rep(0, n-1)
                for (i in 1:(n-1)){
                  tmp1[i] <- t(centeredObs[ , i]) %*% solve(tmpVar[ , , i], centeredObs[ , i])
                  tmp2[i] <- log(det(tmpVar[ , , i]))
                }
                return(sum(tmp1+tmp2))
              }
          }
          )

### Method: getSufficientStat
### Returns the requested sufficient statistics if it exists.

setMethod(
          "getSufficientStat",
          c("OUModel", "character"),
          function(object, stat, ...) {
            if(exists(stat, object@sufficientStat)) {
              return(get(stat, object@sufficientStat))
            } else {
              stop(paste("Sufficient statistics", stat, "does not exist."))
            }
          }
          )
