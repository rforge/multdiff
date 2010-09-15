
### Begin: MultDiffModel

setClass(
         Class = "MultDiffModel",
         representation = representation(
           data = "ContinuousProcess",
           parameters = "list",
           "VIRTUAL"
           )
         )


setMethod(
          "getData",
          "MultDiffModel",
          function(object){
            return(object@data)
          }
          )


setMethod(
          "getParameters",
          "MultDiffModel",
          function(object){
            return(object@parameters)
          }
          )


### End: MultDiffModel-class

### Begin: Loss function

setMethod(
          "loss",
          "MultDiffModel",
          function(object, parameters, lossType){
            data <- as.matrix(getValue(object@data))
            if (length(data)==0)
              stop("'object' must contain 'data'")
            if (missing(parameters))
              {parameters <- object@parameters}
            #validateParameters(object,parameters)
            pos <- getPosition(object@data)
            n <- length(pos)
            p <- dim(parameters$A)[1]
            delta <- pos[2:n]-pos[1:(n-1)]
            if (lossType == 1)
              {
                tmpMean <- do.call(rbind,condMeanVar(object, parameters, x=data[1:(n-1),], t=delta)$condMean)
                return(sum((data[2:n,]-tmpMean)^2)/2)
              }
            if (lossType == 2)
              {
                tmpMeanVar <- condMeanVar(object, parameters, x=data[1:(n-1),], t=delta,var=TRUE)
                tmpMean <- do.call(rbind,tmpMeanVar$condMean)
                tmpVar <- t(sapply(tmpMeanVar$condVar,function(y){return(y)}))
                centeredObs <- data[2:n,]-tmpMean
                wCenteredObs <- t(apply(cbind(tmpVar,centeredObs),1,function(y){
                  solve(matrix(y[1:(p*p)],nrow=p,ncol=p),matrix(y[(p*p+1):(p*(p+1))]))
                }
                                        )
                                  )
                firstSum <- sum(
                                apply(cbind(centeredObs,wCenteredObs),1,function(y){
                                  matrix(y[1:p],nrow=1)%*%matrix(y[(p+1):(2*p)])
                                }
                                      )
                                )
                secondSum <- sum(apply(tmpVar,1,function(y){
                  log(det(matrix(y,nrow=p,ncol=p)))
                }
                                       )
                                 )
                return((firstSum + secondSum)/2)
              }
          }
          )
            
            


### end: Loss function


#setMethod("fitMultDiffModel",
#          "ContinuousProcess",
#          function(object,modelClass,lossType){
#            if(modelClass == "OUModel" && lossType == 1){
#              model <- new("OUModel", data=object)
#              p <- dim(getValue(object))[2]
#              
#              
#              l <- function(x){
#                loss(model,parToList(model,c(x,as.numeric(diag(1,p)))),lossType)+lambda*sum(abs(x))
#              }
#
#              L <- function(y){
#                loss(model,parToList(model,c(x,as.numeric(diag
#
#              optim(0,
#                optim(rep(1,p+p*p),fn=l)
#              
#                
#              }}}
#          )
              
