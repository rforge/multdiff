
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

### Begin: Loss function (NB!!!! OLD!!!! Will be changed!!!)

setMethod(
          "loss",
          "MultDiffModel",
          function(object, parameters, lossType){
            data <- as.matrix(getValue(object@data))
            if (length(data)==0)
              stop("'object' must contain data")
            if (missing(parameters))
              {parameters <- object@parameters}
            pos <- getPosition(object@data)
            n <- length(pos)
            p <- dim(parameters$A)[1]
            delta <- pos[2:n]-pos[1:(n-1)]

            tmpData <- data[1:(n-1),]
            tmp <- condMeanVar(object=object,parameters=parameters,x=tmpData,t=delta)
            tmpMean <- t(sapply(tmp$condMean,function(y){return(y)}))
            tmpVar <- t(sapply(tmp$condVar,function(y){return(y)}))
            centeredObs <- data[2:n,]-tmpMean
            if (lossType == 1)
              {
                return(sum((centeredObs)^2)/2)
              }
           if (lossType == 2)
             {
               tmpData <- cbind(centeredObs,tmpVar)
               wCenteredObs <- t(apply(tmpData,1,function(y)
                                     {
                                       solve(matrix(y[(p+1):((p+1)*p)],nrow=p),matrix(y[1:p]))
                                     }
                                     )
                                 )
               tmpData <- cbind(centeredObs,wCenteredObs)
               theSum <- sum(
                             apply(tmpData,1,function(y)
                                   {
                                     matrix(y[1:p],nrow=1)%*%matrix(y[(p+1):(2*p)])
                                   }
                                   )
                             )
               logDet <- sum(apply(tmpVar,1,function(y)
                               {
                                 log(det(matrix(y,nrow=p)))
                               }
                               )
                             )
               return((theSum+logDet)/2)
               

             }
          }
          )

### End: Loss function


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
              
