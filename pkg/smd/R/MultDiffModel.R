
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

setMethod(
          "loss",
          "MultDiffModel",
          function(object, type){
            data <- as.matrix(getValue(object@data))
            if (length(data)==0)
              stop("'object' must contain data")
            parameters <- object@parameters
            pos <- getPosition(object@data)
            n <- length(pos)
            p <- dim(parameters$A)[1]
            delta <- pos[2:n]-pos[1:(n-1)]

            tmpData <- data[1:(n-1),]
            tmp <- condMeanVar(object=object,x=tmpData,t=delta)
            tmpMean <- t(sapply(tmp$condMean,function(y){return(y)}))
            tmpVar <- t(sapply(tmp$condVar,function(y){return(y)}))
            centeredObs <- data[2:n,]-tmpMean
            if (type == 1)
              {
                return(sum((centeredObs)^2)/2)
              }
           if (type == 2)
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


#setMethod(
#          "loss",
#          "MultDiffModel",
#          function(object, type){
#            data <- as.matrix(getValue(object@data))
#            if (length(data)==0)
#              stop("'object' must contain data")
#            parameters <- object@parameters
#            pos <- getPosition(object@data)
#            n <- length(pos)
#            p <- dim(parameters$A)[1]
#            delta <- pos[2:n]-pos[1:(n-1)]
#
#            tmpData <- cbind(data[1:(n-1),],delta)
#            tmp <- apply(tmpData,1,function(y)
#                         {
#                           condMeanVar(object,y[1:p],y[p+1])
#                         }
#                         )
#              
#            tmpMean <- t(sapply(tmp,function(y){y$condMean}))
#            tmpVar <- sapply(tmp,function(y){y$condVar})
                        #condMean <- t(sapply(1:n,function(y)
            #                     {
            #                       tmp <- tmpFunctions[[f[y]]]
            #                       tmp <- tmp$alpha + tmp$beta %*% t(x[y,,drop=FALSE])
            #                       return(tmp)
            #                     }
            #                     )
            #              )
            #tmpVar <- matrix(sapply(tmp,function(x){x$condVar})[,2],nrow=4)
 #           centObs <- data[2:n,]-tmpMean

#
#            if (type = "type1"){
#              print(sum((centObs)^2)/2)
#            }
#
#            if (type = "type2"){
#              diff
#   }
#            
#          }
#          )
