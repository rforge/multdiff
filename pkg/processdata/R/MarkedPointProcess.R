setMethod("markedPointProcess", c("data.frame", "ContinuousProcess"),
          function(pointData, continuousData, markVar = 'markType', ...) {            

            ## Initialize an empty object - bypasses validity
            ## check. Validity function called before returning object
            ## at the end of the contructor.
            pointProcess <- new("MarkedPointProcess") 
            as(pointProcess, "ContinuousProcess") <- continuousData
            if(pointProcess@positionVar %in% names(pointData)) {
              position <- pointData[ ,pointProcess@positionVar]
            } else {
              stop(paste("pointData needs a column named", pointProcess@positionVar))
            }
            if(continuousData@idVar %in% names(pointData)) {
              id <- pointData[ ,pointProcess@idVar]
             } else {
              id <- rep(getId(pointProcess)[1], dim(pointData)[1])
            }
            ord <- tapply(position, id, order)
            ord <- unlist(lapply(levels(id), function(i) which(id == i)[ord[[i]]]), use.names=FALSE)
            pointData <-  pointData[ord, , drop = FALSE]
            id <- id[ord]
            position <- position[ord]
            
            if(markVar %in% names(pointData)) {
              markType <- pointData[ , markVar]
            } else {
              markType <-  rep("point", dim(pointData)[1])
            }
            markValue <- pointData[ ,!(names(pointData) %in%
                                       c(pointProcess@idVar,
                                         pointProcess@positionVar,
                                         markVar)),
                                   drop = FALSE]

            pointProcessEnv <- new.env(parent = .GlobalEnv)
            pointProcessEnv$id <- factor(id, levels = levels(getId(pointProcess)))
            pointProcessEnv$position <- as.numeric(position)
            pointProcessEnv$markType <- factor(markType)
            pointProcessEnv$markValue <- markValue
            
            pointProcess@pointProcessEnv <- pointProcessEnv
            pointProcess@markVar <- markVar
            iPointSubset(pointProcess) <- -1L
            jPointSubset(pointProcess) <- -1L
            colNames(pointProcess) <- c(colNames(pointProcess),
                                        levels(getMarkType(pointProcess)),
                                        names(getMarkValue(pointProcess)))

            lockEnvironment(pointProcess@pointProcessEnv, bindings = TRUE)

            
            ## Compute a combined set of evaluation points and keep track of which
            ## points are points from the marked point process.
            pointId <- getPointId(pointProcess)
            contId <- getId(pointProcess)
            pointPosition <- getPointPosition(pointProcess)
            contPosition <- getPosition(pointProcess)
            evalPosition <- data.frame(id = factor(c(as.character(pointId), as.character(contId)), levels = levels(contId)),
                                       position = c(pointPosition, contPosition),
                                       pointer = c(seq_along(pointPosition), rep(Inf, length(contPosition))))

            evalPosition <- split(evalPosition, evalPosition$id)
            ## Dropping duplicated entries within 'id'.
            dup <- lapply(evalPosition, function(e) duplicated(e$position))
            evalPosition <- sapply(names(evalPosition), function(i) evalPosition[[i]][!(dup[[i]] & evalPosition[[i]]$pointer == Inf), ], simplify = FALSE)

            ## Order within 'id'.
            ord <- lapply(evalPosition, function(e) order(e$position))
            evalPosition <- sapply(names(evalPosition), function(i) evalPosition[[i]][ord[[i]], ],  simplify = FALSE)

            contPosition <- split(data.frame(pointer = seq_along(contPosition), position = contPosition), contId)

            ii <- sapply(names(evalPosition),
                         function(i) {
                           findInterval(evalPosition[[i]]$position,
                                        c(contPosition[[i]]$position, Inf),
                                        all.inside=TRUE)
                         },
                         simplify = FALSE)

            ii <- unlist(lapply(names(evalPosition), function(i) contPosition[[i]]$pointer[ii[[i]]]), use.names = FALSE)

            evalPosition <- do.call("rbind", evalPosition)                 
           
            pointProcess@pointPointer <- which(evalPosition$pointer < Inf)
          
            valueEnv <- new.env(parent = .GlobalEnv)
            valueEnv$id <- evalPosition$id
            valueEnv$position <- evalPosition$position
            valueEnv$value <- getValue(pointProcess)[ii, , drop=FALSE]
            valueEnv$i <- seq_along(evalPosition$id)
            valueEnv$j <- seq_len(dim(valueEnv$value)[2])
            pointProcess@valueEnv <- valueEnv
            lockEnvironment(pointProcess@valueEnv, bindings = TRUE)
            validObject(pointProcess)
            return(pointProcess)
            
          }
          )

setMethod("markedPointProcess", c("data.frame", "data.frame"),
          function(pointData, continuousData, markVar = 'markType', ...) {
            callGeneric(pointData, continuousProcess(continuousData, ...),  markVar = markVar, ...)
          }
          )

setMethod("markedPointProcess", c("data.frame", "vector"),
          function(pointData, continuousData, positionVar = 'time', idVar = 'id', markVar = 'markType', ...) {
            if(!(idVar %in% names(pointData))) {
            continuousData <- data.frame(continuousData)
            names(continuousData) <- positionVar
            callGeneric(pointData = pointData,
                        continuousData = continuousData,
                        positionVar = positionVar, idVar = idVar, markVar = markVar, ...)
          } else {
            id <- factor(rep(levels(pointData[, idVar]), each = length(continuousData))) 
            continuousData <- rep(continuousData, length(levels(pointData[, idVar])))
            continuousData <- data.frame(id, continuousData)
            names(continuousData) <- c(idVar, positionVar)
            callGeneric(pointData = pointData,
                        continuousData = continuousData,
                        positionVar = positionVar, idVar = idVar, markVar = markVar, ...)
            }
          }
          )

setMethod("markedPointProcess", "data.frame",
          function(pointData, positionVar = 'time', idVar = 'id', markVar = 'markType',...) {
            if(!(positionVar %in% names(pointData)))
              stop(paste("pointData must have a column named", positionVar))
            if(!(idVar %in% names(pointData))) {
            callGeneric(pointData = pointData,
                        continuousData = pointData[ , positionVar, drop = FALSE],
                        positionVar = positionVar, idVar = idVar, markVar = markVar, ...)
          } else {
            callGeneric(pointData = pointData,
                        continuousData = pointData[ , c(idVar, positionVar)],
                        positionVar = positionVar, idVar = idVar, markVar = markVar, ...)
            }
            
          }
          )

setMethod("markedPointProcess", "vector",
          function(pointData, positionVar = 'time', idVar = 'id', markVar = 'markType',...) {
           pointData <- data.frame(pointData)
           names(pointData) <- positionVar
           callGeneric(pointData = pointData,
                       positionVar = positionVar, idVar = idVar, markVar = markVar, ...)         
          }
          )

setMethod("integrator", "MarkedPointProcess",
          function(object, f = 1, jumpVar = '', result = 'JumpProcess', ...)
          {

            id <- getPointId(object)
            jump <- NULL
            
            ## The counting process
            process <- rep(0,length(getId(object)))
            process[getPointPointer(object)] <- 1
            process <- tapply(process, getId(object), cumsum)
            
            if(is.function(f))
              f <- f(getPointPosition(object), ...)

            if(length(f) == 1 && f == 1 && jumpVar %in% colnames(getMarkValue(object))) {
              jump <- getMarkValue(object)[ , jumpVar]
            } else if(length(f) > 1 || any(f != 1)) {
              if(jumpVar %in% colnames(getMarkValue(object))) {
                jump <- f*getMarkValue(object)[ , jumpVar]
              } else {
                jump <- f*rep(1, length(id))
              }
            }

            if(is.null(jump)){
              jump <- rep(1, length(id))
              process <- unlist(process, use.names = FALSE)
            } else {
              compound <- tapply(jump, id, function(s) cumsum(c(0,s)))
              process <- unlist(lapply(levels(id), function(i) compound[[i]][process[[i]]+1]))
            }

            if(result == 'numeric') 
              return(process)

            ## TODO: could this be done more smooth, without dangerous
            ## manipulation of the internals of the object?
             
            CP <- as(object[ , colNames(object) %in% colnames(getUnitData(object))], "ContinuousProcess")
            valueEnv <- new.env(parent = .GlobalEnv)
            valueEnv$id <- getId(CP)
            valueEnv$position <- getPosition(CP)
            valueEnv$value <- Matrix(process, dimnames = list(NULL, "integrated"))
            valueEnv$i <- seq_along(valueEnv$id)
            valueEnv$j <- seq_len(dim(valueEnv$value)[2])
            CP@valueEnv <- valueEnv
            CP@colNames <- c(colNames(CP), "integrated")
            CP@jSubset <- 1L
            validObject(CP)
            MP <- data.frame(id, getPointPosition(object), "integrated", jump)
            names(MP) <- c(object@idVar, object@positionVar, object@markVar, "jump")
                                        
            return(jumpProcess(MP, CP, idVar = object@idVar, positionVar = object@positionVar, markVar = object@markVar))
          }
          )          
                 
setMethod("iPointSubset", "MarkedPointProcess",
          function(object) {
            i <- seq_along(object@pointProcessEnv$id)
            if(!isTRUE(object@iPointSubset == -1L))
              i <- i[object@iPointSubset]
            return(i)
          }
          )

setMethod("jPointSubset", "MarkedPointProcess",
          function(object) {
            j <- seq_len(dim(object@pointProcessEnv$markValue)[2])
            if(!isTRUE(object@jPointSubset == -1L))
              j <- j[object@jPointSubset]
             
            return(j)
          }
          )

setReplaceMethod("iPointSubset", c(object = "MarkedPointProcess", value = "ANY"),
                 function(object, value) {
                   if(length(value) == length(object@pointProcessEnv$id)) {
                     object@iPointSubset <- -1L
                   } else {              
                     object@iPointSubset <- value
                   }
                   return(object)
                 }
                 )

setReplaceMethod("jPointSubset", c(object = "MarkedPointProcess", value = "ANY"),
                 function(object, value) {
                   if(length(value) == dim(object@pointProcessEnv$markValue)[2]) {
                     object@jPointSubset <- -1L
                   } else {              
                     object@jPointSubset <- value
                   }
                   return(object)
                 }
                 )
              

setMethod("getPointId", "MarkedPointProcess",
          function(object, ...){
            if(isTRUE(object@iPointSubset == -1L)) {
              value <- object@pointProcessEnv$id
            } else {
              value <- object@pointProcessEnv$id[iPointSubset(object), drop = TRUE]
            }
            return(value)
          }
          )


setMethod("getPointPosition", "MarkedPointProcess",
          function(object, ...){
            if(isTRUE(object@iPointSubset == -1L)) {
              value <- object@pointProcessEnv$position
            } else {
              value <- object@pointProcessEnv$position[iPointSubset(object)]
            }
            return(value)
          }
          )

setMethod("getPointTime", "MarkedPointProcess",
          function(object, ...) {
            getPointPosition(object, ...)
          }
          )


setMethod("getMarkType", "MarkedPointProcess",
          function(object, ...){
            if(isTRUE(object@iPointSubset == -1L)) {
              value <- object@pointProcessEnv$markType
            } else {
              value <- object@pointProcessEnv$markType[iPointSubset(object), drop = TRUE]
            }
            
            return(value)
          }
          )

setMethod("getMarkValue", "MarkedPointProcess",
          function(object, ...){
            if(isTRUE(object@iPointSubset == -1L && object@jPointSubset == -1L)) {
                value <- object@pointProcessEnv$markValue
            } else if(isTRUE(object@iPointSubset == -1L)) {
              value <- object@pointProcessEnv$markValue[ ,jPointSubset(object), drop = FALSE]
            } else if(isTRUE(object@jPointSubset == -1L)) {
              value <- object@pointProcessEnv$markValue[iPointSubset(object), , drop = FALSE]
            } else {
              value <- object@pointProcessEnv$markValue[iPointSubset(object), jPointSubset(object), drop = FALSE]
            }
              
              return(value)
            }
          )

setMethod("getPlotData", "MarkedPointProcess",
          function(object, y = '@mark', nPoints = 200, allUnitData = FALSE, ...){
            if(length(getMarkType(object)) == 0) {

              plotData <- callGeneric(object = as(object, "ContinuousProcess"),
                                      nPoints = nPoints,
                                      allUnitData = allUnitData, ...)
              
            } else {
              
              pointPlotData <- data.frame(id = getPointId(object),
                                          position = getPointPosition(object),
                                          variable = factor(getMarkType(object)))
              
              plotData <- callGeneric(object = as(object, "ContinuousProcess"),
                                      nPoints = nPoints,
                                      allUnitData = allUnitData,
                                      selectPoints = getPointPointer(object), ...)
              
              if(isTRUE(y %in% names(getMarkValue(object))))
                pointPlotData <- cbind(pointPlotData, getMarkValue(object))

              if(isTRUE(allUnitData))
                pointPlotData <- cbind(pointPlotData, getUnitData(object)[as.numeric(getPointId(object)), , drop = FALSE])

              pointPlotData$type <- as.factor('Track')

              if(is.numeric(y)) {
                pointPlotData$value <- y
                plotData@breaks <- c(breaks, y)
                plotData@labels <- c(labels, as.character(y))
              } else {
                pointPlotData$value <- as.factor(y)
              }
              
              if(y == "@mark") 
                pointPlotData$value <- pointPlotData$variable
              if(y == object@idVar) 
                pointPlotData$value <- pointPlotData$id
              if(isTRUE(y %in% names(getMarkValue(object)))) 
                pointPlotData$value <-  getMarkValue(object)[ ,y]
              if(y == "@top")
                pointPlotData$value <- pointPlotData$type
              if(y == "@bottom") {
                pointPlotData$value <- pointPlotData$type
                plotData@position <- "bottom"
              }
              if("value" %in% names(plotData@continuousPlotData)) {
                if(isTRUE(y %in% levels(plotData@continuousPlotData$variable))) 
                  pointPlotData$value <- getValue(object)[getPointPointer(object), y]
              }
              
                  
              plotData@pointPlotData <- pointPlotData
             
            }
            
            return(plotData)
          }  
          )

setMethod("plot", c("MarkedPointProcess", "missing"),
          function(x, y, ...) callGeneric(x = x, y = '@mark', ...)
          )

setMethod("plot", c("MarkedPointProcess", "character"),
          function(x, y, nPoints = 200, ...){
            plotData <- getPlotData(object = x, y = y, nPoints = nPoints, ...)
            return(plot(plotData))
          }
          )

setMethod("subset", "MarkedPointProcess",
          function(x, ... , markSubset) {
            y <- callNextMethod()
            if (missing(markSubset)) 
              r <- TRUE  
            else {
              e <- substitute(markSubset)
              frame <- cbind(data.frame(getMarkType(y)),
                             getMarkValue(y))
              names(frame)[1] <- y@markVar
              r <- eval(e, frame, parent.frame())
              if (!is.logical(r)) 
                stop("'markSubset' must evaluate to logical")
              r <- r & !is.na(r) 
            }
            if(all(r)) {
              return(y)
            } else {
              return(y[-getPointPointer(y)[!r], ])
            }
          }
          )

setMethod("[", c(x = "MarkedPointProcess", i = "integer", j = "missing"),
          function(x, i, j, ... , drop = FALSE) {
            as(x, "ContinuousProcess") <- callGeneric(as(x, "ContinuousProcess"), i, , )
            i <- getPointPointer(x) %in% iSubset(x)
            iPointSubset(x) <- iPointSubset(x)[i]
            x@pointPointer <- match(getPointPointer(x)[i], iSubset(x))

            return(x)
          }
          )

setMethod("[", c(x = "MarkedPointProcess", i = "missing", j = "integer"),
          function(x, i, j, ... , drop = FALSE) {
            as(x, "ContinuousProcess")  <- callGeneric(as(x, "ContinuousProcess"), ,j, )
            j <- colnames(getMarkValue(x)) %in% colNames(x)
            i <- getMarkType(x) %in% colNames(x)
            if(!any(i) && drop) {
              x <- as(x,"ContinuousProcess")
            } else {
              jPointSubset(x) <- jPointSubset(x)[j]
              iPointSubset(x) <- iPointSubset(x)[i]
            }
            return(x)
          }
          )


setMethod("getPointPointer", c("MarkedPointProcess", "missing"),
          function(object, mark, ...) {
              value <- object@pointPointer
            return(value)
          }
          )

setMethod("getPointPointer", c("MarkedPointProcess", "character"),
          function(object, mark, ...) {
              value <- object@pointPointer[getMarkType(object) %in% mark]
            return(value)
          }
          )

setMethod("showData", "MarkedPointProcess",
          function(object, ....) {
            summaryById <- callNextMethod()
            structure <- summaryById$structure
            summaryById <- summaryById$summary
            if(length(getMarkType(object)) > 0) {
              id <- factor(getPointId(object), levels = rownames(summaryById))
              splitEntries <- split(seq_along(id), id)
              names(splitEntries) <- NULL
              pointSummary <- do.call("rbind", lapply(splitEntries, function(e) summary(getMarkType(object)[e])))
              colnames(pointSummary) <- paste("#", colnames(pointSummary), sep = "")
              summaryById <- cbind(summaryById,
                                   pointSummary)
              
              if(dim(getMarkValue(object))[2] > 0) {
                firstEntries <- sapply(splitEntries , function(e) e[1])
                summaryById <- cbind(summaryById,
                                     as.matrix(getMarkValue(object))[firstEntries, , drop = FALSE])
              }
            }
            
            markLevels <- levels(getMarkType(object))
            if(length(markLevels) > 4)
              markLevels <- c(markLevels[1:3], "...", markLevels[length(markLevels)])
            if(length(markLevels) > 1) {
              structure <- paste(structure, length(levels(getMarkType(object))),
                               " marks.\n   ", object@markVar,": ",
                               paste(markLevels, collapse = " "), "\n\n", sep="")
            }
            
            return(list(summary = summaryById, structure = structure))            
          }
          )

setMethod("updateProcessObject", "ContinuousProcess",
          function(object, ...) {
            pointProcess <- new("MarkedPointProcess",
                                unitData = object@unitData,
                                metaData = object@metaData,
                                iSubset = object@iSubset,
                                jSubset = object@jSubset,
                                idVar = object@idVar,
                                colNames = object@colNames,
                                positionVar = object@positionVar,
                                valueEnv = object@valueEnv,
                                factors = numeric(),
                                iPointSubset = object@iPointSubset,
                                jPointSubset = object@jPointSubset,
                                markVar =  object@markVar,
                                pointPointer = object@pointPointer,
                                pointProcessEnv = object@pointProcessEnv 
                                ) 
         
            validObject(pointProcess)
            return(pointProcess)
            
            }
          )
