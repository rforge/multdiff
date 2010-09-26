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
            
            ## ## Order within 'id'.
            
            ## ord <- tapply(evalPosition$position,
            ##               evalPosition$id,
            ##               order)
            
            ## ord <- unlist(lapply(levels(evalPosition$id),
            ##                      function(id) which(evalPosition$id == id)[ord[[id]]]),
            ##               use.names=FALSE)
            
            ## ## Dropping duplicated entries.
            ## dup <- unlist(tapply(evalPosition$position,
            ##                      evalPosition$id,
            ##                      duplicated),
            ##               use.names = FALSE)        
            
            ## evalPosition <- evalPosition[unlist(tapply(seq_along(evalPosition$id),
            ##                                            evalPosition$id,
            ##                                            function(s) s))[!dup], ]

            pointProcess@pointPointer <- which(evalPosition$pointer < Inf)

            
            ## ii <- sapply(levels(evalPosition$id),
            ##                    function(i) {
            ##                      findInterval(evalPosition[ord, 'position'][evalPosition[ord,"id"] == i],
            ##                                   c(contPosition[contId == i], Inf),
            ##                                   all.inside=TRUE)},
            ##             simplify=FALSE)

            ## ii <- as.numeric(unlist(lapply(levels(evalPosition$id),
            ##                                function(i) {
            ##                                  which(contId == i)[ii[[i]]]
            ##                                }),
            ##                         use.names=FALSE))
          
            valueEnv <- new.env(parent = .GlobalEnv)
            valueEnv$id <- evalPosition$id
            valueEnv$position <- evalPosition$position
            valueEnv$value <- getValue(pointProcess)[ii, , drop=FALSE]
            rownames(valueEnv$value) <- NULL
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
            callGeneric(pointData, continuousProcess(continuousData, ...), ...)
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
            plotData <- callGeneric(as(object, "ContinuousProcess"), nPoints = nPoints, allUnitData = allUnitData, ...)
            plotPointData <- data.frame(id = getPointId(object),
                                        position = getPointPosition(object),
                                        variable = factor(getMarkType(object)))
            if(isTRUE(y %in% names(getMarkValue(object))))
              plotPointData <- cbind(plotPointData, getMarkValue(object))

            if(isTRUE(allUnitData))
              plotPointData <- cbind(plotPointData, getUnitData(object)[as.numeric(getPointId(object)), , drop = FALSE])

            plotData$contVariable <- plotData$variable
            plotPointData$contVariable <- as.factor('points')
            
            if("value" %in% names(plotData)) { 
              low <- signif(min(plotData$value),1)
              high <- signif(max(plotData$value),1)
              pLevels <- levels(plotData$variable)
            } else {
              low <- 0
              high <- 1
            }

            if(y == "@mark") {
                variableNum <- as.numeric(plotPointData$variable)
                valueMark <- variableNum*(high-low)/(max(variableNum)+1)
                plotPointData$value <- valueMark + high
              } else if(y == object@idVar) {
                variableNum <- as.numeric(plotPointData$id)
                valueId <- variableNum*(high-low)/(max(variableNum)+1)
                plotPointData$value <- valueId + high
              } else if(isTRUE(y %in% names(getMarkValue(object)))) {
                plotPointData$value <-  getMarkValue(object)[ ,y]
              } else if("value" %in% names(plotData)) {
                if(y == "@top") {
                  plotPointData$value <- high + (high-low)/4
                  plotPointData$contVariable <- factor(pLevels[1], levels = pLevels)
                } else if(y == "@bottom") {
                  plotPointData$value <- low - (high-low)/4
                  plotPointData$contVariable <- factor(pLevels[length(pLevels)], levels = pLevels)
                } else if(isTRUE(y %in% colnames(getValue(object)))) {
                  plotPointData$value <- getValue(object)[getPointPointer(object), y]
                  plotPointData$contVariable <- as.factor(y)
                }
              } else {
                if(!is.numeric(y)) 
                  plotPointData$value <- y
              }
            
            return(list(plotData = plotData, plotPointData = plotPointData))
          }
          )

setMethod("plot", c("MarkedPointProcess", "missing"),
          function(x, y, ...) callGeneric(x = x, y = '@mark', ...)
          )

setMethod("plot", c("MarkedPointProcess", "character"),
          function(x, y, nPoints = 200, ...){
            if(length(getMarkType(x)) == 0) {
              ## There is only continuous process data in the data set.
              p <- callGeneric(as(x, "ContinuousProcess"), nPoints = nPoints, ...)
            } else {
              plotData <- getPlotData(x, y, values = values, nPoints = nPoints, ...)
              plotPointData <- plotData$plotPointData
              plotData <- plotData$plotData
              facetFormula <- as.formula(paste(x@idVar, "~ ."))


              ## Setting up breaks, labels and limits
              
              if("value" %in% names(plotData)) { 
                low <- signif(min(plotData$value),1)
                high <- signif(max(plotData$value),1)
                breaks <- pretty(c(low,high), n = 4)
              } else {
                low <- 0
                high <- 1
                breaks <- numeric()
                limits <- NULL
              }
              
              if(y == "@mark") {
                varLevels <- levels(plotPointData$variable)
                labels <- c(as.character(breaks), varLevels)
                breaks <- c(breaks, high + seq_along(varLevels)*(high-low)/(length(varLevels)+1))
                limits <- c(1,2)
              } else if(y == x@idVar) {
                idLevels <- levels(plotPointData$id)
                labels <- c(as.character(breaks), idLevels)
                breaks <- c(breaks, high + seq_along(idLevels)*(high-low)/(length(idLevels)+1))
                limits <- c(1,2)
              } else if(y == "@top") {
                labels <- c(as.character(breaks), "points")
                breaks <- c(breaks, high + (high-low)/4)
              } else if(y == "@bottom") {
                labels <- c(as.character(breaks), "points")
                breaks <- c(breaks, low - (high-low)/4)
              } else if(isTRUE(y %in% names(getMarkValue(x)))) {
                breaks <- pretty(c(low,high,plotPointData$value), n = 5)
                labels <- as.character(breaks)
              } else if(isTRUE(y %in% colnames(getValue(x)))) {
                breaks <- pretty(c(low,high,plotPointData$value), n = 5)
                labels <- as.character(breaks)
              } else {
                breaks <- c(breaks, y)
                labels <- as.character(breaks)
              }              
              
              if("value" %in% names(plotData)) {
                ## There is continuous process as well as point process data in the data set. 
                group <- paste(x@idVar, ":variable", sep = "")
                p <-  ggplot(data = plotData,  aes_string(x = "position", y = "value", colour = "variable", group = group)) +
                  facet_grid(facetFormula, scales = "free_y") +
                    scale_x_continuous(name = x@positionVar) +
                        geom_point(data = plotPointData) +
                          geom_line() +
                            scale_y_continuous(breaks = breaks, name = "", labels = labels)
                
                
              } else {
                ## There is only point process data in the data set.
                p <- ggplot(data = plotPointData, aes(x = position, y = value, colour = variable)) +
                  scale_y_continuous(breaks = breaks, name = "", labels = labels, limits = limits) + 
                    scale_x_continuous(name = x@positionVar) + 
                      facet_grid(facetFormula, scales = "free_y") +
                        geom_point()              
              }

            }
              
            return(p)
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
              value <- object@pointPointer[getMarkType(object) == mark]
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
