## setAs("MarkedPointProcess", "ContinuousProcess",
##       def = function(from) {
##         jSubsetFrom <- jSubset(from)
##         to <- new("ContinuousProcess",
##                   equiDistance <- from@equiDistance,
##                   iSubset <- from@iSubset,
##                   type <- from@type,
##                   positionVar <- from@positionVar,
##                   metaData <- from@metaData,
##                   dimensions <- from@dimensions,
##                   colNames <- from@colNames,
##                   valueEnv <- from@valueEnv,
##                   idVar <- from@idVar,
##                   iUnitSubset <- from@iUnitSubset,
##                   jSubset <- jSubsetFrom[from@type[jSubsetFrom] %in% c("unit", "factor", "numeric")]
##                   )
##         return(to)
##       },
##       replace = function(from, value) {
##         if(from@valueEnv == value@valueEnv) {
##           from@iSubset <- value@iSubset
##           colNames(from) %in% colNames(value)

##         } else {
##           stop(paste("Cannot replace ContinuousProcess part of", object, "with", value))
##         }
##         return(from)
##       }
##       )


setMethod("markedPointProcess", c("data.frame", "ContinuousProcess"),
          function(pointData, continuousData, markVar = 'markType', coarsen = NULL, ...) {            

            ## Initialize an empty object - bypasses validity
            ## check. Validity function called before returning object
            ## at the end of the contructor.
            pointProcess <- new("MarkedPointProcess") 
##            as(pointProcess, "ContinuousProcess") <- continuousData
            
            if(continuousData@idVar %in% names(pointData)) {
              id <- pointData[ , continuousData@idVar]
              if(!identical(levels(id), levels(getId(continuousData))))
                id <- factor(id, levels = levels(getId(continuousData)))
            } else {
              id <- factor(rep(getId(continuousData)[1], dim(pointData)[1]))
            }

            ## Sorting the rows if needed.
            if(is.unsorted(id)) {
              ord <- order(id)              
              pointData <-  pointData[ord, , drop = FALSE]
              id <- id[ord]
            }

            if(continuousData@positionVar %in% names(pointData)) {
              position <- pointData[ , continuousData@positionVar]
            } else {
              stop(paste("pointData needs a column named", pointProcess@positionVar))
            }
            
            ip <- split(seq_along(id), id)
            if(any(sapply(names(ip), function(i) is.unsorted(position[ip[[i]]])))) {
              ord <- sapply(names(ip),
                            function(i) order(position[ip[[i]]]), simplify = FALSE)
              ord <- unlist(lapply(names(ip),
                                   function(i) ip[[i]][ord[[i]]]), use.names=FALSE)
              pointData <-  pointData[ord, , drop = FALSE]
              position <- position[ord]
              id <- id[ord]
              ip <- split(seq_along(id), id)
            }
                        
            if(markVar %in% names(pointData)) {
              markType <- pointData[ , markVar]
            } else {
              markType <-  rep("point", dim(pointData)[1])
            }
            markValue <- pointData[ ,!(names(pointData) %in%
                                       c(continuousData@idVar,
                                         continuousData@positionVar,
                                         markVar)),
                                   drop = FALSE]

            pointProcessEnv <- new.env(parent = .GlobalEnv)
            pointProcessEnv$id <- id
            pointProcessEnv$markType <- as.factor(markType)
            for(v in names(markValue)) {
              assign(paste("var_", v, sep = ""),
                     markValue[,v], envir = pointProcessEnv)
            }
            
           
            pointProcess@iPointSubset <- -1L
            pointProcess@jPointSubset <- -1L
            

            
            ## Compute a combined set of evaluation points and keep track of which
            ## points are points from the marked point process.
            ## pointId <- getPointId(pointProcess)
            ## contId <- getId(pointProcess)
            ## pointPosition <- getPointPosition(pointProcess)
            ## contPosition <- getPosition(pointProcess)
            ## id = factor(c(as.character(pointId), as.character(contId)), levels = levels(contId)),
            ##                            position = c(pointPosition, contPosition),
            ##                            pointer = c(seq_along(pointPosition), rep(Inf, length(contPosition))))

            ## evalPosition <- split(evalPosition, evalPosition$id)
            ## ## Dropping duplicated entries within 'id'.
            ## dup <- lapply(evalPosition, function(e) duplicated(e$position))
            ## evalPosition <- sapply(names(evalPosition), function(i) evalPosition[[i]][!(dup[[i]] & evalPosition[[i]]$pointer == Inf), ], simplify = FALSE)

            ## ## Order within 'id'.
            ## ord <- lapply(evalPosition, function(e) order(e$position))
            ## evalPosition <- sapply(names(evalPosition), function(i) evalPosition[[i]][ord[[i]], ],  simplify = FALSE)

            ## contPosition <- split(data.frame(pointer = seq_along(contPosition), position = contPosition), contId)

            ## ii <- sapply(names(evalPosition),
            ##              function(i) {
            ##                findInterval(evalPosition[[i]]$position,
            ##                             c(contPosition[[i]]$position, Inf),
            ##                             all.inside=TRUE)
            ##              },
            ##              simplify = FALSE)

            ## ii <- unlist(lapply(names(evalPosition), function(i) contPosition[[i]]$pointer[ii[[i]]]), use.names = FALSE)
            ic <- split(iSubset(continuousData), getId(continuousData))         
            contPosition <- getPosition(continuousData)
            ii <- sapply(names(ic),
                         function(i) { 
                           findInterval(position[ip[[i]]],
                                        contPosition[ic[[i]]])
                         }, simplify = FALSE                         
                         )
            ## Set any 0-pointer to point to the left most end point.
            ## Computes the pointers in terms of the full data set.
            l <- cumsum(c(0L, sapply(ic, length)[-length(ic)]))
            leftPoints <- which(unlist(ii) == 0L)
            ii <- as.integer(unlist(lapply(seq_along(l), function(i) ii[[i]] + l[[i]])), use.names = FALSE)
            ii[leftPoints] <- ii[leftPoints] + 1L

            if(is.null(coarsen)) {
              ## Finding exactly duplicated position entries.
              iii <- !(position %in% contPosition[ii])
              ## Setting pointer information.
              pointPointer <-  ii + cumsum(iii) ## which(orderii > length(i))
              pointPointer[leftPoints] <- pointPointer[leftPoints] - 1L
              ## Creating copy of the continuous data part if necessary.
              if(any(iii)) {
                
                i <- c(unlist(ic, use.names = FALSE), ii[iii])
                iSubset(continuousData) <- sort(i)
                as(pointProcess, "ContinuousProcess") <- continuousProcess(continuousData)
                ## evalPosition <- do.call("rbind", evalPosition)                 
                ## TODO: Move pointer to environment, remove points from pointprocess.
                pointProcess@valueEnv$position[pointPointer] <- position
              } else {
                as(pointProcess, "ContinuousProcess") <- continuousData
              }
            } else if(coarsen == 'right') {
              ## Setting pointer information
              pointPointer <-  pmin(ii + 1L, dim(continuousData)[1])
              pointPointer[leftPoints] <- pointPointer[leftPoints] - 1L
              as(pointProcess, "ContinuousProcess") <- continuousData
            } else if(coarsen == 'left') {
              ## Setting pointer information
              pointPointer <-  ii
              as(pointProcess, "ContinuousProcess") <- continuousData
            } else {
              stop("Argument 'coarsen' can only be 'NULL', 'left' or 'right'")
            }

            pointProcessEnv$pointPointer <- pointPointer
            pointProcess@pointPointer <- -1L
            pointProcess@pointProcessEnv <- pointProcessEnv
            lockEnvironment(pointProcess@pointProcessEnv, bindings = TRUE)
            
            pointProcess@markColNames <-levels(getMarkType(pointProcess))
            pointProcess@markValueColNames <- names(markValue)
         
            ## newContData <- c(list(evalPosition$id,
            ##                       evalPosition$position,
            ##                       getValue(continuousData)[ii, , drop=FALSE]),
            ##                  lapply(getFactors(continuousData), function(fac) fac[ii]))
            ## names(newContData)[1:3] <- c(continuousData@idVar,
            ##                              continuousData@positionVar,
            ##                              "value")
          
            ## as(pointProcess, "ContinuousProcess") <- continuousProcess(newContData)

            
            lockEnvironment(pointProcess@valueEnv, bindings = TRUE)
            validObject(pointProcess)
            return(pointProcess)
            
          }
          )

setMethod("markedPointProcess", c("MarkedPointProcess", "missing"),
          function(pointData, continuousData, ...) {
            continuousData <- continuousProcess(pointData)
            tmp <- data.frame(getPointId(pointData),
                              getPointPosition(pointData),
                              getMarkType(pointData))
            names(tmp) <- c(pointData@idVar,
                            pointData@positionVar,
                            "markType")
            pointData <- cbind(tmp, getMarkValue(pointData))
            if(dim(pointData)[1] > 0) {
              return(callGeneric(pointData, continuousData, ...))
            } else {
              warning("No point process data in marked point process, returning a 'ContinuousProcess'", call = FALSE)
              return(continuousData)
            }
               
          }
          )


setMethod("markedPointProcess", c("data.frame", "data.frame"),
          function(pointData, continuousData, ..., markVar = 'markType', coarsen = NULL) {
            callGeneric(pointData = pointData,
                        continuousData = continuousProcess(continuousData, ...),
                        markVar = markVar, coarsen = coarsen)
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

setMethod("colNames", c("MarkedPointProcess", "missing"),
          function(object, ...) {
            colnames <- c(object@markColNames, object@markValueColNames)
            if(!identical(object@jPointSubset[1], -1L)) 
              colnames <- colnames[object@jPointSubset]

            colnames <- c(callNextMethod(object), colnames)
            return(colnames)
          }
          )


setMethod("colNames", c("MarkedPointProcess", "character"),
          function(object, type, ...) {
            colnames <- callGeneric(object = object, ...)
            if(type == "mark") {
              colnames <- colnames[colnames %in% object@markColNames]
            } else if(type == "markValue") {
              colnames <- colnames[colnames %in% object@markValueColNames]             
            } else {
              colnames <- callGeneric(as(object, "ContinuousProcess"), type = type, ...) 
            }
            return(colnames)
          }
          )

setMethod("dim", "MarkedPointProcess",
          function(x) {
            d <- callNextMethod(x)
          
            if(identical(x@jPointSubset[1], -1L)) {
              d2 <- length(x@markColNames) + length(x@markValueColNames) 
            } else {
              d2 <- length(x@jPointSubset)
            }
            d[2] <- d[2] + d2

            return(d)
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
             
            ## CP <- as(object[ , colNames(object) %in% colnames(getUnitData(object))], "ContinuousProcess")
            ## valueEnv <- new.env(parent = .GlobalEnv)
            ## valueEnv$id <- getId(CP)
            ## valueEnv$position <- getPosition(CP)
            ## valueEnv$value <- Matrix(process, dimnames = list(NULL, "integrated"))
            ## valueEnv$i <- seq_along(valueEnv$id)
            ## valueEnv$j <- seq_len(dim(valueEnv$value)[2])
            ## CP@valueEnv <- valueEnv
            ## CP@colNames <- c(colNames(CP), "integrated")
            ## CP@jSubset <- 1L
            ## validObject(CP)
            CP <- continuousProcess(c(as(object[ , FALSE], "list"),
                                      list(integrated = process)),
                                    idVar = object@idVar,
                                    positionVar = object@positionVar,
                                    unitData = getUnitData(object))
            MP <- data.frame(id, getPointPosition(object), "integrated", jump)
            names(MP) <- c(object@idVar, object@positionVar, "markType", "jump")
                                        
            return(jumpProcess(MP, CP, idVar = object@idVar,
                                      positionVar = object@positionVar))
          }
          )          
                 
setMethod("iPointSubset", "MarkedPointProcess",
          function(object) {
            if(identical(object@iPointSubset[1], -1L)) {
              i <- seq_along(object@pointProcessEnv$id)
            } else {
              i <- object@iPointSubset
            }
            return(i)
          }
          )

setMethod("jPointSubset", "MarkedPointProcess",
          function(object) {
            j <- seq_len(length(object@markColNames) + length(object@markValueColNames))
            if(!isTRUE(object@jPointSubset == -1L))
              j <- j[object@jPointSubset]
             
            return(j)
          }
          )

setReplaceMethod("iPointSubset", c(object = "MarkedPointProcess", value = "ANY"),
                 function(object, value) {
                   if(length(value) == length(object@pointProcessEnv$id) &&
                      identical(value, seq_along(object@pointProcessEnv$id))) {
                     object@iPointSubset <- -1L
                   } else {
                     object@iPointSubset <- value
                     ## We find the, potentially, new pointers. 
                     object@pointPointer <- match(object@pointProcessEnv$pointPointer[value],
                                                  iSubset(object))
                   }
                   return(object)
                 }
                 )

setReplaceMethod("jPointSubset", c(object = "MarkedPointProcess", value = "ANY"),
                 function(object, value) {
                   d2 <- length(object@markColNames) + length(object@markValueColNames)
                   if(length(value) == d2 && identical(value, seq_len(d2))) {
                     object@jPointSubset <- -1L
                   } else {              
                     object@jPointSubset <- value
                   }
                   return(object)
                 }
                 )
              

setMethod("getPointId", "MarkedPointProcess",
          function(object, drop = TRUE, ...){
            if(isTRUE(object@iPointSubset == -1L)) {
              value <- object@pointProcessEnv$id
            } else {
              value <- object@pointProcessEnv$id[iPointSubset(object), drop = TRUE]
              
              if(!drop)
                value <- factor(value, levels = levels(getId(object)))
            }
            return(value)
          }
          )


setMethod("getPointPosition", "MarkedPointProcess",
          function(object, ...){
            return(getPosition(object)[getPointPointer(object, ...)])
          }
          )

setMethod("getPointTime", "MarkedPointProcess",
          function(object, ...) {
            getPointPosition(object, ...)
          }
          )


setMethod("getMarkType", "MarkedPointProcess",
          function(object, drop = TRUE, ...){
            if(isTRUE(object@iPointSubset == -1L)) {
              value <- object@pointProcessEnv$markType
            } else {
              value <- object@pointProcessEnv$markType[iPointSubset(object), drop = TRUE]
              if(!drop)
                value <- factor(value, levels = colNames(object, "mark"))
            }
            
            return(value)
          }
          )

setMethod("getMarkValue", "MarkedPointProcess",
          function(object, ...){
            j <- colNames(object, "markValue")
            if(length(j) > 0) {
              value <- as.data.frame(getColumns(object, j, drop = FALSE))
            } else {
              value <- data.frame()[seq_along(getPointId(object)),]
            }
            ## if(isTRUE(object@iPointSubset == -1L && object@jSubset == -1L)) {
            ##     value <- object@pointProcessEnv$markValue
            ## } else if(isTRUE(object@iPointSubset == -1L)) {
            ##   value <- object@pointProcessEnv$markValue[ ,jSubset(object), drop = FALSE]
            ## } else if(isTRUE(object@jSubset == -1L)) {
            ##   value <- object@pointProcessEnv$markValue[iPointSubset(object), , drop = FALSE]
            ## } else {
            ##   value <- object@pointProcessEnv$markValue[iPointSubset(object), jSubset(object), drop = FALSE]
            ## }
              
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
                pointPlotData$value <-  getColumns(object, y)
              if(y == "@top") {
                pointPlotData$value <- pointPlotData$type
                plotData@position <- "top"
              }
              if(y == "@bottom") {
                pointPlotData$value <- pointPlotData$type
                plotData@position <- "bottom"
              }
              if("value" %in% names(plotData@continuousPlotData)) {
                if(y == "@variable") {
                  pointPlotData$value <- numeric(dim(pointPlotData)[1])
                  for(variable in levels(pointPlotData$variable)) {
                    pointPlotData$value[pointPlotData$variable == variable] <-
                      getColumns(object[getPointPointer(object, variable),
                                        colNames(as(object, "ContinuousProcess"))],
                                        variable)
                  }
                } else if(isTRUE(y %in% levels(plotData@continuousPlotData$variable))) {
                    pointPlotData$value <- getColumns(object, y)[getPointPointer(object)]
                  }
                
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
            return(plot(plotData, ...))
          }
          )

setMethod("subset", "MarkedPointProcess",
          function(x, ... , pointSubset) {
            y <- callNextMethod()
            if (missing(pointSubset)) 
              r <- TRUE  
            else {
              e <- substitute(pointSubset)
              variables <- all.vars(e)
              tmpEnv <- new.env(parent = .GlobalEnv)

              for(v in variables) {
                assign(v, getColumns(x, v), envir = tmpEnv) 
              }
             
              r <- eval(e, tmpEnv, parent.frame())
              if (!is.logical(r)) 
                stop("'pointSubset' must evaluate to logical")
              r <- r & !is.na(r) 
            }
            if(!all(r)) {
              iPointSubset(y) <- iPointSubset(y)[r]
            } 
            return(y)
          }
          )


setMethod("getColumns", c("MarkedPointProcess", "character"),
          function(object, j, drop = TRUE) {
            checkColumns <- j %in% colNames(object)
            if(!all(checkColumns)) 
              stop(paste(c("No column '", j[!checkColumns][1], "' in the object."),
                         collapse = ""), call. = FALSE)
            
            contj <- j %in% colNames(as(object, "ContinuousProcess"))
            
            if(drop && length(j) == 1) {
              if(contj) {
                column <- callNextMethod(object, j, drop = TRUE)
              } else {
                if(j %in% object@markColNames) {
                  column <- getPointPosition(object, j)
                } else {
                  column <- get(paste("var_", j, sep = ""),
                                envir = object@pointProcessEnv)
                }
                if(!identical(object@iPointSubset[1], -1L))
                  column <- column[object@iPointSubset]
              }
              } else {
                column <- callNextMethod(object, j[contj], drop = drop)
                for(jj in j[!contj]) {
                  if(jj %in% object@markColNames) {
                    column[[jj]] <- getPointPosition(object, jj)
                  } else {
                    column[[jj]] <- get(paste("var_", jj, sep = ""),
                                        envir = object@pointProcessEnv)
                  }
                  if(!identical(object@iPointSubset[1], -1L))
                    column[[jj]] <- column[[jj]][object@iPointSubset]
                }
              }
                        
            return(column)
          }
          )

setMethod("[", c(x = "MarkedPointProcess", i = "integer", j = "missing"),
          function(x, i, j, ... , drop = FALSE) {
            as(x, "ContinuousProcess") <- callGeneric(as(x, "ContinuousProcess"), i, , )
##            i <- x@pointProcessEnv$pointPointer %in% iSubset(x)
##            x@pointPointer <- match(getPointPointer(x)[i], iSubset(x))
            iPointSubset(x) <- iPointSubset(x)[getPointPointer(x) %in% iSubset(x)]
            if(drop) {
              marks <- colNames(x, "mark")
              dropCol <- colNames(x) %in% marks[!(marks %in% levels(getMarkType(x)))]
              if(any(dropCol)) 
                x <- x[ , !dropCol, drop = TRUE]
            }
              
            return(x)
          }
          )

setMethod("[", c(x = "MarkedPointProcess", i = "missing", j = "character"),
          function(x, i, j, ... , drop = FALSE) {
            if(drop && length(j) == 1) {
              x <- getColumns(x, j)
            } else {
              as(x, "ContinuousProcess") <- callGeneric(as(x, "ContinuousProcess"), , j)             
              i <- getMarkType(x) %in% j
              if(!any(i)) {
                if(drop) {
                  x <- as(x, "ContinuousProcess")
                } else {
                  iPointSubset(x) <- integer()
                  jPointSubset(x) <- integer()
                }
              } else {
                iPointSubset(x) <- iPointSubset(x)[i]
                jPointSubset(x) <- which(c(x@markColNames, x@markValueColNames) %in% j)
              }
            }
            
            return(x)
          }
          )

setMethod("getPointPointer", c("MarkedPointProcess", "missing"),
          function(object, mark, ...) {
            if(identical(object@pointPointer[1], -1L)) {
              value <- object@pointProcessEnv$pointPointer
            } else {
              value <- object@pointPointer
            }
            return(value)
          }
          )

setMethod("getPointPointer", c("MarkedPointProcess", "character"),
          function(object, mark, ...) {
            getPointPointer(object)[getMarkType(object) %in% mark]
          }
          )

setMethod("object.size", "ProcessData", 
          function(x) {
            size <- sum(sapply(ls(x@vpointProcessEnv),
                               function(name) object.size(get(name, envir = x@vpointProcessEnv))))
            size <- size + object.size(as(x, "ContinuousProcess"))
            return(structure(size, class = "object_size"))
          }
          )

setMethod("summarizeData", "MarkedPointProcess",
          function(object, ....) {
            summaryById <- callNextMethod()
            if(length(colNames(object, "mark")) > 0) {
              id <- getPointId(object, drop = FALSE)
              splitEntries <- split(seq_along(id), id)
              names(splitEntries) <- NULL
              pointSummary <- do.call("rbind", lapply(splitEntries, function(e) summary(getMarkType(object, drop = FALSE)[e])))
              colnames(pointSummary) <- paste("#", colnames(pointSummary), sep = "")
              summaryById <- cbind(summaryById, pointSummary)
              sumVal <- list()
              for(j in colNames(object, "markValue")) {
                column <- getColumns(object, j)
                sumVal[[paste("mean(", j ,")", sep ="")]] <- 
                  sapply(splitEntries, function(e) mean(column[e]))
              }
              if(length(sumVal) > 0) 
                summaryById <- cbind(summaryById,
                                     as.data.frame(sumVal, optional = TRUE))
              
            }            
            return(summaryById)            
          }
          )

## setMethod("updateProcessObject", "ContinuousProcess",
##           function(object, ...) {
##             pointProcess <- new("MarkedPointProcess",
##                                 unitData = object@unitData,
##                                 metaData = object@metaData,
##                                 iSubset = object@iSubset,
##                                 jSubset = object@jSubset,
##                                 idVar = object@idVar,
##                                 colNames = object@colNames,
##                                 positionVar = object@positionVar,
##                                 valueEnv = object@valueEnv,
##                                 factors = numeric(),
##                                 iPointSubset = object@iPointSubset,
## #                                jPointSubset = object@jPointSubset,
##                                 pointPointer = object@pointPointer,
##                                 pointProcessEnv = object@pointProcessEnv 
##                                 ) 
         
##             validObject(pointProcess)
##             return(pointProcess)
            
##             }
##           )

setMethod("unsubset", "MarkedPointProcess",
          function(x, ...) {
            as(x, "ContinuousProcess") <- callGeneric(as(x, "ContinuousProcess"))
            
            x@iPointSubset <- -1L
            x@jPointSubset <- -1L
            x@pointPointer <- -1L
            return(x)
          }
          )
