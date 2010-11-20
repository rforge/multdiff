setMethod("continuousProcess", "data.frame",
          function(continuousData, unitData = data.frame(), metaData = list(), positionVar = 'time', idVar = 'id', ...){
            if(!is.data.frame(continuousData))
              stop("The 'initdata' needs to be a data frame")
            
            colNames <- names(continuousData)
            if(!(idVar %in% colNames)) {
              id <- factor(rep(1, dim(continuousData)[1]))
            } else {
              id <- factor(continuousData[ ,idVar])
            }
            if(!(positionVar %in% colNames)) {
              position <- as.numeric(sapply(levels(id), function(i) seq_along(id[id == i])))
            } else {
              position <- continuousData[ ,positionVar]
            }
            if(dim(unitData)[2] == 0)
              unitData <- data.frame(row.names = levels(id))
            if(idVar %in% names(unitData)) {
              rownames(unitData) <- unitData[ , names(unitData) == idVar]
              unitData <- unitData[ , names(unitData) != idVar, drop = FALSE]
            }

            unitData <- unitData[levels(id), , drop = FALSE] 

            if(any(rownames(unitData) != levels(id))) {
              ord <- match(levels(id), rownames(id))
              unitData <- unitData[ord, ]
            }
            
            if(any(tapply(position, id, is.unsorted))) {
              ord <- tapply(position, id, order)
              ord <- unlist(lapply(levels(id),
                                   function(i) which(id == i)[ord[[i]]]),
                            use.names=FALSE)
              value <- continuousData[ord, !(colNames %in% c(idVar,positionVar)), drop = FALSE]
            } else {
              value <- continuousData[ , !(colNames %in% c(idVar,positionVar)), drop = FALSE]
            }

            valueEnv <- new.env(parent = .GlobalEnv)
            valueEnv$id <- id
            valueEnv$position <- position
            if(all(sapply(value, class) %in% c("integer", "numeric"))){
              valueEnv$value <- Matrix(as.matrix(value))
              factors <- numeric()
            } else {
              valueEnv$value <- model.Matrix(~., data = value,
                                             sparse = TRUE,
                                             row.names = FALSE)
              assign <- attr(valueEnv$value, "assign")
              intercept <- which(assign == 0)
              valueEnv$value <- valueEnv$value[, -intercept]
              assign <- assign[-intercept]
              factors <- which(assign %in% which(!(sapply(value, class) %in% c("integer", "numeric"))))
            }
              
            valueEnv$i <- seq_along(id)
            valueEnv$j <- seq_len(dim(valueEnv$value)[2])

            colNames <- c(names(unitData), colnames(valueEnv$value))
            
            new("ContinuousProcess",
                unitData = unitData,
                metaData = metaData,
                iSubset = -1L,
                jSubset = -1L,
                idVar = idVar,
                colNames = colNames,
                positionVar = positionVar,
                valueEnv = valueEnv,
                factors = factors)
          }
          )

setMethod("dim", "ContinuousProcess",
          function(x) {
            d1 <- length(getId(x))
            d2 <- length(colNames(x))
            return(c(d1,d2))
          }
          )

setMethod("colNames", "ContinuousProcess",
          function(object, ...) {
            return(object@colNames)
          }
          )

setReplaceMethod("colNames", c(object = "ContinuousProcess", value = "character"),
                 function(object, value) {
                   object@colNames <- value
                   return(object)
                 }
                 )

setMethod("iSubset", "ContinuousProcess",
          function(object) {
            i <- object@valueEnv$i
            if(!isTRUE(object@iSubset == -1L))
              i <- i[object@iSubset]
            return(i)
          }
          )

setMethod("jSubset", "ContinuousProcess",
          function(object) {
            j <- object@valueEnv$j
            if(!isTRUE(object@jSubset == -1L))
              j <- j[object@jSubset]
             
            return(j)
          }
          )

setReplaceMethod("iSubset", c(object = "ContinuousProcess", value = "ANY"),
                 function(object, value) {
                   if(length(value) == length(object@valueEnv$id)) {
                     object@iSubset <- -1L
                   } else {              
                     object@iSubset <- value
                   }
                   
                   return(object)
                 }
                 )

setReplaceMethod("jSubset", c(object = "ContinuousProcess", value = "ANY"),
                 function(object, value) {
                   if(length(value) == dim(object@valueEnv$value)[2]) {
                     object@jSubset <- -1L
                   } else {              
                     object@jSubset <- value
                   }
                   
                   return(object)
                 }
                 )

setMethod("[", c(x = "ContinuousProcess", i = "integer", j = "missing"),
          function(x, i, j, ... , drop = FALSE) {
            iSubset(x) <- iSubset(x)[i]
            as(x, "ProcessData") <- callGeneric(as(x, "ProcessData"), levels(getId(x)), , drop = drop)
            return(x)
          }
          )

setMethod("[", c(x = "ContinuousProcess", i = "numeric", j = "missing"),
          function(x, i, j, ... , drop = FALSE) {
            i <- as.integer(i)
            x <- callGeneric(x, i, , drop = drop)
            return(x)
          }
          )

setMethod("[", c(x = "ContinuousProcess", i = "logical", j = "missing"),
          function(x, i, j, ... , drop = FALSE) {
            i <- iSubset(x)[i]
            x <- callGeneric(x, i, , drop = drop)
            return(x)
          }
          )

setMethod("[", c(x = "ContinuousProcess", i = "missing", j = "integer"),
          function(x, i, j, ... , drop = FALSE) {
            colNames <- colNames(x)[j]
            jj <- colnames(getValue(x)) %in% colNames
            jSubset(x) <- jSubset(x)[jj]            
            j <- colNames(as(x, "ProcessData")) %in% colNames
            as(x, "ProcessData")  <- callGeneric(as(x, "ProcessData"), ,j, drop = drop)
            colNames(x) <- colNames
            return(x)
          }
          )

setMethod("[", c(x = "ContinuousProcess", i = "missing", j = "numeric"),
          function(x, i, j, ... , drop = FALSE) {
            j <- as.integer(j)
            x <- callGeneric(x, , j, drop = drop)
            return(x)
          }
          )

setMethod("[", c(x = "ContinuousProcess", i = "missing", j = "logical"),
          function(x, i, j, ... , drop = FALSE) {
            j <- seq_along(colNames(x))[j]
            x <- callGeneric(x, , j, drop = drop)
            return(x)
          }
          )

setMethod("[", c(x = "ContinuousProcess", i = "missing", j = "character"),
          function(x, i, j, ... , drop = FALSE) {
            j <- colNames(x) %in% j
            x <- callGeneric(x, , j, drop = drop)
            return(x)
          }
          )

setMethod("[", "ContinuousProcess",
          function(x, i, j, ... , drop = FALSE) {

            if(!(class(i) %in% c("logical", "numeric", "integer")))
              stop("i must be a 'logical' or a 'numeric' vector.")
            if(!(class(j) %in% c("logical", "numeric", "integer", "character")))
              stop("j needs to be a 'logical', a 'numeric' or a 'character' vector.")
            
            x <- callGeneric(x, i, , drop = drop)
            x <- callGeneric(x, , j, drop = drop)
            return(x)
          }
          )

setMethod("[", c("ContinuousProcess", "missing", "missing", "missing"),
          function(x, i, j, ... , drop) {
            return(x)
          }
          )

setMethod("subset", "ContinuousProcess",
          function(x, subset, select, ...) {
            if (missing(subset)) 
              r <- TRUE 
            else {
              e <- substitute(subset)
              unitVar <- names(getUnitData(x))
              variables <- all.vars(e)
              unitVariables <- variables[variables %in% unitVar]
              variables <- variables[!(variables %in% unitVariables)]
              variables <- variables[!(variables %in% c(x@idVar, x@positionVar))]
              frame <- data.frame(getId(x), getPosition(x))
              frame <- cbind(frame,
                             as.data.frame(as.matrix(getValue(x)[ ,variables])),
                             getUnitData(x)[as.numeric(getId(x)), unitVariables])
              
              ## The construction assumes that the level of id and
              ## column id in unitData are in the same order. This is
              ## assumed and enforced by the validity check. 

              names(frame) <- c(x@idVar, x@positionVar, variables, unitVariables)
              
              r <- eval(e, frame, parent.frame())
              
              if (!is.logical(r)) 
                stop("'subset' must evaluate to logical")
              r <- r & !is.na(r)
            }
            if (missing(select)) 
              vars <- TRUE
            else {
              nl <- as.list(seq_along(colNames(x)))
              names(nl) <- colNames(x)
              vars <- eval(substitute(select), nl, parent.frame())
            }
            return(x[r, vars])
          }
          )

setMethod("getId", "ContinuousProcess",
          function(object, ...) {
            if(isTRUE(object@iSubset == -1L)) {
              value <- object@valueEnv$id
            } else {
              value <- object@valueEnv$id[iSubset(object), drop =TRUE]
            }
            
            return(value)
           }
          )

setMethod("getFactors", "ContinuousProcess",
          function(object, ...) {
            factors <- numeric()
            tryCatch(factors <- object@factors,
                     error = function(e) {
                       warning("Object is deprecated and lacks the factor slot. Either recreate the\n object from scratch or try the 'updateProcessObject' function.", call. = FALSE)}
                     )
            if(!isTRUE(object@jSubset == -1L)) {
              factors <- which(jSubset(object) %in% factors)
            }
            
            return(factors)
           }
          )

setMethod("getValue", "ContinuousProcess",
          function(object, ...) {
            if(isTRUE(object@iSubset == -1L && object@jSubset == -1L)) {
              value <- object@valueEnv$value
            } else if(isTRUE(object@iSubset == -1L)) {
              value <- object@valueEnv$value[ ,jSubset(object), drop = FALSE]
            } else if(isTRUE(object@jSubset == -1L)) {
              value <- object@valueEnv$value[iSubset(object), , drop = FALSE]
            } else {
              value <- object@valueEnv$value[iSubset(object), jSubset(object), drop = FALSE]
            }
            
            return(value)
          }
          )

setMethod("getPlotData", "ContinuousProcess",
          function(object, nPoints = 200, allUnitData = FALSE, selectPoints = NULL, ...){
            factorPlotData <- data.frame()
            factors <- getFactors(object)
            if(length(factors) > 0) {

              patterns <- apply(getValue(object)[, factors, drop = FALSE] != 0, 2, function(x) {
                tmp <- tapply(seq_along(x), getId(object), function(j) {
                  z <- x[j]
                  n <- length(z)
                  i <- which(z[-1] != z[-n])
                  names(i) <- NULL
                  i[!z[i]] <-  i[!z[i]] + 1
                  if(z[1]) 
                    i <- c(1,i)
                  if(z[n])
                    i <- c(i,n)
                  if(length(i) == 0)
                    i <- NA
                  return(getPosition(object)[j][i])
                })
                attributes(tmp) <- NULL
                names(tmp) <- levels(getId(object))
                return(tmp)
              })              
              
              uniqueNames <- unlist(lapply(patterns, function(y) {
                lapply(y, function(z) {
                  i <- seq(1,length(z),2)
                  return(rep(z[i], each = 2, length.out = length(z)))
                })
              }), use.names = FALSE)

              patterns <- melt(patterns)
              names(patterns) <- c("position", "id", "variable")
              patterns$id <- factor(patterns$id, levels = levels(getId(object)))
              patterns$variable <- factor(patterns$variable)              
              patterns$group <- paste(patterns$id,
                                      patterns$variable,
                                      uniqueNames, sep = "")

              factorPlotData <- patterns[!is.na(patterns$position), ]
              factorPlotData$value <- factorPlotData$variable
              
              if(isTRUE(allUnitData))            
                factorPlotData <- cbind(plotData, getUnitData(object)[as.numeric(plotData$id), , drop = FALSE])

              factorPlotData$type <- as.factor("Track")
              
            }
            
            iList <- unlist(tapply(seq_along(getId(object)),
                                     getId(object),
                                     list,
                                     simplify = FALSE),
                            recursive = FALSE)
            
            i <- unlist(lapply(iList,
                               function(ii) {
                                 if(is.null(ii))
                                   return(NULL)
                                 l <- length(ii)
                                   as.integer(seq(ii[1], ii[l], length.out = min(l,nPoints)))
                               }))
            
            if(!is.null(selectPoints))
              i <- unique(sort(c(i, selectPoints)))
            
            object <- object[i, ]
            tmp <- as.matrix(getValue(object))
            if(length(factors) > 0) 
              tmp <- tmp[, -factors, drop = FALSE]

            rownames(tmp) <- NULL
            measureVar <- colnames(tmp)
            tmp <- cbind(tmp, data.frame(iSubset = I(i)))
            
            continuousPlotData <- data.frame(getId(object))
            names(continuousPlotData)[1] <- object@idVar
            
            if(isTRUE(allUnitData))            
              continuousPlotData <- cbind(continuousPlotData, getUnitData(object)[as.numeric(getId(object)), , drop = FALSE])
            
            continuousPlotData <- melt(cbind(continuousPlotData,
                                             position = getPosition(object),
                                             as.data.frame(tmp)),
                                       measure.vars = measureVar)
            
            continuousPlotData$type <- as.factor("Continuous")

            if("value" %in% names(continuousPlotData)) {
              limits <- range(continuousPlotData$value)
              breaks <- pretty(limits, 4)
              labels <- as.character(breaks)
            } else {
              limits <- c(-1,0)
              breaks <- numeric()
              labels <- character()
            }
            
            plotData <- new("ProcessPlotData",
                            continuousPlotData = continuousPlotData,
                            factorPlotData = factorPlotData,
                            pointPlotData = data.frame(),
                            position = "top",
                            limits = limits,
                            breaks = breaks,
                            labels = labels,
                            idVar = object@idVar,
                            positionVar = object@positionVar)
            
            return(plotData)
          }
          )


setMethod("plot", "ContinuousProcess",
          function(x, nPoints = 200, ...){
            plotData <- getPlotData(object = x, nPoints = nPoints, ...)
            return(plot(plotData))
          }
          )

setMethod("getPosition", "ContinuousProcess",
          function(object, ...) {
            if(isTRUE(object@iSubset == -1L)) {
              value <- object@valueEnv$position
            } else {
              value <- object@valueEnv$position[iSubset(object)]
            }
            
            return(value)
           }
          )

setMethod("getTime", "ContinuousProcess",
          function(object, ...) {        
           return(getPosition(object, ...))
           }
          )

setMethod("showData", "ContinuousProcess",
          function(object, ...) {
            unitData <- getUnitData(object)
            id <- factor(getId(object))
            if(length(id) > 0) {
              splitEntries <- split(seq_along(id), id)
              ranges <- lapply(splitEntries, function(e) signif(range(getPosition(object)[e]), 3))
              if(is.list(ranges)) {
                positionSummary <- as.data.frame(sapply(ranges, function(r) paste("[", r[1],";",r[2],"]", sep="")), stringsAsFactors = FALSE) 
              } else {
                positionSummary <- data.frame(1)[FALSE,]
              }
              
              summaryById <- as.numeric(table(id))
              firstEntries <- sapply(splitEntries, function(e) e[1])
              values <- as.matrix(getValue(object))[firstEntries, , drop = FALSE]
              summaryById <- cbind(summaryById,
                                   positionSummary,
                                   unitData,
                                   as.data.frame(values))
              colnames(summaryById)[1:2] <- c("grid points", paste(object@positionVar,"range"))
            } else {
              summaryById <- unitData
            }
            
            idLevels <- levels(id)
            if(length(idLevels) > 4)
              idLevels <- c(idLevels[1:3], "...", idLevels[length(idLevels)])
            if(length(idLevels) > 1) {
              structure <- paste(length(levels(id)), " units.\n   ", object@idVar,": ", paste(idLevels, collapse = " "), "\n\n", sep="")
            } else {
              structure = ""
            }
                                             
            return(list(summary = summaryById, structure = structure))
          }
          )

setMethod("show", "ContinuousProcess",
          function(object) {
            d <- dim(object)
            op <- options()
            options(list(quote = FALSE, digits = 3))
            cat(paste("Class:", class(object), "\n"))
            cat(paste("Dimensions:", d[1], "x", d[2], "\n\n"))
            cat(showData(object)$structure)
            print(showData(object)$summary)
            options(op)
            return(invisible(object))
          }
          )

## TODO: Improve the summary function. Report ranges. 

setMethod("summary", "ContinuousProcess",
           function(object) {
             list("Summary of variables" =
                  summary(getUnitData(object)),
                  "Summary of values" =
                  apply(getValue(object), 2, summary))
           }
           )

setMethod("updateProcessObject", "ContinuousProcess",
          function(object, ...) {

            continuousProcess <- new("ContinuousProcess",
                                     unitData = object@unitData,
                                     metaData = object@metaData,
                                     iSubset = object@iSubset,
                                     jSubset = object@jSubset,
                                     idVar = object@idVar,
                                     colNames = object@colNames,
                                     positionVar = object@positionVar,
                                     valueEnv = object@valueEnv,
                                     factors = numeric())
            validObject(continuousProcess)
            return(continuousProcess)
            }
          )
