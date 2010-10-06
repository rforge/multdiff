setMethod("jumpProcess", "MarkedPointProcess",
          function(object, ...) {
            .Object <- new("JumpProcess")
            as(.Object, "MarkedPointProcess") <- object
            .Object@jumpVar <- "jump"
            markLevels <- levels(getMarkType(.Object))
            jumpLevels <- markLevels %in% colnames(getValue(.Object))
            valueOrder <- colnames(getValue(.Object))[colnames(getValue(.Object)) %in% markLevels[jumpLevels]]
            if(length(jumpLevels) > 0 && any(markLevels[jumpLevels] != valueOrder)) {
              pointProcessEnv <- new.env(parent = .GlobalEnv)
              pointProcessEnv$id <- getPointId(.Object)
              pointProcessEnv$position <- getPointPosition(.Object)
              pointProcessEnv$markType <- factor(getMarkType(.Object), levels = c(valueOrder, markLevels[!jumpLevels]))
              pointProcessEnv$markValue <- getMarkValue(.Object)
                                                 
              .Object@pointProcessEnv <- pointProcessEnv
            }

            .Object@colNames <- unique(colNames(.Object))
            validObject(.Object)
            return(.Object)
          }
          )

setMethod("jumpProcess", "data.frame",
          function(object, continuousData, ...) {
            callGeneric(markedPointProcess(object, continuousData), ...)
          }
          )


setMethod("getPlotData", "JumpProcess",
          function(object, nPoints = 200, ...) {
            plotData <- getPlotData(as(object, "MarkedPointProcess"), y = object@jumpVar, nPoints = nPoints, ...)
            marks <- getMarkType(object)
            plotPointData <- plotData$plotPointData
            plotData <- plotData$plotData
            jumpVariables <- plotPointData$variable %in% levels(plotData$variable)
            plotPointData <- plotPointData[jumpVariables, ]
            pointee <- getPointPointer(object)[jumpVariables]
            variable <- as.character(plotPointData$variable)
            plotPointData$value <-  as.numeric(sapply(seq_along(plotPointData$id), function(i) getValue(object)[pointee[i], variable[i]]))
            plotPointData$iSubset <- pointee
            plotPointDataLeft <- plotPointData
            plotPointDataLeft$iSubset <- plotPointDataLeft$iSubset - 1
            if(object@jumpVar %in% colNames(object)) {
              plotPointDataLeft$value <- plotPointDataLeft$value -  as.numeric(plotPointDataLeft[ ,object@jumpVar])
            }
            plotData <- rbind(plotData,
                              plotPointData[ , colnames(plotData), drop = FALSE],
                              plotPointDataLeft[ , colnames(plotData), drop = FALSE])
            plotData <- plotData[order(plotData$variable),]            
            continuousComp <- list()
            for(var in levels(plotData$variable)) {
              continuousComp[[var]] <- sapply(plotData$iSubset[plotData$variable == var], function(i) sum(pointee[marks == var] <= i))
            }
            plotData$continuousComp <- as.factor(unlist(continuousComp))
            plotPointData$continuousComp <- plotData$continuousComp[plotPointData$iSubset]
            plotPointDataLeft$continuousComp <- plotData$continuousComp[plotPointDataLeft$iSubset]
            plotPointData <- list(Right = plotPointData,
                             Left = plotPointDataLeft)
           
                         
            return(list(plotData = plotData, plotPointData = plotPointData))
          }
          )
              
          

setMethod("plot", c("JumpProcess", "missing"),
          function(x, nPoints = 200, ...){
            plotData <- getPlotData(x, nPoints = nPoints, ...)

            group <- paste(x@idVar, ":continuousComp:variable", sep = "")

            
            p <-  ggplot(data = plotData$plotData, aes_string(x = "position", y = "value", colour = "variable", group = group)) +
              facet_grid(id ~ ., scales = "free_y") +
                scale_x_continuous(name = x@positionVar) +
                  scale_y_continuous(name = "") +
                        geom_line()

            if(dim(plotData$plotPointData$Right)[1] > 0) {
              p <- p +    geom_point(data = plotData$plotPointData$Right, shape = 19, size = 1.5, aes(group = NULL)) +
                geom_point(data = plotData$plotPointData$Left, shape = 1, size = 1.5, aes(group = NULL))
            }
            
            return(p)
          }
          )

setMethod("showData", "JumpProcess",
          function(object, ...) {
            summaryById <- callNextMethod()
            structure <- strsplit(summaryById$structure, "\n")[[1]]
           if (length(structure) > 3) 
             structure <- paste(structure[-(length(structure) - c(0, 1, 2))], "\n")
            summaryById <- summaryById$summary

            return(list(summary = summaryById, structure = structure))            
          }
          )
          
