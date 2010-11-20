setMethod("plot", c("ProcessPlotData", "missing"),
          function(x, y, ...) {

            ### Setting up y-values for discrete variables
            value <- factor()
            factorpoint <- numeric()
            
            if(dim(x@factorPlotData)[1] > 0){
              value <- x@factorPlotData$value
              factorpoint <- rep(1, length(value))
            }

            if(dim(x@pointPlotData)[1] > 0){
              if(is.factor(x@pointPlotData$value)){
                levels <- unique(c(levels(value), levels(x@pointPlotData$value)))
                value <- factor(c(as.character(value), as.character(x@pointPlotData$value)), levels = levels)
                factorpoint <- c(factorpoint, rep(2, length(x@pointPlotData$value)))
              }
            }

            factorLevels <- levels(value)
            factorRange <-  length(factorLevels)*(x@limits[2] - x@limits[1])/(length(factorLevels) + 3)
            offset <- x@limits[2]
            if(x@position == "bottom") {
              factorRange <- -factorRange
              offset <- x@limits[1]
            }
            
            x@labels <- c(x@labels, factorLevels)
            x@breaks <- c(x@breaks, offset + factorRange*seq_along(factorLevels)/(length(factorLevels) + 1))

            if(factorpoint[1] == 1){
              x@factorPlotData$value <- offset + factorRange*as.numeric(value[factorpoint == 1])/(length(factorLevels) + 1)
            }

            if(factorpoint[length(factorpoint)] == 2){
              x@pointPlotData$value <- offset + factorRange*as.numeric(value[factorpoint == 2])/(length(factorLevels) + 1)
            }
            
            facetFormula <- as.formula(paste(x@idVar, "~ ."))

            p <- ggplot(x@continuousPlotData, aes(x = position)) +
              facet_grid(facetFormula) +
                scale_x_continuous(x@positionVar) +
                  scale_y_continuous(breaks = x@breaks,
                                     name = "",
                                     labels = x@labels)
            
            if("value" %in% names(x@continuousPlotData)){
              group <- paste(x@idVar, ":variable", sep = "")
              p <- p + geom_line(aes_string(y = "value",
                                            colour = "variable",
                                            group = group))
            }

            if(dim(x@factorPlotData)[1] > 0){
              p <- p + geom_line(data = x@factorPlotData,
                                 aes_string(y = "value",
                                            colour = "variable",
                                            group = "group"),
                                 size = 3)
            }

             if(dim(x@pointPlotData)[1] > 0){
              p <- p + geom_point(data = x@pointPlotData,
                                  aes_string(y = "value",
                                             colour = "variable"))
            }

            return(p)
          }
          )
