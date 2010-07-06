setClass("ProcessData",
         representation(
                        markedPointProcess = "MarkedPointProcess",
                        continuousProcess = "ContinuousProcess",
                        pointPosition = "data.frame"
           )
         )

setMethod("initialize","ProcessData",
          function(.Object,
                   markedPointProcess=new("MarkedPointProcess"),
                   continuousProcess=new("ContinuousProcess")){
            ## Compute a combined set of evaluation points and keep track of which
            ## points are points from the marked point process.
            
            evalPosition <- data.frame(id = factor(c(as.character(getId(markedPointProcess)),
                                         as.character(getId(continuousProcess))),
                                         levels=levels(getId(markedPointProcess))),
                                       position = c(getPosition(markedPointProcess),
                                         getPosition(continuousProcess)),
                                       pointer=c(seq(along=getId(markedPointProcess)),
                                         rep(Inf,length(getId(continuousProcess)))))
            
            dup <- unlist(tapply(evalPosition$position,evalPosition$id,duplicated),use.names=FALSE)
            evalPosition <- evalPosition[c(rep(TRUE,length(getId(markedPointProcess))),!dup[-c(seq(along=getId(markedPointProcess)))]),]

            ord <- tapply(evalPosition$position,evalPosition$id,order)
            ord <- unlist(lapply(levels(evalPosition$id),
                                 function(id) which(evalPosition$id==id)[ord[[id]]]),use.names=FALSE)

             where <- which(evalPosition$pointer[ord] < Inf)
            .Object@pointPosition <- data.frame(pointer=evalPosition$pointer[ord][where],where=where)

            j <- sapply(levels(evalPosition$id),
                               function(i) {
                                 findInterval(evalPosition[ord,"position"][evalPosition[ord,"id"]==i],
                                              c(getPosition(continuousProcess)[getId(continuousProcess)==i],Inf),
                                              all.inside=TRUE)},
                        simplify=FALSE)

            j <- unlist(lapply(levels(evalPosition$id),
                        function(i) {
                          which(getId(continuousProcess)==i)[j[[i]]]
                        }),
                        use.names=FALSE)
           

            value <- Matrix(as.matrix(getValue(continuousProcess))[j,,drop=FALSE])
           .Object@continuousProcess <- new("ContinuousProcess",id=evalPosition[ord,"id"],position=evalPosition[ord,"position"],value=value)

           
            .Object@markedPointProcess <- markedPointProcess
            
            return(.Object)
          }
          )

setMethod("[",c("ProcessData","ANY","missing"),
          function(x,i,j,drop=FALSE){
            cp <- x@continuousProcess[i]
            mpp <- x@markedPointProcess[x@pointPosition$pointer[x@pointPosition$where %in% i]]
            new("ProcessData",
                markedPointProcess=mpp,
                continuousProcess=cp)
          }
          )

            
setMethod("getContinuousProcess","ProcessData",
          function(object,pos){
            return(object@continuousProcess)
          }
          )

setMethod("getMarkedPointProcess","ProcessData",
          function(object){
            return(object@markedPointProcess)
          }
          )

setMethod("getPointPosition","ProcessData",
          function(object){
            return(object@pointPosition)
          }
          )

setMethod("getMarkTypePosition","ProcessData",
          function(object,type){
            return(object@pointPosition$where[object@markedPointProcess@markType[object@pointPosition$pointer]==type])
          }
          )           


setMethod("plot",signature(x="ProcessData",y="missing"),
          function(x,y,by=1,...){
            if(dim(getValue(x@continuousProcess))[2]==0){
              returnPlot <- plot(x@markedPointProcess) 
            } else {
              plotMarkedData <- getPlotData(x@markedPointProcess)
              plotContinuousData <- getPlotData(x@continuousProcess,by=by)
              
              plotMarkedData$processType <- factor(rep("points",length.out=dim(plotMarkedData)[1]),
                                                   levels=c("points","continuous"))
              plotContinuousData$processType <- factor(rep("continuous",length.out=dim(plotContinuousData)[1]),
                                                       levels=c("points","continuous"))

             # plotData <- as.data.frame(table(rbind(plotMarkedData[,c("id","processType")],plotContinuousData[,c("id","processType")])))

              low <- signif(min(plotContinuousData$value),1)
              high <- signif(max(plotContinuousData$value),1)
              prettyBreaks <- pretty(c(low,high))

              numtype <- as.numeric(plotMarkedData$type)
              breaks <- c(prettyBreaks,high + seq(1,max(numtype))*(high-low)/max(numtype))
              numtype <- numtype*(high-low)/max(numtype)
              plotMarkedData$numtype <- numtype + high
             
              returnPlot <-  ggplot() +
                facet_grid(id ~ .,scales="free_y")  +
                  geom_point(data=plotMarkedData,aes(x=position,y=numtype,colour=type)) +
                    geom_line(data=plotContinuousData,aes(x=position,y=value,colour=variable)) +
                      scale_y_continuous(breaks=breaks,name="",labels=c(as.character(prettyBreaks),levels(plotMarkedData$type))) +
                        scale_colour_discrete("Variable",breaks=c(levels(plotMarkedData$type),levels(plotContinuousData$variable)))
            }
            return(returnPlot)
          }
          )

## Possible arguments to the following function are subset, mppOnly and cpOnly

setMethod("subset","ProcessData",
          function(x,...,mppOnly=FALSE,cpOnly=FALSE){
            continuousProcess <- if(!mppOnly) subset(x@continuousProcess,...) else x@continuousProcess
            markedPointProcess <- if(!cpOnly) subset(x@markedPointProcess,...) else x@markedPointProcess

            return(new("ProcessData",markedPointProcess,continuousProcess))
          }
          )
            
            
          
