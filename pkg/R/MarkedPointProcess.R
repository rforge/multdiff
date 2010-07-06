setClass("MarkedPointProcess",
         representation(
                        id="factor",
                        position="numeric",
                        markType="factor",
                        markValue="data.frame"
                        ),
          validity=function(object){
           len <- length(object@id)
           if(len != length(object@position)) stop("mismatch in the size of the slots 'id' and 'position'")
           if(len != length(object@markType)) stop("mismatch in the size of the slots 'id' and 'markType'")
           if(dim(object@markValue)[1] > 0 & dim(object@markValue)[1] != len) stop("mismatch in the size of the slot 'id' and the dimension of 'markValue'")
           return(TRUE)
         }
         )

setMethod("initialize","MarkedPointProcess",
          function(.Object,
                   id=factor(),
                   position=numeric(),
                   markType=factor(),
                   markValue=data.frame()){
            .Object@id <- factor(id)
            if(length(id) == 0) .Object@id <- factor(rep(1,length(position))) 
            .Object@position <- position
            .Object@markType <- factor(markType)
            if(!is.data.frame(markValue)) stop("The 'markValue' needs to be a data frame")
            .Object@markValue <- markValue
            
            validObject(.Object)
            return(.Object)
          }
          )

setMethod("mpp2cp","MarkedPointProcess",
         function(object,position,process="countingProcess",...){
            id <- factor(rep(levels(object@id),sapply(position,length)))
            index <- lapply(position,function(x) seq(1,length(x)))
            aggIndex <- c(0,cumsum(sapply(index,length)))
            index <- lapply(seq(along=index),function(i) index[[i]]+aggIndex[[i]])
            if(length(levels(object@id)) != length(position)) stop("Not the same number of position vectors and id levels")
            countingProcess <- Matrix(0,nrow=aggIndex[length(aggIndex)],
                                      ncol=length(levels(object@markType)))
              
            
            for(i in seq(along=levels(object@markType))){
              for(j in seq(along=levels(object@id))){
                select <- (object@id == levels(object@id)[j]) & (object@markType == levels(object@markType)[i]) 
                cumPos <- sort(object@position[select])
                count <- apply(sapply(position[[j]], ">=",cumPos),2,sum)
                if(process=="countingProcess") countingProcess[index[[j]],i] <- count
                if(process=="pointProcess") countingProcess[index[[j]],i] <- c(0,diff(count))
              }
            }

            colnames(countingProcess) <- levels(object@markType)
            
            tmp <- new(Class="ContinuousProcess",
                       id=id,
                       position=unlist(position),
                       value=countingProcess)

            return(tmp)
          }
          )

setMethod("[",c("MarkedPointProcess","ANY","missing"),
          function(x,i,j,drop=FALSE){
           new("MarkedPointProcess",
               id=factor(x@id[i]),
               position=x@position[i],
               markType=factor(x@markType[i]),
               markValue=x@markValue[i,,drop=drop]) 
          }
          )

setMethod("getPlotData","MarkedPointProcess",
          function(object,...){
            tmpPlotData <- data.frame(id=factor(object@id),
                                      position=object@position,
                                      type=object@markType)
           
            return(tmpPlotData)
          }
          )

setMethod("plot",signature(x="MarkedPointProcess",y="missing"),
          function(x,y,...){
            plotData <- getPlotData(x,...)
            
            ggplot(data=plotData,aes(x=position,y=type,colour=type)) +
                   geom_point() +
                     scale_colour_discrete(breaks=levels(plotData$type)) +
                       scale_y_discrete(limits=levels(plotData$type)) + 
                         facet_grid(id ~ .)
          }
          )
                 

setMethod("getId","MarkedPointProcess",
          function(object,id=NULL){
            if(is.null(id)) select <- seq(along=object@id)  
            if(!is.null(id)) select <- object@id %in% id
            return(object@id[select])
          }
          )



setMethod("getPosition","MarkedPointProcess",
          function(object,id=NULL){
            if(is.null(id)) select <- seq(along=object@id)  
            if(!is.null(id)) select <- object@id %in% id
            return(object@position[select])
          }
          )


setMethod("getMarkType","MarkedPointProcess",
          function(object,id=NULL){
            if(is.null(id)) select <- seq(along=object@id)  
            if(!is.null(id)) select <- object@id %in% id
            return(object@markType[select])
          }
          )

setMethod("getMarkValue","MarkedPointProcess",
          function(object,id=NULL){
            if(dim(object@markValue)[1]>0){
              if(is.null(id)) select <- rep(TRUE,length.out=length(object@id))  
              if(!is.null(id)) select <- object@id %in% id
              result <- subset(object@markValue,select)
            } else {
              result <- c()
            }
            return(result)
          }
          )

setMethod("subset","MarkedPointProcess",
          function(x,subset,...){
            tmp <- data.frame(id=x@id,
                           position=x@position,
                           markType=x@markType)
            
            if(dim(x@markValue)[1] > 0) {
              tmp <- cbind(tmp,x@markValue)
            }

            if (missing(subset)) 
              r <- TRUE
            else {
              e <- substitute(subset)
              r <- eval(e, tmp, parent.frame())
              if (!is.logical(r)) 
                stop("'subset' must evaluate to logical")
              r <- r & !is.na(r)
            }
                      
            
            x@id <- factor(tmp$id[r])
            x@position <- tmp$position[r]
            x@markType <- factor(tmp$markType[r])
            if(dim(x@markValue)[1] > 0) {
              x@markValue <- tmp[r,-c(1,2,3),drop=FALSE]
            }
            return(x)
          }
          )
