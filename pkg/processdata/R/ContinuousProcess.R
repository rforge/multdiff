setClass("ContinuousProcess",
         representation(
                        id="factor",
                        position="numeric",
                        value="Matrix"
                        ),
         validity=function(object){
           len <- length(object@id)
           if((len != length(object@position)) || (len != dim(object@value)[1]))
             stop("mismatch in the size of slots")
           if(any(unlist(tapply(object@position,object@id,is.unsorted))))
             stop("Evaluation positions for continuous process must be sorted within each id")
           return(TRUE)
         }
         )


continuousProcess <- function(initdata){
  if(!is.data.frame(initdata)) stop("The 'initdata' needs to be a data frame")
  colNames <- names(initdata)
  if(!("id" %in% colNames)) id <- factor(rep(1,dim(initdata)[1])) else id <- initdata$id
  ind <- !(colNames %in% c("id","position"))                                  
  new("ContinuousProcess",id=id,position=initdata$position,value=initdata[,ind])
}

setMethod("initialize","ContinuousProcess",
          function(.Object,
                   id=factor(),
                   position=numeric(),
                   value=Matrix(ncol=0,nrow=0)){
            .Object@id <- factor(id)
            if(length(id) == 0) .Object@id <- factor(rep(1,length(position))) 
            .Object@position <- position
            if(is.data.frame(value)) {
              .Object@value <- Matrix(as.matrix(value),dimnames=dimnames(value))
            } else {
              .Object@value <- Matrix(value)
            }             
            validObject(.Object)
            return(.Object)
          }
          )
            

setMethod("show","ContinuousProcess",
          function(object){
            tmp <- list(id=summary(object@id),
                        positions=list(summary=summary(object@position),
                          head=head(object@position)),
                        processes=list(columnwiseSummary=apply(object@value,2,summary),
                          head=head(object@value)))
            print(tmp)
          }
          )

setMethod("[",c("ContinuousProcess","ANY","missing"),
          function(x,i,j,drop=FALSE){
           new("ContinuousProcess",
               id=factor(x@id[i]),
               position=x@position[i],
               value=x@value[i,]) 
          }
          )

setMethod("getPlotData","ContinuousProcess",
          function(object,by=1,...){
            select <- seq(1,length(object@id),by=by)        
            tmp <- as.matrix(object@value[select,])
            row.names(tmp) <- NULL
            tmpPlotData <- cbind(id=factor(object@id[select]),
                                 position=object@position[select],
                                 melt(as.data.frame(tmp),id.vars=NULL))
            return(tmpPlotData)
          }
          )

setMethod("plot",signature(x="ContinuousProcess",y="missing"),
          function(x,y,by=1,...){

            ggplot(data=getPlotData(x,by=by,...),
                   aes(x=position,y=value,colour=variable)) +
                     facet_grid(id ~ .) + geom_line()
          }
          )

setMethod("getValue","ContinuousProcess",
          function(object,...){
            return(object@value)
          }
          )
            
setMethod("getPosition","ContinuousProcess",
          function(object,...){
              return(object@position)
          }
          )

setMethod("getId","ContinuousProcess",
          function(object,...){
            return(object@id)
          }
          )

setMethod("subset","ContinuousProcess",
          function(x,subset,...){
            if (missing(subset)) 
              r <- TRUE
            else {
              e <- substitute(subset)
              r <- eval(e, data.frame(id=x@id,position=x@position), parent.frame())
              if (!is.logical(r)) 
                stop("'subset' must evaluate to logical")
              r <- r & !is.na(r)
            }
            x@id <- factor(x@id[r])
            x@position <- x@position[r]
            x@value <- x@value[r,,drop=FALSE]
            return(x)
          }
          )

