setClass("ProcessData",
         representation(
                        ## metaData is a list, in general unrestricted.
                        ## The intention is that metaData contains information
                        ## on the experiment or procedure by which the
                        ## data is obtained.
                        metaData = "list",
                        
                        ## unitData holds general information for each unit
                        ## for which we have observations.
                        unitData = "data.frame"
                        ),
         validity = function(object) {
           
           return(TRUE)
         }
         ## validity = function(object) {
         ##   if(!("id" %in% names(object@unitData)))
         ##     stop("The 'id' column for the 'unitData' data frame in 'ProcessData' is not defined")
         ##   if(class(object@unitData$id) !=  "factor")
         ##     stop("The 'id' column for the 'unitData' data frame in 'ProcessData' must be a 'factor'")
           
         ##   return(TRUE)
         ## }
         )

setClass("ContinuousProcess",
         representation(
                        ## Two vectors indexing a subset of the full data set
                        ## in the valueEnv environment.
                        iSubset = "integer",
                        jSubset = "integer",
                        
                        ## Names of the id variable and position variable.
                        idVar = "character",
                        positionVar = "character",
                        
                        ## Column names. 
                        colNames = "character",                      

                        ## The 'valueEnv' contains the following three
                        ## components of the data:
                        ## id = "factor",
                        ## position = "numeric",
                        ## value = "Matrix",
                        ## i = "numeric",
                        ## j = "numeric"
                        valueEnv = "environment"                 
                        ),
         contains = "ProcessData",
         validity = function(object){
           ## Checks that the environment contains the correct components.
           if(!is.factor(object@valueEnv$id))
             stop(paste("No", object@idVar, "or", object@idVar, "is not a factor."))
           if(!is.numeric(object@valueEnv$position))
             stop(paste("No", object@positionVar, "or", object@positionVar, "is not a numeric."))
           if(!is(object@valueEnv$value, "Matrix"))
             stop("No 'value' or 'value' is not a Matrix.")

           ## Checks that assumed orderings are correct.
           if(any(levels(getId(object)) != rownames(getUnitData(object))) || length(levels(getId(object))) != dim(getUnitData(object))[1]){
             stop(paste("The levels of", object@idVar, "and the row names of unitData are either not of equal length or in the same order."))
           }
           if(any(unlist(tapply(getPosition(object),getId(object),is.unsorted))))
             stop(paste(object@positionVar,"not sorted within", object@idVar))
            
           ## Checks that the components in the environment have the correct
           ## format. 
           len <- length(object@valueEnv$id)
           if((len != length(object@valueEnv$position)) || (len != dim(object@valueEnv$value)[1]))
             stop("Mismatch in the size of slots for 'ContinuousProcess'.")
          
           return(TRUE)
         }
         )


setClass("MarkedPointProcess",
         representation(
                        iPointSubset = "integer",
                        jPointSubset = "integer",
                        markVar = "character",
                        pointPointer = "integer",
                        
                        ## The 'pointProcessEnv' environment contains
                        ## the following four components:
                        ## id = "factor",
                        ## position = "numeric",
                        ## markType = "factor",
                        ## markValue = "data.frame"
                        pointProcessEnv = "environment"
                        ),
         contains = "ContinuousProcess",
         validity = function(object) {
           len <- length(object@pointProcessEnv$id)
           if(len != length(object@pointProcessEnv$position))
             stop("Sizes of the slots 'id' and 'position' do not match.")
           if(len != length(object@pointProcessEnv$markType))
             stop("Sizes of the slots 'id' and 'markType' do not match.")
           if(dim(object@pointProcessEnv$markValue)[1] > 0 & dim(object@pointProcessEnv$markValue)[1] != len)
             stop("Size of the slot 'id' and the dimension of 'markValue' do not match.")
           if(any(levels(getPointId(object)) != rownames(getUnitData(object)))) 
             stop(paste("The point process levels of", object@idVar, "and the row names of unitData are not of equel length or in the same order."))
           if(any(unlist(tapply(getPointPosition(object),getPointId(object),is.unsorted))))
             stop(paste(object@positionVar,"for the point process data not sorted within", object@idVar))
           
           return(TRUE)
         }
         )

setClass("JumpProcess",
         representation(jumpVar = "character"),
         contains = "MarkedPointProcess",
         validity = function(object) {
           return(TRUE)
         }
         )

## TODO: Validity checks, constructors etc. 
