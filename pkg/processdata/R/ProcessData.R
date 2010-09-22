setMethod("dim", "ProcessData",
          function(x) {
            return(dim(x@unitData))
          }
          )
                     
setMethod("colNames", "ProcessData",
          function(object, ...) {
            return(names(object@unitData))
          }
          )

setMethod("getUnitData", "ProcessData",
          function(object, ...) {
            return(object@unitData)
          }
          )

setMethod("[", "ProcessData",
          function(x, i, j, ..., drop = FALSE) {
            x@unitData <- x@unitData[i, j, drop = drop]
            return(x)
          }
          )

setMethod("[", c(x = "ProcessData", i = "ANY", j = "missing"),
          function(x, i, j, ..., drop = FALSE) {
            j <- seq_len(dim(as(x, "ProcessData"))[2])
            x <- callGeneric(x, i, j, drop)
            return(x)
          }
          )

setMethod("[", c(x = "ProcessData", i = "missing", j = "ANY"),
          function(x, i, j, ..., drop = FALSE) {
            i <- seq_len(dim(as(x, "ProcessData"))[1])
            x <- callGeneric(x, i, j, drop)
            return(x)
          }
          )

setMethod("subset", "ProcessData",
          function(x, subset, select, ...) {
            if (missing(subset)) 
              r <- TRUE
            else {
              e <- substitute(subset)
              r <- eval(e, x@unitData, parent.frame())
              if (!is.logical(r)) 
                stop("'subset' must evaluate to logical")
              r <- r & !is.na(r)
            }
            if (missing(select)) 
              vars <- TRUE
            else {
              nl <- as.list(seq(1L,dim(x@unitData)[2]))
              names(nl) <- colnames(x@unitData)
              vars <- eval(substitute(select), nl, parent.frame())
            }
            
            return(x[r, vars])
          }
          )

          
setReplaceMethod("setUnitData", c("ProcessData", "data.frame"),
          function(object, value) {
            object@unitData <- value
            validObject(object)
            return(object)
          }
          )
