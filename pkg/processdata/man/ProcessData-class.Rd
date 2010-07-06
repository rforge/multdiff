\name{ProcessData-class}
\docType{class}
\alias{ProcessData-class}
\alias{getContinuousProcess,ProcessData-method}
\alias{getContinuousProcess}
\alias{getMarkTypePosition,ProcessData-method}
\alias{getMarkTypePosition}
\alias{getMarkedPointProcess,ProcessData-method}
\alias{getMarkedPointProcess}
\alias{getPointPosition,ProcessData-method}
\alias{getPointPosition}
\alias{initialize,ProcessData-method}
\alias{plot,ProcessData,missing-method}
\alias{subset,ProcessData-method}

\title{Class "ProcessData"}
\description{	 

An object of class \code{ProcessData} combines continuous 
process data and marked point process data into a single, 
coherent object. Note that whether this is meaningful is up to
the user -- it assumes a common time/position axis. One important
gain is consistent subsetting behavior across the two different
data types. 

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ProcessData", markedPointProcess, continuousProcess)}.

The object consists of two parts. One object of class \code{ContinuousProcess} 
and one object of class \code{MarkedPointProcess}


}
\section{Slots}{
	 \describe{
    \item{\code{markedPointProcess}:}{Object of class \code{"MarkedPointProcess"}}
    \item{\code{continuousProcess}:}{Object of class \code{"ContinuousProcess"}}
    \item{\code{pointPosition}:}{Object of class \code{"data.frame"} - contains the positions of the point process and the corresponding discrete marks.}
  }
}
\section{Methods}{
  \describe{
    \item{getContinuousProcess}{\code{signature(object = "ProcessData")}: Returns the continuous process.}
    \item{getMarkTypePosition}{\code{signature(object = "ProcessData")}: Returns the point process positions for a specific discrete mark value.}
    \item{getMarkedPointProcess}{\code{signature(object = "ProcessData")}: Returns the marked point process.}
    \item{getPointPosition}{\code{signature(object = "ProcessData")}: Returns all the (time) points of observations.}
    \item{initialize}{\code{signature(.Object = "ProcessData")}.}
    \item{plot}{\code{signature(x = "ProcessData", y = "missing")}: Produces a combined plot of the continuous process(es) 
and the marked point process.}
    \item{subset}{\code{signature(x = "ProcessData")}: Selects a subset of the data by passing 
the arguments to the subset method for marked point process, see
\code{\linkS4class{MarkedPointProcess}}, and the subset method for
continuous processes, see \code{\linkS4class{ContinuousProcess}}. One of
two logical parameters, \code{mppOnly} and \code{cpOnly}, can be
specified, which controls if the subset method should be invoked for
the marked point process or the continuous process only. Default is that
\code{mppOnly=FALSE} and \code{cpOnly=FALSE}.
}
	 }
       }
       
\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
	 \code{\linkS4class{ContinuousProcess}} and  \code{\linkS4class{MarkedPointProcess}}. 
}
\examples{
showClass("ProcessData")
}
\keyword{classes}
