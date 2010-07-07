\name{MarkedPointProcess-class}
\docType{class}
\alias{MarkedPointProcess-class}
\alias{getId,MarkedPointProcess-method}
\alias{getId}
\alias{getMarkType,MarkedPointProcess-method}
\alias{getMarkType}
\alias{getMarkValue,MarkedPointProcess-method}
\alias{getMarkValue}
\alias{getPlotData,MarkedPointProcess-method}
\alias{getPlotData}
\alias{getPosition,MarkedPointProcess-method}
\alias{getPosition}
\alias{mpp2cp,MarkedPointProcess-method}
\alias{mpp2cp}
\alias{plot,MarkedPointProcess,missing-method}
\alias{subset,MarkedPointProcess-method}

\title{Class "MarkedPointProcess"}
\description{	

An object of class \code{MarkedPointProcess} is a container for data
obtained as observations of a one-dimensional marked point process. 

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MarkedPointProcess", ...)}.

The object contains time points or positions in one-dimensional space of
discrete events and their associated mark information. Discrete marks are given a
special treatment but any type of mark value is possible.  


}
\section{Slots}{
	 \describe{
    \item{\code{id}:}{Object of class \code{"factor"} - the
      identification of the point process that the observation belongs to.}
    \item{\code{position}:}{Object of class \code{"numeric"} - the
      position or time of the observation.}
    \item{\code{markType}:}{Object of class \code{"factor"} - the
      (discrete) mark of the observation.}
    \item{\code{markValue}:}{Object of class \code{"data.frame"} - any
      addition mark values of the observation.}
  }
}
\section{Methods}{
  \describe{
    \item{getId}{\code{signature(object = "MarkedPointProcess")}:
      Returns the identification variable \code{id}.}
    \item{getMarkType}{\code{signature(object = "MarkedPointProcess")}:
      Returns the mark variable \code{markType}.}
    \item{getMarkValue}{\code{signature(object = "MarkedPointProcess")}:
      Returns the additional mark values \code{markValue}.}
    \item{getPlotData}{\code{signature(object = "MarkedPointProcess")}:
      Returns the data frame that is used for plotting.}
    \item{getPosition}{\code{signature(object = "MarkedPointProcess")}:
      Returns the time or position, \code{position}, of the observed discrete events.}
    \item{mpp2cp}{\code{signature(object = "MarkedPointProcess")}:
      Converts the object of class \code{MarkedPointProcess} to an
      object of class \code{\linkS4class{ContinuousProcess}}.}
    \item{plot}{\code{signature(x = "MarkedPointProcess", y =
	"missing",usedLayer=geom_point())}: Plots the data.}
    \item{subset}{\code{signature(x = "MarkedPointProcess")}: Selects a
      subset of the full data set by either specifying a vector of
      logicals or an expression in terms of \code{id}, \code{position}
      or \code{markType} that will evaluate to logical. Returns an
      object of class \code{MarkedPointProcess}.
	 }
       }
       }
\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
  \code{\linkS4class{ProcessData}} 
}

\examples{
showClass("MarkedPointProcess")
}
\keyword{classes}
