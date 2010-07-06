\name{ContinuousProcess-class}
\docType{class}
\alias{ContinuousProcess-class}
\alias{getId,ContinuousProcess-method}
\alias{getPlotData,ContinuousProcess-method}
\alias{getPosition,ContinuousProcess-method}
\alias{getValue,ContinuousProcess-method}
\alias{getValue}
\alias{initialize,ContinuousProcess-method}
\alias{plot,ContinuousProcess,missing-method}
\alias{show,ContinuousProcess-method}
\alias{subset,ContinuousProcess-method}

\title{Class "ContinuousProcess"}
\description{	
  
An object of class \code{ContinuousProcess} is a container for data
obtained as discrete observations of univariate or multivariate continuous 
time processes. 

}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ContinuousProcess", id, position, value)}.
  
  The object contains time points or positions in one-dimensional space 
  and corresponding observations of a univariate or multivariate process. 

}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"factor"} - the
      identification of the continuous process that the observation belongs to.}
    \item{\code{position}:}{Object of class \code{"numeric"} - the
      position or time of the observation. }
    \item{\code{value}:}{Object of class \code{"Matrix"} - the value(s) of the observations 
      stored as a (potentially) sparse matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{getId}{\code{signature(object = "ContinuousProcess")}: Returns
      the identification variable \code{id}. }
    \item{getPlotData}{\code{signature(object = "ContinuousProcess")}: Returns the data frame that is used for plotting. }
    \item{getPosition}{\code{signature(object = "ContinuousProcess")}: Returns the time or position, \code{position}, of the observations. }
    \item{getValue}{\code{signature(object = "ContinuousProcess")}:
      Returns the observations, \code{value}, as a \code{"Matrix"}.}
    \item{initialize}{\code{signature(.Object = "ContinuousProcess")}.}
    \item{plot}{\code{signature(x = "ContinuousProcess", y =
	"missing")}: Plots the observations as continuous process(es).}
    \item{show}{\code{signature(object = "ContinuousProcess")}.}
    \item{subset}{\code{signature(x = "ContinuousProcess")}: Selects a
      subset of the full data set by either specifying a vector of
      logicals or an expression in terms of \code{id}, \code{position}
      or \code{value} that will evaluate to logical. Returns an
      object of class \code{ContinuousProcess}. }
  }
}

\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
  \code{\linkS4class{ProcessData}} 
}

\examples{
showClass("ContinuousProcess")
}
\keyword{classes}
