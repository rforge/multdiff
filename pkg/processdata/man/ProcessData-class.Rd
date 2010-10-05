\name{ProcessData-class}
\docType{class}
\alias{ProcessData-class}
\alias{getUnitData,ProcessData-method}
\alias{setUnitData,ProcessData-method}
\alias{setUnitData}
\alias{colNames,ProcessData-method}
\alias{dim,ProcessData-method}

\title{Class "ProcessData"}
\description{	 

  A base class for holding unit specific informations and meta data
  for more complex classes for process data.
  
}
\section{Objects from the Class}{
  Objects of this class are usually not created. It is essentially just
  a data.frame wrapped up in an S4-object.

  The class serves the purpose of being a building block for the classes 
  \code{ContinuousProcess} and \code{MarkedPointProcess}.
}
\section{Slots}{
  \describe{
    
    \item{\code{metaData}:}{a \code{list}.}

    \item{\code{unitData}:}{a \code{data.frame}.}
  
  }
}
\section{Methods}{
  \describe{

    \item{colNames}{\code{signature(object = "ProcessData")}: Returns
       the column names.}
     
    \item{dim}{\code{signature(x = "ProcessData")}: Returns the
      dimensions.}

    \item{getUnitData}{\code{signature(object = "ProcessData")}: Returns
      the \code{data.frame} of unit data from slot \code{unitData}.}
    
    \item{setUnitData}{\code{signature(object = "ProcessData", value =
      "data.frame")}: Sets the unit data in slot \code{unitData} equal
      to \code{value}.} 

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
