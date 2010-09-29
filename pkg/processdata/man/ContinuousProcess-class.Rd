\name{ContinuousProcess-class}
\docType{class}
\alias{ContinuousProcess-class}
\alias{dim,ContinuousProcess-method}
\alias{colNames,ContinuousProcess-method}
\alias{colNames<-,ContinuousProcess,character-method}
\alias{getId,ContinuousProcess-method}
\alias{getId}
\alias{getPosition}
\alias{getPosition,ContinuousProcess-method}
\alias{getTime}
\alias{getTime,ContinuousProcess-method}
\alias{getValue,ContinuousProcess-method}
\alias{getValue}
\alias{getUnitData,ContinuousProcess-method}
\alias{getUnitData}
\alias{show,ontinuousProcess-method}

\title{Class "ContinuousProcess"}
\description{	
  
An object of class \code{ContinuousProcess} is a container for data
obtained as discrete observations of univariate or multivariate continuous 
time processes. 

}
\section{Objects from the Class}{
  Objects can be created by calls of the function \code{continuousProcess}.
  
  The object contains time points or positions in one-dimensional space 
  and corresponding observations of a univariate or multivariate
  process. The data structure can hold data for several units, or
  individuals.  

  The class extends \code{ProcessData} and is extended by the classes
  \code{markedPointProcess} and \code{JumpProcess}. These classes make
  up an infrastructure for storing and handling general stochastic
  process data. 

}
\section{Slots}{
  \describe{
    \item{\code{colNames}:}{ a \code{character}.}

    \item{\code{idVar}:}{ a \code{character}.}

    \item{\code{iSubset}:}{ an \code{integer} vector. References for subsetting.}
    
    \item{\code{jSubset}:}{ an \code{integer} vector. References for subsetting.}

    \item{\code{positionVar}:}{ a \code{character}.}
    
    \item{\code{valueEnv}:}{ an \code{environment}.}
  }
}
\section{Methods}{
  \describe{

     \item{colNames}{\code{signature(object = "ContinuousProcess")}: Returns
       the column names.}

     \item{colNames<-}{\code{signature(object = "ContinuousProcess",
       value="character")}: Sets the column names.}
     
     \item{dim}{\code{signature(x = "ContinuousProcess")}: Returns the
       dimensions.}
    
    \item{getId}{\code{signature(object = "ContinuousProcess")}: Returns
      the identification variable with name in slot \code{idVar}. }

    \item{getPosition}{\code{signature(object = "ContinuousProcess")}: Returns the time/position of the observations. }

    \item{getTime}{\code{signature(object = "ContinuousProcess")}: Same
    as \code{getPosition}.}
    
    \item{getValue}{\code{signature(object = "ContinuousProcess")}:
      Returns the observations as a \code{"Matrix"}.}

     \item{getUnitData}{\code{signature(object = "ContinuousProcess")}:
      Returns the unit data as a \code{"data.frame"} with column names
      the unit identifikations.}
    
  }

    Other methods are documented in their own help pages.
}

\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
  \code{\linkS4class{ProcessData}},  \code{\link{continuousProcess}}. 

}

\examples{
showClass("ContinuousProcess")
CP <- continuousProcess(data.frame(id = rep(c(1,2), each = 10),
                                   time = rep(seq(0.1,1,0.1), 2),
                                   value = rnorm(20)))  
subset(CP, id == 1)
subset(CP, time  >= 0.5)
}
\keyword{classes}
