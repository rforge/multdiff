\name{continuousProcess}
\alias{continuousProcess}
\alias{continuousProcess,data.frame-method}
\title{Continuous Process}
\description{
  Creates a \code{ContinuousProcess} object from 
  data frame(s).
}
\usage{
\S4method{continuousProcess}{data.frame}(continuousData, unitData = data.frame(), metaData = list(), positionVar = 'time', idVar = 'id', ...)
}
\arguments{
  
  \item{continuousData}{ a data frame containing the process data. An
    optional column with name \code{positionVar} (default "time")
    specifies the x-coordinates where we have observations, and an
    optional column with name \code{idVar} (default "id") specifies the
    identification variable for the division of the data into units.
    Additional columns must be numeric and are taken as the coordinates
    of observations of the process.}

  \item{unitData}{an optional data frame containing further unit
    specific variables.}

  \item{metaData}{an unstructured \code{list} containing additional
    information on the data in the resulting \code{ContinuousProcess}
    object.}

  \item{positionVar}{a \code{character} specifying the name of the
    variable holding the x-coordinates for the observations. Default
    value \code{'time'}.}

  \item{idVar}{a \code{character} specifying the name of the
    identification variable for each unit. Default value 'id'.}

  \item{...}{ other arguments.}

}

\value{ An object of class \code{ContinuousProcess}.  }

\details{

  If the data frame \code{continuousData} lacks either of the
  \code{positionVar} or \code{idVar} columns they are
  automatically generated. If \code{positionVar} is missing the
  observations are assumed equidistant at time points 1, 2 etc within
  each unit. If \code{idVar} is missing it is assumed that there is only
  one unit. 
  
}

\author{
  Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}
  }

\seealso{
\code{\linkS4class{ContinuousProcess}}  
}
\examples{

CP <- data.frame(time = 1:10, value = rnorm(10))
continuousProcess(CP)

CP <- data.frame(id = c(rep("A",5), rep("B", 5)), value = rnorm(10))
continuousProcess(CP)

}
\keyword{manip}
