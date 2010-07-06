\name{processdata-package}
\alias{processdata-package}
\alias{processdata}
\docType{package}
\title{Process Data}
\description{

  This package implements data structures for stochastic processes,
  which consists of either discretely observed continuous
  time processes, marked point processes or a combination.

}
\details{
\tabular{ll}{
Package: \tab processdata\cr
Type: \tab Package\cr
Version: \tab 0.5 \cr
Date: \tab 2010-07-06 \cr
License: \tab GPL version 2 or newer \cr
LazyLoad: \tab yes\cr
}

The class \code{\linkS4class{ContinuousProcess}} implements a data
structure for observations from a potentially multivariate continuous
time process observed at discrete time points, which are not necessarily
equidistant. The class \code{\linkS4class{MarkedPointProcess}}
implements a data structure for observations of a marked point process.
The class \code{\linkS4class{ProcessData}} combines the two classes
above for combined observations of continuous time processes and marked
point processes in one data structure. 


}
\author{
  Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}.
  
  Maintainer: \email{Niels.R.Hansen@math.ku.dk}
}


\examples{

### Setting up a simple continuous time process with equidistant 
### sampled observations.

exampleCP <- new("ContinuousProcess",
                 position = seq(0,1,0.01),
                 value = Matrix(rnorm(202,c(rep(1,101),rep(-1,101))),ncol=2))

\dontrun{plot(exampleCP)}

### Setting up a simple point process with two marks

examplePP <- new("MarkedPointProcess",
                 position = c(0.1,0.12,0.35,0.5,0.55,0.56,0.7,0.75,0.9,0.91),
                 markType = c("A","A","B","B","A","A","B","A","B","A"))

\dontrun{plot(examplePP)}

### Combining the two data structures in a single object

exampleProcess <- new(Class="ProcessData",
                      markedPointProcess=examplePP,
                      continuousProcess=exampleCP)

\dontrun{plot(exampleProcess)}


}
\keyword{models}
