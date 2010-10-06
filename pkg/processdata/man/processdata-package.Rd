\name{processdata-package}
\alias{processdata-package}
\alias{processdata}
\docType{package}
\title{Process Data}
\description{

  This package implements data structures for stochastic processes.
  The data can be in the form of either discretely observed continuous
  time processes, marked point processes, jump processes or a combination.

}
\details{
\tabular{ll}{
Package: \tab processdata\cr
Type: \tab Package\cr
Version: \tab 0.7 \cr
Date: \tab 2010-10-07 \cr
License: \tab GPL version 2 or newer \cr
LazyLoad: \tab yes\cr
}

The class \code{\linkS4class{ContinuousProcess}} implements a data
structure for observations from a potentially multivariate continuous
time process observed at discrete time points, which are not necessarily
equidistant. The class \code{\linkS4class{MarkedPointProcess}} extends
this class to include observations of marked point processes. A further
extension is the class \code{\linkS4class{JumpProcess}}, which uses the
combination of marked points and continuous process observations to
represent data from jump processes.



}
\author{
  Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}.
  
  Maintainer: \email{Niels.R.Hansen@math.ku.dk}
}


\examples{

### Setting up a simple continuous time process with equidistant 
### sampled observations.

CPframe <- data.frame(time = seq(0,1,0.01),
                      val1 = rnorm(101,log(1:100)),
                      val2 = cumsum(rnorm(101,-0.1)))

CP <- continuousProcess(CPframe)

\dontrun{
plot(CP)}

### Setting up a simple point process with two marks

PPframe <- data.frame(
 time = c(0.1,0.12,0.35,0.5,0.55,0.56,0.7,0.75,0.9,0.91),
 markType = c("A","A","B","B","A","A","B","A","B","A")
)

PP <- markedPointProcess(PPframe, seq(0,1,0.01))


\dontrun{
plot(PP)}

### Combining the two data structures in a single object

process <- markedPointProcess(PPframe, CPframe)

\dontrun{
plot(process)}


}
\keyword{models}
