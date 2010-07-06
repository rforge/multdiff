\name{continuousProcess}
\Rdversion{1.1}
\alias{continuousProcess}
\title{Continuous Process
}
\description{
  Creates continuous process of type \code{ContinuousProcess} from a
  data frame.
}
\usage{
continuousProcess(initdata)
}
\arguments{
  \item{initdata}{
  A data frame holding a column named \code{position} and optionally one
  named \code{id}. Additional columns must be numeric and are taken as
  the observations of the process.
  }
}

\value{
An object of class \code{ContinuousProcess}.
}

\author{
  Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}
  }

\seealso{
\code{\linkS4class{ContinuousProcess}}  
}
\examples{
CP <- data.frame(position=1:10,value=rnorm(10))
continuousProcess(CP)
}
\keyword{manip}
