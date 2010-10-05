\name{ExampleData}
\alias{example}
\alias{contExam}
\alias{pointExam}
\alias{jumpExam}


\docType{data}
\title{
Three small example data sets.
}
\description{

  The data sets are small, simulated data sets used to illustrate the
  functionality of the \code{processdata} data structure, subsetting and
  plotting methods.

}
\usage{example}

\format{ The three data sets \code{contExam}, \code{pointExam},
  \code{jumpExam} in \code{example} are of class \code{\linkS4class{ContinuousProcess}},
  \code{\linkS4class{MarkedPointProcess}} and
  \code{\linkS4class{JumpProcess}}, respectively. 
 
}
\seealso{\code{\link{exampleDataFrames}}}
\examples{
data(example)
contExam
pointExam

}
\keyword{datasets}