\name{ExampleData}
\alias{contExam}
\alias{pointExam}
\alias{jumpExam}
\alias{unitData}
\alias{continuousData}
\alias{pointData}

\docType{data}
\title{
Three small example data sets.
}
\description{

  The data sets are small, simulated data sets used to illustrate the
  functionality of the \code{processdata} data structure, subsetting and
  plotting methods.

}
\usage{

data(example, package = "processdata")
data(exampleDataFrames, package = "processdata")

}
\format{ The three data sets \code{contExam}, \code{pointExam},
  \code{jumpExam} in \code{example} are of class \code{\linkS4class{ContinuousProcess}},
  \code{\linkS4class{MarkedPointProcess}} and
  \code{\linkS4class{JumpProcess}}, respectively. The three data frames
  \code{unitData}, \code{continuousData} and \code{pointData} contains
  the data for the \code{contExam} and \code{pointExam}, and can be used
  for the contruction of these.
 
}

\examples{
data(example)
contExam
pointExam


data(exampleDataFrames)
continuousProcess(continuousData, unitData)
markedPointProcess(pointData, continuousData, unitData = unitData)


}
\keyword{datasets}