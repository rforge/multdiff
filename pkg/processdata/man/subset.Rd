\name{subset}
\title{subset}
\alias{subset,ProcessData-method}
\alias{[,ProcessData-method}
\alias{[,ProcessData,ANY,ANY,ANY-method}
\alias{[,ProcessData,ANY,missing,ANY-method}
\alias{[,ProcessData,missing,ANY,ANY-method}
\alias{subset,ContinuousProcess-method}
\alias{[,ContinuousProcess-method}
\alias{[,ContinuousProcess,ANY,ANY,ANY-method}
\alias{[,ContinuousProcess,missing,missing,missing-method}
\alias{[,ContinuousProcess,integer,missing,ANY-method}
\alias{[,ContinuousProcess,numeric,missing,ANY-method}
\alias{[,ContinuousProcess,logical,missing,ANY-method}
\alias{[,ContinuousProcess,missing,integer,ANY-method}
\alias{[,ContinuousProcess,missing,logical,ANY-method}
\alias{[,ContinuousProcess,missing,character,ANY-method}
\alias{[,ContinuousProcess,missing,numeric,ANY-method}
\alias{subset,MarkedPointProcess-method}
\alias{[,MarkedPointProcess-method}
\alias{[,MarkedPointProcess,ANY,ANY,ANY-method}
\alias{[,MarkedPointProcess,missing,missing,missing-method}
\alias{[,MarkedPointProcess,integer,missing,ANY-method}
\alias{[,MarkedPointProcess,numeric,missing,ANY-method}
\alias{[,MarkedPointProcess,logical,missing,ANY-method}
\alias{[,MarkedPointProcess,missing,integer,ANY-method}
\alias{[,MarkedPointProcess,missing,logical,ANY-method}
\alias{[,MarkedPointProcess,missing,character,ANY-method}
\alias{[,MarkedPointProcess,missing,numeric,ANY-method}
\alias{subset,JumpProcess-method}
\alias{[,JumpProcess-method}
\alias{[,JumpProcess,ANY,ANY,ANY-method}

\description{
  Methods for subsetting and using \code{subset} or \code{[} for process data objects.
}

\usage{
\S4method{[}{ProcessData}(x, i, j, ..., drop = FALSE)

\S4method{[}{ContinuousProcess}(x, i, j, ..., drop = FALSE)

\S4method{[}{MarkedPointProcess}(x, i, j, ..., drop = FALSE)

\S4method{[}{JumpProcess}(x, i, j, ..., drop = FALSE)

\S4method{subset}{ProcessData}(x, subset, select, ...)

\S4method{subset}{ContinuousProcess}(x, subset, select, ...)

\S4method{subset}{MarkedPointProcess}(x, subset, select, markSubset, ...)

\S4method{subset}{JumpProcess}(x, subset, select, ...)


}

\arguments{

  \item{x}{ a process object.}

  \item{i}{ an \code{integer}, \code{numeric}, \code{logical} or
  \code{character}.}

  \item{j}{ an \code{integer}, \code{numeric}, \code{logical} or
  \code{character}.}

  \item{drop}{ a \code{logical}. Detault value \code{FALSE}.}

  \item{subset}{ an expression given in terms of the variables in the
     data set that evaluates to \code{logical}. See details for further
     explanations.}

   \item{select}{ an expression that evaluates to \code{logical} in a
    frame consisting of the column names.}

   \item{markSubset}{ an expression given in terms of the mark
   variables in the data set that evaluates to \code{logical}. See
   details for further explanations.  }

   \item{...}{other arguments.}
}

\details{

TO BE FILLED IN.     

}

\value{

  A process object of the same class as \code{x} unless \code{drop =
  TRUE}. If \code{drop = TRUE} the object is dropped to one further down
  the class hierarchy if appropriate. For instance, subsetting to column entries
  for a \code{MarkedPointProcess} such that all mark columns are lost
  result in an object of class \code{ContinuousProcess} if  \code{drop =
    TRUE}.
 
}

\author{
  Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}
  }

\seealso{
\code{\linkS4class{ContinuousProcess}},
\code{\linkS4class{MarkedPointProcess}}, \code{\linkS4class{JumpProcess}}  

}

\examples{
data(example)
contExam[1:1500, ]
subset(pointExam, time > 4 & gender == "M")

contExam[ ,2:3]

colNames(pointExam)
pointExam[ , c("gender", "foo", "ALPHA")]
pointExam[ , c("gender", "foo"), drop = TRUE]
subset(pointExam, select = gender:BETA)




}
\keyword{manip}
