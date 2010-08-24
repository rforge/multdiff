\name{MultDiffModel-class}
\docType{class}
\alias{MultDiffModel-class}
\alias{getData}
\alias{getData,MultDiffModel-method}
\alias{getParameters}
\alias{getParameters,MultDiffModel-method}
\alias{loss}
\alias{loss,MultDiffModel-method}


\title{Class "MultDiffModel"}
\description{

  \code{MultDiffModel} is a virtual class which represents a
  multivariate diffusion model.
  
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"ContinuousProcess"}}
    \item{\code{parameters}:}{Object of class \code{"list"}}
  }
}
\section{Methods}{
  \describe{
    \item{getData}{\code{signature(object = "MultDiffModel")}: Returns
      the data object of class \code{"ContinuousProcess"} used to
      estimate the model. }
    \item{getParameters}{\code{signature(object = "MultDiffModel")}:
      Returns the parameters of the model. }
    \item{loss}{\code{signature(object = "MultDiffModel")}:
      Computes the loss   }
	 }
}
\author{
  Nina Munkholt, \email{nina.m@math.ku.dk}
  
}
\seealso{
  \code{\linkS4class{OUModel}} 
}
\examples{
showClass("MultDiffModel")
}
\keyword{classes}
