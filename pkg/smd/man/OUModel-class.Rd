\name{OUModel-class}
\docType{class}
\alias{OUModel-class}
\alias{condMeanVar}
\alias{condMeanVar,OUModel-method}

\title{Class "OUModel"}
\description{

    The class "OUModel" represents multivariate
    Ornstein-Uhlenbeck diffusion models.

  }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("OUModel",
    ...)}.
  
  A \eqn{p}-dimensional Ornstein-Uhlenbeck diffusion model is determined
  by three parameters,
  \describe{
      \item{\eqn{A}:}{A \eqn{p\times
	  1}{p x 1} matrix (the constant drift coefficient).}
      \item{\eqn{B}:}{A \eqn{p\times p}{p x p} matrix (the linear drift
	coefficient).}
      \item{\eqn{C}:}{A symmetric, positive semidefinite \eqn{p\times
	  p}{p x p} matrix (the infinitesimal covariance matrix).}
    }

    The diffusion is a solution to a SDE of the form
    \deqn{dX_t = (A + BX_t) dt + D dW_t} where \eqn{\mathbf{W}}{W} is a p-dimensional Brownian
    Motion and \eqn{D} satisfies \eqn{DD^\textrm{T} = C}{D*t(D) = C}.
    }
    
    
\section{Slots}{
  \describe{
    \item{\code{parameters}:}{Object of class \code{"list"} with
      components \code{A}, \code{B} and \code{C} of class
      \code{"matrix"}, corresponding to the parameters of the model.}
    \item{\code{data}:}{Object of class \code{"ContinuousProcess"}
      containing the data used to estimate the diffusion model.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{MultDiffModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{condMeanVar}{\code{signature(object = "OUModel")}: Computes
      the conditional mean and covariance matrix of the model.}}
}
%\references{}
\author{Nina Munkholt, \email{nina.m@math.ku.dk}}
\examples{
showClass("OUModel")
}
\keyword{classes}
