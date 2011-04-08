\name{OUModel-class}
\docType{class}
\alias{OUModel-class}
\alias{gradient}
\alias{gradient,OUModel-method}
\alias{parToList}
\alias{parToList,OUModel-method}
\alias{parToNumeric}
\alias{parToNumeric,OUModel-method}
\alias{validateParameters}
\alias{validateParameters,OUModel-method}



\title{Class "OUModel"}
\description{

    The class \code{OUModel} represents multivariate
    Ornstein-Uhlenbeck diffusion models. 

  }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("OUModel",
    ...)}. 
  
  A \eqn{p}-dimensional Ornstein-Uhlenbeck diffusion model is determined
  by three parameters,
  \describe{
      \item{\eqn{A}:}{A \eqn{p}{p}-vector (the constant drift coefficient).}
      \item{\eqn{B}:}{A \eqn{p\times p}{p x p} matrix (the linear drift
	coefficient).}
      \item{\eqn{C}:}{A symmetric, positive semidefinite \eqn{p\times
	  p}{p x p} matrix (the infinitesimal covariance matrix).}
    }

    The diffusion is a solution to a SDE of the form
    \eqn{dX_t = B(X_t-A) dt + D dW_t} where \eqn{\mathbf{W}}{W} is a p-dimensional Brownian
    Motion and \eqn{D} satisfies \eqn{DD^T = C}{D*t(D) = C}.

    The class extends \code{MultDiffModel}, a virtual class which
    represents multivariate diffusion models.
    }

\section{Slots}{
  \describe{
    
    \item{\code{parameters}:}{ A \code{list} which contains
      components \code{A}, \code{B} and \code{C} corresponding to the
      parameters of the model. \code{A} must be a \code{numeric} vector,
      \code{B} and \code{C} must be of class \code{matrix}.}
    
  }
}

\section{Extends}{
Class \code{\linkS4class{MultDiffModel}} (directly).
}

\section{Methods}{
  \describe{

    \item{gradient}{
      \code{signature(object = "OUModel", parameters = NULL, lossType =
	1, useSufficientStat = TRUE, ...)}:
      Returns the gradient \eqn{\nabla l(B,A)}{\nabla l(B,A)} of the
      type 1 or 2 loss function \eqn{l}{l} as a \code{list} with two elements
      \code{A} and \code{B} (the derivatives with respect to the
      parameters \code{A} and \code{B}, respectively). Uses data from \code{object}. If
      \code{parameters} are missing, the parameters of
      \code{object} are used. If \code{parameters} are not
      missing, they must be specified as a \code{list} of parameters which
      are compatible with \code{object}. \code{lossType} must be either
      \code{1} or \code{2}.
      }

    \item{validateParameters}{\code{signature(object = "OUModel", parameters)}:
      May be used to check that \code{parameters}, a set of parameters corresponding to an
      object of class \code{OUModel}, has the correct form.}

    \item{parToList}{\code{signature(object = "OUModel", x)}: Converts the
      \code{numeric} vector \code{x} into a \code{list} of parameters
      \code{A}, \code{B} and \code{C} which are compatible with an
      object of class \code{OUModel} (by column).}

    \item{parToNumeric}{\code{signature(object = "OUModel",
	parameters)}: Converts \code{parameters}, a \code{list} of parameters
      corresponding to an object of  class \code{OUModel}, into a
      \code{numeric} vector. If \code{parameters} are not specified,
      the parameters of \code{object} are converted.}
  }
}

\author{
  Nina Munkholt, \email{nina.m@math.ku.dk}
}
\seealso{
  \code{\linkS4class{MultDiffModel}} 
}
\examples{

showClass("OUModel")

}

\keyword{classes}
