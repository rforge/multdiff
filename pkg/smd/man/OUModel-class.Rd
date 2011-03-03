\name{OUModel-class}
\docType{class}
\alias{OUModel-class}
\alias{condMeanVar}
\alias{condMeanVar,OUModel-method}
\alias{validateParameters}
\alias{validateParameters,OUModel-method}
\alias{parToList}
\alias{parToList,OUModel-method}
\alias{parToNumeric}
\alias{parToNumeric,OUModel-method}
\alias{gradient}
\alias{gradient,OUModel-method}

%% Explain alpha, beta, Sigma

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
    \eqn{dX_t = B(X_t-A) dt + D dW_t} where \eqn{\mathbf{W}}{W} is a p-dimensional Brownian
    Motion and \eqn{D} satisfies \eqn{DD^\textrm{T} = C}{D*t(D) = C}.
    }
    
    
\section{Slots}{
  \describe{
    \item{\code{parameters}:}{Object of class \code{"list"} which contains
      components \code{A}, \code{B} and \code{C} of class
      \code{"matrix"}, corresponding to the parameters of the model.}
    \item{\code{data}:}{Object of class \code{"ContinuousProcess"}
      containing the data used to estimate the diffusion model.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{MultDiffModel}"} (directly).
}
\section{Methods}{
  \describe{
    \item{condMeanVar}{\code{signature(object = "OUModel")}: Computes
      the conditional mean and covariance matrix of \eqn{X_{t_0+t} given
	X_{t_0} = x}{X_{t_0+t} given X_{t_0} = x}.}
    \item{validateParameters}{\code{signature(object = "OUModel")}:
  May be used to check that a set of parameters corresponding to an
  object of class \code{"OUModel"} has the correct form.}
\item{parToList}{\code{signature(object = "OUModel")}: Converts a vector
  into a list of parameters (A,B,C) (by column)}
\item{parToNumeric}{\code{signature(object = "OUModel")}: Converts a
  list of parameters into a vector.}
  }
}
%\references{}
%\usage{
%\S4method{condMeanVar}{OUModel}(object, parameters, x, t, var=FALSE)
%
%If \code{"parameters"} are not specified they are always taken to be the parameters of
%\code{"object"}. \code{"var"} indicates whether or not the conditional variance should
%be calculated.
%
%\S4method{validateParameters}{OUModel}(object, parameters)
%\S4method{parToList}{OUModel}(object, x)
%\S4method{parToNumeric}{OUModel}(object, parameters)
%
%}
\author{Nina Munkholt, \email{nina.m@math.ku.dk}}
\seealso{
  \code{\linkS4class{MultDiffModel}} 
}
\examples{

showClass("OUModel")

}
\keyword{classes}
