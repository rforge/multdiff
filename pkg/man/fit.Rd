\name{fit}
\alias{fit}
\alias{fit.multModel}
\title{Fitting a multivariate process model}
\usage{
fit(model, ...)

\method{fit}{multModel}(model, method = "cd", ...)
}
\arguments{
  \item{model}{}

  \item{method}{a \code{character}. The default,
  \code{"cd"} means a standard coordinate descent
  algorithm.  Alternative, \code{'mf'} means a matrix free
  method.}
}
\value{
a \code{multModel} object with the additional elements
\tabular{ll}{ lambda \tab the lambda sequence used \cr
Blambda \tab the corresponding matrix of estimated
parameters \cr status \tab either 0 (meaning no errors), 1
(convergence), 2 (MLE not computed) or 3 (other errors)
indicating different problems or errors \cr msg \tab errors
encountered or other notable conditions }
}
\description{
Estimation of parameters in a multivariate process model
using l1-penalized squared error loss and a coordinate wise
descent algorithm.
}
\seealso{
\code{\link{multModel}}, \code{\link{coordinateDescent}},
\code{\link{coordinateDescentMF}}
}

