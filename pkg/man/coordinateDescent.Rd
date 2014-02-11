\name{coordinateDescent}
\alias{coordinateDescent}
\title{l1-penalized least squares coordinate descent}
\usage{
coordinateDescent(y, xi, dxi, p, beta, lambda = NULL, nlambda = 50,
  lambda.min.ratio = 0.001, penalty.factor = rep(1, p), rho = 0.9,
  c = 1e-04, reltol = sqrt(.Machine$double.eps), trace = 0, N = 1000)
}
\arguments{
  \item{y}{a \code{numeric}. The response vector.}

  \item{xi}{a \code{function}. The mean value map. A vector
  function of \code{beta}.}

  \item{dxi}{a \code{function}. The derivative of the mean
  value map. A function of \code{beta} that returns a
  matrix.}

  \item{p}{a \code{numeric}. The length of the parameter
  vector.}

  \item{beta}{a \code{numeric}. The vector of initial
  parameter values. If \code{p} is missing \code{beta} must
  be specified. When \code{p} is specified, the value of
  the \code{beta} argument is ignored, and \code{beta =
  rep(0, p)}.}

  \item{lambda}{a \code{numeric}. A sequence of penalties.
  The default value \code{NULL} implies an automatic
  computation of a suitable sequence.}

  \item{nlambda}{a \code{numeric}. The number of lambda
  values.}

  \item{lambda.min.ratio}{a \code{numeric}. Determines the
  smallest value of lambda relative to the largest
  (automatically) determined value of lambda.}

  \item{penalty.factor}{a \code{numeric}. A vector of
  penalty weight factors.}

  \item{rho}{a \code{numeric}. Step length control in the
  backtracking.}

  \item{c}{a \code{numeric}. Sufficient decrease control in
  the backtracking.}

  \item{reltol}{a \code{numeric}. Controls the convergence
  criterion.}

  \item{trace}{a \code{numeric}. Values above 0 prints
  increasing amounts of trace information.}

  \item{N}{a \code{numeric}. The maximal number of
  interations in the loops.}
}
\value{
A \code{list} of length 3. The first element is the
sequence \code{lambda}, and the second, \code{beta}, is the
matrix of parameter estimates. Each column in \code{beta}
corresponds to an entry in \code{lambda}. The third element
is \code{status}, where 0 means convergence, and 1 means
termination due to maximal number of interations reached
(for some lambda).
}
\description{
Computes l1-penalized least squares estimates using a
coordinate wise descent algorithm of Gauss-Newton quadratic
approximations.
}
\details{
The function computes a matrix of parameter estimates. Each
column corresponds to a value of the penalty parameter in
\code{lambda}. The estimates are computed in decreasing
order of the penalty parameters, and for each column the
previous is used as a warm start.

The algorithm relies on iterative optimization of the
l1-penalized Gauss-Newton quadratic approximation using a
standard coordinate wise descent algorithm. An outer
backtracking step is added to ensure that the algorithm
takes descent steps.

This function relies on two auxiliary functions. The
functions \code{xi} and \code{dxi}, which are the mean
value map and its derivative, respectively.

The function returns a list with the vector of lambda
values as the first entry and the estimated beta parameters
as a matrix in the second entry. Each column in the matrix
corresponds to one lambda value.
}
\author{
Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}
}
\seealso{
\code{\link{coordinateDescentMF}}.
}

