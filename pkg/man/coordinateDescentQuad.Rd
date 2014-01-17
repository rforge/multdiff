\name{coordinateDescentQuad}
\alias{coordinateDescentQuad}
\title{Optimize l1-penalized loss using a coordinate descent algorithm}
\usage{
coordinateDescentQuad(beta, f, gr, quad, lambda = 1, penalty.factor = rep(1,
  p), rho = 0.9, c = 1e-04, reltol = sqrt(.Machine$double.eps),
  trace = 0)
}
\arguments{
  \item{beta}{a \code{numeric}. The vector of initial
  parameter values}

  \item{f}{a \code{function}. The loss function. A function
  of \code{beta}}

  \item{gr}{a \code{function}. The gradient of the loss. A
  function of \code{beta}}

  \item{quad}{a \code{function}. The coordinate wise
  quadratic approximation term.  A function of the
  coordinate index and \code{beta}}

  \item{lambda}{a \code{numeric}. A sequence of penalties.
  Default value \code{NULL} implies a automatic computation
  of a suitable sequence}

  \item{penalty.factor}{a \code{numeric}. A vector of
  weight factors}

  \item{rho}{a \code{numeric}. Step length control in the
  backtracking}

  \item{c}{a \code{numeric}. Sufficient decrease control in
  the backtracking}

  \item{reltol}{a \code{numeric}. Controls the convergence
  criterion}

  \item{trace}{a \code{numeric}. Values above 0 prints
  increasing amounts of trace information.}
}
\value{
A \code{list} of length 2. The first entry contains
\code{lambda} and the second the matrix of parameter
estimates. Each column in the matrix corresponds to an
entry in \code{lambda}.
}
\description{
The function computes a matrix of parameter estimates. Each
column corresponds to a value of the penalty parameter in
\code{lambda}. The estimates are computed in decreasing
order of the penalty parameters, and for each column the
previous is used as a warm start.
}
\details{
The algorithm relies on iterative optimization of an
l1-penalized quadratic approximation of the loss using a
standard coordinate descent algorithm. An outer
backtracking step is added to ensure that the algorithm
takes descent steps.

This function relies on three auxiliary functions. The loss
function \code{f}, its gradient \code{gr} and a third
function, \code{quad}, that computes the coefficient of the
quadratic approximation for the coordinate wise
optimization.

The function returns a list with the vector of lambda
values as the first entry and the estimated beta parameters
as a matrix in the second entry. Each column
}

