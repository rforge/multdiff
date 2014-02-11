\name{riskHat}
\alias{riskHat}
\title{Estimation of the risk.}
\usage{
riskHat(model, par, df, sigma, ...)
}
\arguments{
  \item{model.}{The model object for which the risk is to
  be computed.}

  \item{par}{a \code{numeric}. The estimated parameter
  vector.}

  \item{df}{a \code{numeric}. The effective degrees of
  freedom.}

  \item{sigma}{a \code{numeric}. The variance parameter.}
}
\value{
a numeric.
}
\description{
Computation of an estimate of the risk based on the
effective degrees of freedom. The estimate can be used for
model selection.
}
\author{
Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}
}
\seealso{
\code{\link{df}}, \code{\link{risk}}
}

