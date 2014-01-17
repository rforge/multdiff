\name{condMeanVar}
\alias{condMeanVar}
\alias{condMeanVar,OUModel-method}
\title{Conditional Mean and Variance}
\description{
  \code{condMeanVar} is a method which returns the conditional mean and, optionally, the covariance matrix of
      \eqn{X_{t_0+t}|X_{t_0} = x_{t_0}}{X_{t_0+t} given X_{t_0} =
	x_{t_0}} for \eqn{N}{N} pairs of \eqn{x_{t_0}}{x_{t_0}} and
      \eqn{t}{t}.
    }
    \usage{
\S4method{condMeanVar}{OUModel}(object, parameters,
	x, t, var = FALSE, ...)
}

\arguments{

  \item{object}{ An object of class \code{OUModel}. }

  \item{parameters}{ A \code{list} of parameters which are compatible
    with \code{object}. If \code{parameters} are not specified, the parameters from \code{object} are used.}

  \item{x}{  If \code{x} is
    specified, \code{t} must be specified as well (and vice
    versa).
    
    \code{x} must be a \eqn{N\times p}{N x p} \code{matrix}, with rows
    representing \eqn{x_{t_0}}{x_{t_0}}-vectors.
    
  }

  \item{t}{ If \code{t} is
    specified, \code{x} must be specified as well (and vice
    versa).

    \code{t} may be a
    \code{numeric} vector of length \eqn{N}{N} with \eqn{t}{t}-values
    corresponding to the rows in \code{x}. \code{t} may also be a \code{numeric} (of length
    1), in which case the same value of \eqn{t}{t} is used for all
    calculations.
    }

  \item{var}{ If \code{var = FALSE} only the
    conditional means are computed. If \code{var = TRUE} both the
    conditional means and conditional covariance matrices are computed.}

  \item{...}{ Other arguments. }
}

\details{
  If \code{x} and \code{t} are missing, the data and time-distance
    between observations from \code{object} are used, and the
    conditional means and, optionally, the conditional covariance
    matrices of \eqn{X_{t_2}|X_{t_1} = x_{t_1}}{X_{t_2}|X_{t_1} =
      x_{t_1}} up to \eqn{X_{t_n}|X_{t_{(n-1)}} =
      x_{t_{(n-1)}}}{X_{t_n}|X_{t_(n-1)} = x_{t_(n-1)}} are computed.
    }

  \value{
    If \code{var = FALSE}, the conditional means are returned as columns
    in a \eqn{p\times N}{p x N} \code{matrix}.

    If \code{var = TRUE}, the conditional means and covariance matrices
    are returned in a list with two elements. The first element,
    \code{"condMean"}, contains the conditional mean vectors as columns
    in a \eqn{p\times N}{p x N} \code{matrix}. The second element,
    \code{"condVar"}, contains the conditional covariance matrices in a
    \eqn{p\times p \times N}{p x p x N} \code{array}. 
  }

  \author{Nina Munkholt, \email{nina.m@math.ku.dk}}
  
\seealso{
  \code{\linkS4class{MultDiffModel}}, \code{\linkS4class{OUModel}}
}
