\name{MultDiffModel-class}
\docType{class}
\alias{MultDiffModel-class}
\alias{getData}
\alias{getData,MultDiffModel-method}
\alias{getParameters}
\alias{getParameters,MultDiffModel-method}
\alias{getSufficientStat}
\alias{getSufficientStat,MultDiffmodel-method}
\alias{loss}
\alias{loss,MultDiffModel-method}


\title{Class "MultDiffModel"}
\description{

  \code{MultDiffModel} is a virtual class which represents a
  multivariate diffusion model.
  
}
\section{Objects from the class}{
  A virtual class: No objects may be created from it.

  The class is extended by the class \code{OUModel}, the objects of
  which represent multivariate Ornstein-Uhlenbeck diffusion models.
  
}
\section{Slots}{
  \describe{
    
    \item{\code{data}:}{ An object of class \code{ContinuousProcess}
      containing the data used to estimate the model and additional
      information about the structure of the data.}
    
    \item{\code{parameters}:}{ A \code{list} of the parameters of the model. }

    \item{\code{sufficientStat}:}{ An \code{environment}. }

    \item{\code{hasSufficientStat}:}{ An object of type
      \code{logical}. Flag to indicate if sufficient statistics are computed. }
  }
}
\section{Methods}{
  \describe{
    
    \item{getData}{
      \code{signature(object = "MultDiffModel", ...)}:
      Returns the slot \code{data} as an object of class
      \code{"ContinuousProcess"}.
    }
    
    \item{getParameters}{
      \code{signature(object = "MultDiffModel", ...)}:
      Returns the slot \cr \code{parameters} as a \code{list}.
    }

    \item{getSufficientStat}{
      \code{signature(object = "MultDiffModel", stat, ...)}:
      ### Description missing ###
      }
    
    \item{loss}{
      \code{signature(object = "MultDiffModel", parameters = NULL,
	lossType = 1, ...)}: Returns the value of the
      type 1 or 2 loss function. Uses \code{data} from \code{object}. If
      \code{parameters} are missing, the calculations are based on
      \code{parameters} from \code{object}. If \code{parameters} are not
      missing, they must be specified as a \code{list} of parameters which
      are compatible with \code{object}. \code{losstype} must be either
      \code{1} or \code{2}.
    }
    
  }
 
}

\author{Nina Munkholt, \email{nina.m@math.ku.dk}, Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
  \code{\linkS4class{OUModel}} 
}

\examples{
showClass("MultDiffModel")
}

\keyword{classes}
