\name{OUModel-class}
\docType{class}
\alias{OUModel-class}

\title{Class "OUModel"}
\description{

  \code{OUModel} is a class which represents multivariate
  Ornstein-Uhlenbeck diffusion models. %% Explain A, B, D...
  
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OUModel", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{A}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{B}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{D}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{data}:}{Object of class \code{"ContinuousProcess"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{MultDiffModel}"}, directly.
}
\section{Methods}{
No methods defined with class "OUModel" in the signature.
}
%%\references{}
\author{
Nina Munkholt, \email{nina.m@math.ku.dk}
}
\examples{
showClass("OUModel")
}
\keyword{classes}
