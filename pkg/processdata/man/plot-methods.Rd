\name{plot-methods}
\alias{plot,ContinuousProcess-method}
\alias{plot,ContinuousProcess,ANY-method}
\alias{plot,JumpProcess-method}
\alias{plot,JumpProcess,missing-method}
\alias{plot,MarkedPointProcess-method}
\alias{plot,MarkedPointProcess,missing-method}
\alias{plot,MarkedPointProcess,character-method}
\alias{getPlotData,ContinuousProcess-method}
\alias{getPlotData,JumpProcess-method}
\alias{getPlotData,MarkedPointProcess-method}
\alias{getPlotData,MarkedPointProcess,character-method}
\alias{getPlotData}
\title{Methods for 'plot' in Package 'processdata'}
\description{
  Methods for plotting process data contained in objects of class
  \code{ContinuousProcess}, \code{MarkedPointProcess} or \code{JumpProcess}.
}
\usage{
\S4method{plot}{ContinuousProcess}(x, nPoints = 200, ...)

\S4method{plot}{MarkedPointProcess,character}(x, y, nPoints = 200, ...)

\S4method{plot}{JumpProcess}(x, nPoints = 200, ...)

\S4method{getPlotData}{ContinuousProcess}(object, nPoints = 200, allUnitData = FALSE, ...)

\S4method{getPlotData}{MarkedPointProcess,character}(object, y, nPoints = 200, allUnitData = FALSE, ...)

\S4method{getPlotData}{JumpProcess}(object, nPoints = 200, ...)
}
\arguments{

  \item{x}{ process object to be plotted.  }

  \item{y}{ a \code{character}. Specifies how the points are plotted for
  a \code{MarkedPointProcess} object. Possible values are '@mark',
  '@top', '@bottom', a name of a continuous process, a name of a mark
  value, the unit identification variable name or a numeric. See Details
  for further explanations. If missing, the default value is '@mark'.}
  
  \item{object}{ process object for which to compute a data set for
    plotting.  }

  \item{nPoints}{ a \code{numeric}. Specifies the number of observation
    points for each process to include in the data set that is
    plotted. Default value is 200.  }
  
  \item{allUnitData}{ a \code{logical}. Specifies if non-process data
    should be included in the data set for plotting. Default value is \code{FALSE}.  }

  \item{...}{ other arguments.}
}

\value{
  \code{plot} returns a \code{\link{ggplot}} object.

  \code{getPlotData} returns a \code{data.frame} for objects of class
  \code{ContinuousProcess}. For objects of class
  \code{MarkedPointProcess} and \code{JumpProcess} the value is a
  list with two components named
  \code{plotData} and \code{plotPointData}, each a \code{data.frame}. 
  
}

\details{

  The underlying plotting engine is \code{ggplot} from the \pkg{ggplot2}
  package. The result of a call to \code{plot} is thus a \code{ggplot}
  object, which can be altered and plotted like any other \code{ggplot}
  object.
  
  The default plot is a plot of the (multivariate) process data held in
  the process object \code{x} divided according to the observation units
  -- the individuals. This is most suitable with relatively few units.
  
  With many units and \code{p} the value of the call, the plot

  \code{p + facet_grid(. ~ .)}

  can be a more useful alternative. For \code{ContinuousProcess} objects
  one can also consider 

  \code{p + facet_grid(variable ~ .)}

  and for \code{MarkedPointProcess} objects

  \code{p + facet_grid(contVariable ~ .)}

  For \code{JumpProcess} objects one can in principle proceed as for
  \code{ContinuousProcess} objects. However, due to the discontinuous
  nature of a jump process, the curve for each unit becomes
  disconnected. A better alternative, assuming that the unit
  identification variable is \code{id}, is

    \code{p + facet_grid(variable ~ .) + aes(color = id)}
 
   See \code{demo(examples, package = "ppstat")} for examples.

   The \code{y} argument controls the y-values used for plotting
   the point process part of the data set when plotting a
   \code{MarkedPointProcess} object. Values '@top' and '@bottom' place
   all points above or below, respectively, of the plot of the
   \code{ContinuousProcess} part of the object. A numeric value of
   \code{y}, e.g. 0, places all points at this y-value. A name of a
   mark value places the points at their corresponding mark value and a
   name of a continuous process places the points on the continous
   process.  The default, '@mark', places the points at different levels
   on the y-axis corresponding to the different marks, and the unit
   identification variable name places the points at different levels on the
   y-axis corresponding to the different units. 
   
}

\seealso{

   \code{\linkS4class{ProcessData}}, \code{\linkS4class{ContinuousProcess}},
   \code{\linkS4class{MarkedPointProcess}},
   \code{\linkS4class{JumpProcess}}, 
   \code{\link{ggplot}}

 }

 \author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}
 
\keyword{methods}
\keyword{aplot}