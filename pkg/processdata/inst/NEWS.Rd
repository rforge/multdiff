\name{NEWS}
\title{processdata News}

\section{Changes in processdata version 0.7.5}{
  \subsection{New features}{
    \itemize{
      \item A vignette, overview, has been written to introduce the
      usage of the package, its classes and methods. 
      \item A \code{getColumns} function has been implemented for column
      wise extraction of data, and other methods for data extraction are
      now based on this function. 
      \item The \code{getValue} function is deprecated and replaced by
      \code{getNumerics} and \code{getFactors}. A call of
      \code{getValue} is now identical to a call of \code{getNumerics}
      and returns a \code{matrix} instead of a \code{Matrix}.
      \item The \code{colNames} has a new argument, \code{type}, to
      extract column names for \code{'unit'}, \code{'factor'},
      \code{'numeric'}, \code{'mark'} and \code{'markValue'} types of
      columns only. 
      \item The \code{ContinuousProcess} class has been modified to hold
      columns of discrete data (termed factors) in addition to numeric data. 
      \item The \code{plot} method for a \code{ContinuousProcess}
      handles the factor columns in a sensible way. 
      \item More flexible constructors are now available for
      constructing \code{ContinuousProcess} and
      \code{MarkedPointProcess} objects from many combinations of data
      formats.
      \item A \code{coarsen} argument has been introduced for the
      constructors of \code{MarkedPointProcess} objects. This allows for
      coarsing event points to the observervation grid instead of adding
      additional points to the grid.
      \item \code{str} methods and \code{object.size} methods have been
      implemented.
      \item The \code{show} methods have been modified to print out only
      a limited amount of structural information, and the \code{summary}
      methods have taken over the role of computing and showing more
      (summarized) information on the content of the objects. 
    }
  }
   \subsection{Other Changes}{
    \itemize{
      \item The internal representations for all classes have been
      reimplemented and all objects from prior versions must be
      recreated from the original data.
      \item The implementation does not use the \pkg{Matrix} class
      anymore.
      \item The \code{plot} methods are now based on a formal S4 class instead
      of a list, and a \code{plot} method is implemented for this
      class. The class is, however, currently not exported. 
      \item Many parts of the code have been clarified and optimized,
      for instance the implementation of the \code{subset} method.
      \item A battery of unit tests based on \pkg{RUnit} has been added. 
    }
    }
}


