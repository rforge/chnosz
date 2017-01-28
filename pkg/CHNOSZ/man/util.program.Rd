\name{util.program}
\alias{util.program}
\alias{caller.name}
\alias{palply}

\title{Functions for Various Programming Tasks}
\description{
  Get name of calling function; send messages to standard output; use multiple processors if parallel package is loaded.
}

\usage{
  caller.name(n)
  palply(varlist, X, FUN, ...)
}

\arguments{
  \item{n}{numeric, number of frame to go up}
  \item{...}{equivalent to the same argument in \code{\link[parallel]{parLapply}}}
  \item{varlist}{character, names of variables to export using \code{\link[parallel]{clusterExport}}}
  \item{X}{vector, argument for \code{\link{lapply}} or \code{parLapply}}
  \item{FUN}{function, argument for \code{lapply} or \code{parLapply}}
}

\details{

\code{caller.name} returns the name of the calling function \code{n} frames up (i.e., for \code{n} equal to 2, the caller of the function that calls this one).
If called interactively, returns \code{\link{character}()}.

\code{palply} is a wrapper function to run \code{\link[parallel]{parLapply}} if length of \code{X} > \code{\link{thermo}$opt$paramin} and multiple cores are available, otherwise it runs \code{\link{lapply}}.
Note that \code{parLapply} is called with \code{methods} set to FALSE.
If lots of package startup messages are created when running \code{\link[parallel]{makeCluster}} (which is called by \code{palply}), it can probably be stopped by adding a test for \code{\link{interactive}} sessions around any \code{\link{library}} commands in the \code{\link{Rprofile}}.

}

\seealso{
\code{\link{read.fasta}}, \code{\link{count.aa}}, \code{\link{affinity}}, \code{\link{equil.boltzmann}} and \code{\link{equil.reaction}} for functions that use \code{palply}. Tests are in \file{tests/test-util.program.R}, and a \dQuote{real life} application is in \file{demos/density.R}.
}


\examples{
caller.name()  # character(0)
afun <- function() caller.name()
afun()         # character(0)
bfun <- function() afun()
bfun()         # "bfun"
}

\keyword{utilities}
