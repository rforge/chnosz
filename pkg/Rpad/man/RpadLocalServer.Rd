\name{RpadLocalServer}
\alias{Rpad}
\alias{startRpadServer}
\alias{stopRpadServer}

\title{ Rpad local server }

\description{      
  Functions to implement Rpad locally.

  Rpad is an interactive, web-based analysis program. Rpad pages
  are interactive workbook-type sheets based on R, an open-source
  implementation of the S language. Rpad is an analysis package, a
  web-page designer, and a gui designer all wrapped in one. Rpad
  makes it easy to develop powerful data analysis applications that
  can be shared with others.

  Rpad is available in two versions: a local
  version and an intranet/internet version. The local version works
  through the user's local installation of R with through the user's web
  browser. The intranet/internet
  version works in client-server fashion with the user accessing a remote
  server through a standard web browser.
}

\usage{
  Rpad(file = "", defaultfile = "LocalDefault.Rpad", port = 8079)
  startRpadServer(defaultfile = "LocalDefault.Rpad", port = 8079)
  stopRpadServer()
}

\arguments{
  \item{file}{ the file to load into the browser. }
  \item{defaultfile}{ the default filename for the minihttpd server to
    serve if the URL does not specify one.}
  \item{port}{ the TCP port of the server (8079 by default).}
}

\details{
  \code{Rpad()} starts the local Rpad server and launches the default
  browser with a default startup page. Use the \code{file} argument to
  specify a different starting page. The default startup page allows the
  user to select any of the Rpad html files in R's current working
  directory. You can also use \code{startRpadServer()} to start the server
  (without launching another browser window)
  and \code{stopRpadServer()} to stop the server.

  The Rpad local server implements a mini-httpd, a minimal web-page server. This mini server is implemented in Tcl/Tk, using the powerful 'socket' command.
  Since it runs in the separate tcltk event loop, it is not blocking R, and it
  runs in the background; the user can still enter commands at the R
  prompt. The user can use Rpad along side of other user interfaces,
  including Sciviews-R, Rgui, ESS, and/or Rcmdr. 

  The mini-httpd server first looks for files
relative to R's current working directory. If it can't find them
there, it looks for files relative to the "basehtml" directory in the
Rpad package directory. This way, the user can store Rpad html files
wherever he wants and not have to worry about carrying around the
javascript, CSS, and other html-related files.

  The original implementation of Rpad uses a classical web server like Apache with perl scripts.
  It is working as a client-server through Intra/Internet (look at http://www.Rpad.org)
  for a live example.

  A number of R utility functions are provided for
accessing directories and URL's that should keep compatibility between
the local version and the client/server version (\code{RpadURL}, \code{RpadBaseURL}, and
\code{RpadBaseFile}).
}

\author{Philippe Grosjean (\email{phgrosjean@sciviews.org}) and Tom Short, EPRI Solutions, Inc., (\email{tshort@eprisolutions.com})}

\note{
  For security reasons, the server can only run for a local client. However,
  it is very easy to eliminate this limitation by hacking the \code{startSocketServer()}
  function in the \pkg{svSocket} package (SciViews bundle).
}

\seealso{ See \code{link{RpadServer}} for information on the
  client/server version of Rpad. For utility functions for compatibility
between the local and client/server versions of Rpad, see \code{RpadURL}, \code{RpadBaseURL}, and
\code{RpadBaseFile}.}

\keyword{ IO }
