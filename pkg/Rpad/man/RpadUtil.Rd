\name{RpadUtil}
\alias{RpadURL}
\alias{RpadBaseURL}     
\alias{RpadBaseFile}    
\alias{RpadIsLocal}     
                   
\title{ Rpad utilities }
\description{      
  Rpad utilities to generate filenames or URL's
}
\usage{
  RpadURL(filename = "")
  RpadBaseURL(filename = "")
  RpadBaseFile(filename = "")
  RpadIsLocal()
}

\arguments{
  \item{filename}{ the name of a file. }
}
\value{
  \code{RpadURL} returns the URL for the given filename: "./filename" for
  the local version of Rpad and "/Rpad/server/dd????????/filename" for
  the server version. Use this to output HTML links for the user.

  \code{RpadBaseURL} returns the base URL: "filename" for the local
  version and "/Rpad/filename" for the client-server version. Use this
  to point the user to data files or other links on the server that is
  somewhere permanent. (The current R working directory is not
  permanent in the client-server version.)

  \code{RpadBaseFile} returns the file name relative to the base R
  directory: "filename" for the local version and "../../filename" for
  the client-server version. Use this in R to read in data files or save
  data files somewhere permanent.

}

\author{ Tom Short, EPRI Solutions, Inc., (\email{tshort@eprisolutions.com}) }

\seealso{\code{Rpad}, \code{RpadHTML}}

\examples{
\dontrun{
  # make some data
  x <- 1:10
  y2 <- x^3
  save(x, y2, file = RpadBaseFile("testdata.RData"))
  # output a link to the user:
  HTMLon()
  cat("<a href='", RpadBaseURL("testdata.RData"), sep="")
  cat("'>Click</a> to download the test data.")
  unlink(RpadBaseFile("testdata.RData"))
}
}

\keyword{ math }
