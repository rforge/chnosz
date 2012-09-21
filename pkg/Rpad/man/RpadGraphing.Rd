\name{RpadGraphing}
\alias{graphoptions}   
\alias{newgraph}   
\alias{showgraph}
\alias{RpadPlotName}
                   
\title{ Rpad graphing utilities }

\description{
  Utilities to make graphing in Rpad R scripts easier.
}

\usage{
  graphoptions(..., reset = FALSE, override.check = TRUE)
  newgraph(name = "", extension = graphoptions()$extension,
     type = graphoptions()$type, res = graphoptions()$res,
     width = graphoptions()$width, height = graphoptions()$height,
     deviceUsesPixels = graphoptions()$deviceUsesPixels,
     pointsize = graphoptions()$pointsize, sublines = graphoptions()$sublines,
     toplines = graphoptions()$toplines, ratio = graphoptions()$ratio,
     leftlines = graphoptions()$leftlines, lwd = graphoptions()$lwd, ...)
  showgraph(name = RpadPlotName(), link = FALSE, ...)
  RpadPlotName()
}

\arguments{
  \item{reset, override.check}{ logical arguments passed to \code{check.options}.  }
  \item{name}{ is the name of the graph with the extension left OFF.}
  \item{extension}{ is the file extension of the device (typically '.png').}
  \item{type}{ is the graphics device, either a ghostscript device as a
    character string, "Rpng" for R's builtin png device, or any R
    graphics device function. For a ghostscript device, common
    possibilities are "png256", "pngalpha" (the default), or "pdf".}
  \item{res}{ is the resolution of the bitmap in dots per inch.}
  \item{width, height}{ are the dimensions of the graph in inches.}
  \item{deviceUsesPixels}{ is a logical specifying whether the R
    graphics driver uses pixels. This only applies if \code{type} is
    used to specify an R graphics driver. It defaults to TRUE.}
  \item{ratio}{ specifies the ratio of the graph if either the width or
    height is not specified.}
  \item{pointsize}{ is the font point size passed to the postscript
    device.}
  \item{sublines, toplines, leftlines}{ specify the dimension of
    the graph along with the outside border. It defaults to fairly tight
    outside dimensions.}
  \item{lwd}{ is the line width set with \code{par}.}
  \item{link}{ is a logical specifying whether to show a link to allow
    the user to download an EPS file of the graph (not available when
    using R's png driver.}
  \item{...}{ in \code{fspdf} and \code{newgraph}, arguments are passed
    to \code{postscript}. In \code{graphoptions}, the arguments are
    assigned as defaults for \code{newgraph}.}
}
\details{
      
The \code{graphoptions}, \code{newgraph}, and \code{showgraph} set of functions
allows quick setup and display of web-friendly graphics. The user can
normally just use any of the plot commands followed by
\code{showgraph}. \code{newgraph} sets up the graphics
device, and it's called when the Rpad package
starts. \code{showgraph} generates the HTML to show the graph and runs
\code{newgraph} to advance
to the next graphics file. The user only runs \code{newgraph} to change
parameters from their defaults. Graphics files are by default named
Rpad\_plot1, Rpad\_plot2, and so on. Named graphs can also be used, but
there's more of a chance that if the user has caching set wrong (or the
server's caching is set wrong) that graphs won't update properly in the
user's browser. With the default sequential numbering of files, caching
problems are less likely. \code{graphoptions} is also available to change the
defaults for subsequent graphs.

Internally, \code{newgraph} uses the postscript device and ghostscript
to generate the bitmap for the browser. The pngalpha
(\code{type="pngalpha"}) device of ghostscript does anti-aliasing for
smoother-looking PNG output. Also, this approach has the side-effect of
creating an EPS file for each graph, so it's easy to add a link to allow
the user to download the EPS file. A second approach is to use R's
\code{png} device directly by specifying \code{type="Rpng"} (but on
linux, this requires an X server, and results are not antialiased). A
third approach allows arbitrary R graphic's devices to be
specified. Simon Urbanek's GDD or Cairo packages can be used to generate
png (for both packages, see http://www.rosuda.org/R/). This allows
high-quality png generation with anti-aliasing, and you don't need
ghostscript or X11. 
}

\value{
  \code{RpadPlotName()} returns the name of the currently active
  plot. None of the other routines return a value: all are run for their
  side effects. 
}

\author{ Tom Short, EPRI, (\email{tshort@epri.com}) }

\seealso{ See also \code{bitmap}, \code{png}, and \code{pdf}. }

\examples{
# make some graphs (a default graphics device is already available)
  x <- 1:10
  y <- x^2
  y2 <- x^3
  if (capabilities("png")) graphoptions(type="Rpng")
  newgraph()
  plot(x, y)  # does the plot
  plot(x, y2) # does the second plot
  HTMLon()    # sets Rpad to HTML output
  showgraph() # closes the device, outputs the HTML for the both
              # images, and creates the next device
  plot(x, y2)
  showgraph()

# graphs with named files:
  newgraph("graph_A")
  plot(x, y)
  showgraph("graph_A")
  newgraph("graph_B", width = 4, height = 6) # also adjust the width and height
  plot(x, y2)
  showgraph("graph_B")

# an example with Simon Urbanek's GDD device:
  if (require(GDD)) {
    graphoptions(type = GDD, width = 4, height = 3) # note the use of inches
    newgraph()
    plot(x, y)
    plot(x, y2)
    showgraph()
  }
}

\keyword{ math }
