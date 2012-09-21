\name{RpadHTML}
\alias{ROutputFormat}
\alias{HTMLon}     
\alias{HTMLoff}    
\alias{H}     
\alias{HtmlTree}     
\alias{HtmlTree.default}     
\alias{Html}     
\alias{HfromHTML}
\alias{HTMLh1}     
\alias{HTMLh2}     
\alias{HTMLh3}     
\alias{HTMLh4}     
\alias{HTMLh5}     
\alias{HTMLargs}      
\alias{HTMLtag}      
\alias{HTMLetag}      
\alias{HTMLradio}      
\alias{HTMLcheckbox}   
\alias{HTMLselect}     
\alias{HTMLinput}      
\alias{HTMLlink}       
\alias{HTMLimg}        
\alias{HTMLembed}
\alias{BR}     
\alias{print.condition}
                   
\title{ Rpad utilities }
\description{      
  Rpad utilities to generate HTML output.
}
\usage{
  ROutputFormat(Format)
  HTMLon()
  HTMLoff()
  HtmlTree(tagName, ...)
  H(tagName, ...)
  \method{HtmlTree}{default}(tagName, ..., standaloneTag = FALSE, collapseContents = TRUE)
  Html(x, ...)
  BR
  HTMLh1(text)
  HTMLh2(text)
  HTMLh3(text)
  HTMLh4(text)
  HTMLh5(text)
  HTMLargs(x)
  HTMLtag(tagName, ...)
  HTMLetag(tagName)
  HTMLradio(variableName, commonName = "radio", text = "", ...)
  HTMLcheckbox(name, text = "", ...)
  HTMLselect(name, text, default = 1, size = 1, id = name, contenteditablewrapper = TRUE, optionvalue = text, ...)
  HTMLinput(name, value = "", rpadType = "Rvariable", contenteditablewrapper = TRUE, ...)
  HTMLlink(url, text, ...)
  HTMLimg(filename = RpadPlotName(), ...)
  HTMLembed(filename, width = 600, height = 600, ...)
  HfromHTML(x)
  print.condition(x, ...)
}

\arguments{
  \item{Format}{ a string specifying the desired output format, like "html", "text", or "none" }
  \item{text}{ specifies the text to be displayed in the HTML output. }
  \item{tagName}{ is a string specifying the HTML nodeName ("H1" or
    "DIV" for example). }
  \item{standaloneTag}{ is a logical that specifies whether the HTML tag
    is a standalone tag, meaning it has no ending tag (<br/> for example). }
  \item{collapseContents}{ is a logical that specifies whether vector
    tag contents are collapsed (merged).}
  \item{variableName, commonName}{ specify attributes of the radio
    element. \code{commonName} specifies the common radio elements that
    are grouped together (the name attribute of the radio element)--this gets translated into the R variable with
    the same name. \code{variableName} is the string result assigned to
    \code{commonName} when this radio item is selected (the value attribute of the radio element).}
  \item{name}{ is the name attribute in HTML that is translated into the
    R variable with the same name.}
  \item{value}{ is the initial value of the input box. }
  \item{url}{ is the URL for the <A> element.}
  \item{contenteditablewrapper}{ is a logical that specifies whether a
    wrapper is placed around the output to disable editing. Internet
    Explorer needs this to make links and input elements active.}
  \item{default, size, id, optionvalue}{ specify parameters of a select
    box. \code{default} specifies which of the options is selected initially.}
  \item{filename, width, height}{ are parameters of an IMG or EMBED object.}
  \item{rpadType}{ is the "rpadType" attribute of an INPUT
    element. Normal values are either "Rvariable" or "Rstring" but other
    values of "rpadType" should also be possible. }
  \item{x}{ == To define == }
  \item{...}{ arguments are passed on as HTML arguments for the given tag.}
}
\details{
\code{HTMLon} and \code{HTMLoff} specify how Rpad
interprets results (either HTML or text). Rpad output sections are normally rendered as plain text with
a fixed-width font, so text script outputs are formatted
properly. The output blocks can also contain HTML codes for
displaying images or for formatting blocks. To change to HTML
formatting, use \code{HTMLon()} (which sends '<htmlon/>' to the output).
To turn HTML mode off, use \code{HTMLoff()} (which sends
'<htmloff/>'). \code{HTMLon} and \code{HTMLoff} only apply to the
existing Rpad input section; they don't carry over into the next.

\code{HtmlTree} and it's shortcut \code{H} generates HTML (or other
XML-like syntax) for an arbitrary tag with the specified name and
arguments.  Unnamed parameters to \code{H} are content, and named
parameters are tag attributes.  By using nesting, it's a good way to
generate HTML (or any XML) with properly formed tags. If \code{tagName}
is NULL, then the contents are created without a surrounding tag. For
vectors, each vector element becomes surrounded by tags, so
\code{H("div", c(1,2,3))} results in
\code{<div>1</div><div>2</div><div>3</div>}.

\code{Html} converts an R object to class \code{HtmlTree}. Methods are
defined for matrices and data frames. For data frames, all ... arguments
are passed to \code{format} before conversion to an HTML
representation. Unnamed parameters to \code{H} that have a method
defined for \code{Html} are automatically converted. For example,
\code{H("div", data.frame(a=1:2, b=1:2))} is the HTML representation of
the data frame wrapped by a div tag.

For conversion of a wider selection of R objects to HTML, use
\code{HTML} from the R2HTML package. \code{HTML} differs from
\code{Html} in that \code{HTML} uses \code{cat} to directly send output
whereas \code{Html} returns a character string. \code{HfromHTML}
captures the output of \code{HTML} and returns it as a character string
of class HtmlTree. This is useful for nesting \code{HTML} inside of
\code{H} trees.


\code{ROutputMode} changes how R behaves with automatic
printing. Possible values are: "text" (the default), "html", or
"none". "text" is like at the command line: values returned in the
script are automatically printed (without an explicit \code{print}
statement) in standard text format. With "html", values returned are
automatically printed, but HTML output is generated by using the
\code{HTML} method from R2HTML (if available) instead of the \code{print}
method. \code{ROutputMode} applies to all subsequent Rpad input section,
including a rollover back to the beginning when a page is run several
times.

}

\value{
  All of the HTML code generation routines return a character string of
  class HtmlTree. With automatic printing, the string is sent to the
  output (with \code{cat}). This has the effect that the HTML is
  displayed in Rpad. Note that if a HTML generation function is used
  inside of a loop or other scenario where the function results are not
  automatically printed, then you need to enclose the function with
  \code{cat} or \code{print}.

  \code{HTMLargs} returns a string with the arguments as "a='arg1'
  b='arg2'", and so on. Note the use of single quotes. This will affect
  how quoted strings should be passed as elements.
}

\author{ Tom Short, EPRI, (\email{tshort@epri.com}) }

\seealso{\code{R2HTML} for other useful HTML routines that can be used in Rpad. }

\examples{
  a <- data.frame(x = 1:10, y = (1:10)^3)

  # generate some GUI elements
  # - normally done in an Rpad input section with rpadRun="init"
  HTMLon() # switch to html mode
  data(state)
  HTMLselect("favoriteAmericanState", state.name) # generate the select box

  H("div", "Hello world") # a simple div: "<div>Hello world</div>"
  H("div", class="helloStyle", "Hello world") # a div with a class
    # --> "<div class='helloStyle'>Hello world</div>"

  # you can nest them:
  H("div", class = "myClass", dojoType = "tree",
    "This is some text in the div ",  
    H("em", "emphasized"), "plus some more",
    "and more",
    H("div", class = "anotherClass",
      "text in the div",
      c(1,5,8)))

  HTMLoff() # switch to text mode
  summary(a)

  # fancy HTML output with R2HTML
  if (require(R2HTML)) {
    HTMLon()
    .HTML.file = "" # not needed in normal use
                    # kludge to pass R CMD check
    HTML(summary(a))
    HTMLoff()
  }

}

\keyword{ math }
