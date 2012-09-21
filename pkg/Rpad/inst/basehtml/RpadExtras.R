
#
#
# NOTE -- many of these may be out of date or broken. 2006/12/15
#




"HTML.list" <- function (x, file = .HTML.file, first = TRUE,  
                         append = TRUE, collapsed = TRUE, align="left", ...)  
# This replaces the HTML.list in R2HTML.  
# It creates an expandable tree that's collapsed by default.  
# This uses Detlef Groth's tree.htc behavior.  
# The tree.htc is attached in the style heading of  
# this html file.  
# Adding an expandable tree is as easy as wrapping standard HTML  
# lists with <span class='tree'>.
# See also the MoreGuiExamples.Rpad example.
{ 
   cat("\n", file = file, append = append, ...) 
   if (first) { # IE needs contenteditable off 
     cat("<div contenteditable='false'><span class='tree'><ul>",  
         file = file, append = TRUE, sep = "\n") 
   } 
   for (i in seq(along=x)) { 
     cat('<li style="list-style-image: url(img/', ifelse(collapsed,"plus","minus"),'.gif);', 
         ' color: rgb(0, 0, 51); cursor: pointer;">', 
         file = file, append = TRUE, sep = "") 
     cat(names(x)[i])  
     cat("<ul ",ifelse(collapsed,"style='display:none'",""),">",  
         file = file, append = TRUE, sep = "\n") 
     HTML(x[[i]], file = file, first = FALSE, collapsed=collapsed, align=align, ...) 
     cat("</ul></li>", file = file, append = TRUE, sep = "\n") 
   } 
   if (first) { 
     cat("</ul></span></div>", file = file, append = TRUE, sep = "\n") 
   }  
} 
"HTML.data.frame.min" <- 
function (x, ...)
# minimal (but smaller output) HTML output of a data frame
{
    x <- format(x)
    x <- as.matrix(x)
    x[is.na(x)] <- " "
    x[is.nan(x)] <- " "
    x[grep("^ *(NA|NaN) *$", x)] <- " "
    header <- paste("<table", HTMLargs(...),">")
    rowheader <- paste(paste(c('<tr><th>',rep.int('<th>', ncol(x) - 1)),
                             c(as.matrix(colnames(x))),
                             c(rep.int('</th>', ncol(x) - 1),'</th></tr>'),
                             sep = ""),
                       collapse="")
    body <- paste(paste(c('<tr><td>', rep.int('<td>', ncol(x) - 1)),
                        c(t(x)),
                        c(rep.int('</td>', ncol(x) - 1),'</td></tr>'),
                        sep = ""),
                  collapse="")
    trailer <- "</table>" 
    str <- paste(header, rowheader, body, trailer, sep="")
    class(str) <- "HTMLchunk"
    str
}
"HTML.individual.row.of.df" <- 
function (x, ...)
# make an HTML table for each row of 'x'
{
    x <- format(x)
    x <- as.matrix(x)
    x[is.na(x)] <- " "
    x[is.nan(x)] <- " "
    x[grep("^ *(NA|NaN) *$", x)] <- " "
    header <- paste('"<table ', HTMLargs(...),">")
    rowheader <- paste(paste(c('<tr><th>',rep.int('<th>', ncol(x) - 1)),
                             c(as.matrix(colnames(x))),
                             c(rep.int('</th>', ncol(x) - 1),'</th></tr>'),
                             sep = ""),
                       collapse="")
    trailer <- '</table>",'
    paste(paste(c(paste(header,rowheader,'<tr><td>'), rep.int('<td>', ncol(x) - 1)),
                c(t(x)),
                c(rep.int('</td>', ncol(x) - 1),paste('</td></tr>',trailer)),
                sep = ""),
          collapse="")
}
"runjavascript" <-
function(filename) 
  # loads a javascript file and runs it
  # must have HTML enabled
  #cat('&amp;nbsp;&lt;SCRIPT TYPE="text/javascript" DEFER src="data.js"&gt;&lt;SCRIPT&gt;')<br>
  # this replaces the above line, which doesn't work under IE.
  # see https://lists.latech.edu/pipermail/javascript/2004-January/006879.html
  cat('&nbsp;<SCRIPT TYPE="text/javascript" DEFER>q="\u0022";(e=document.createElement("script")).src="',
      filename,
      '";e.type="text/javascript"; document.getElementsByTagName("head")[0].appendChild(e);</SCRIPT>',
      sep="")

"HTMLjs" <-
function(js) {
  # Runs a javascript snippet. HTML must be enabled.
  str <- paste(
    HTMLtag("script", type="text/javascript"),
    js,
    HTMLetag("script"))
  class(str) <- "HTMLchunk"
  str
}

"HTMLsetInnerHtml" <-
function(id, html) {
  # Sets the innerHTML of the element "id" to "html".
  # Runs a bit of javascript to do it. HTML must be enabled.
  HTMLjs(paste("top._Rpad_editor._doc.getElementById('", id, "').innerHTML = '", html, "'", sep=""))
}

"HTMLbutton" <-
function(label = "Calculate", js = "top.Rpad_calculate()") 
  # Other useful js parameters:
  #   js = "top.Rpad_calculate_next(this)"  # calculate the next Rpad block
  #   js = "top.R_run_commands('put commands here')"
  HTMLtag("input", onclick = paste("javascript:", js, sep=""),
          value = label, type = "button")

# This is a file upload block of HTML.
# It is used with the fileupload.pl on the server side.

# The process is fairly convoluted. Xmlhttp can't be used to upload files,
# so we have to use a hidden iframe to "post" the uploaded file.
# When the file has finished uploading, perl sends back a snippet
# of javascript that turns around and executes an R command (normally to
# load or do some manipulation on the uploaded file).
# For an even fancier file upload widget, see
# http://blog.joshuaeichorn.com/archives/2005/05/01/ajax-file-upload-progress/
#
# The functions RresponseToUpload & HTMLupload go together.
#
# This is the default R command after a function is uploaded.
# Redefine it for other functionality.
"RresponseToUpload" <-
function(name) {
  if (grep(".RData$", name))
    load(name, envir = globalenv())
}
# This function generates the HTML for the form.
# It can be used dynamically, or you can paste the
# generated HTML into a form.
"HTMLupload" <- function(dir = NULL) {
  # "dir" is relative to the Rpad base directory.
  # if "dir" is null, R's current working directory is used.
  activedir = rev(strsplit(RpadURL(), '/')[[1]])[1]
  if (is.null(dir)) dir = activedir
  else dir = paste(activedir, "/../../", dir, sep="")
  str <- 
  paste(HTMLtag('span', contenteditable='false'),
        HTMLtag('form', method="post", enctype="multipart/form-data", target="target_iframe",
                action="server/fileupload.pl"),
        HTMLtag('input', type="file", name="uploadfilename"),
        HTMLtag('input', name="activedir", value=dir, type="hidden"),
        HTMLtag('input', type="submit", value="Submit file"),
        HTMLtag('iframe', id="target_iframe", name="target_iframe",
                src="", style="width:1px;height:1px;border:0"),
        HTMLetag('iframe'),
        HTMLetag('form'),
        HTMLetag('span'))
  class(str) <- "HTMLchunk"
  str
}
#Example usage:
# HTMLupload(dir="server")

