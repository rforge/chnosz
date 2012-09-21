
# Rpad utility functions for running Rpad locally.
# Here we use a local Tcl httpd server to receive Rpad commands.

"processRpadCommands" <-
function() {
  require("tcltk")
  commands <- tclvalue(.Tcl("set user(R_commands)"))
  textcommands <- textConnection(commands)
  .dev.active <- dev.cur()
  if (exists("RpadPlotParams", envir = .RpadEnv))
    dev.set( get("RpadPlotParams", envir = .RpadEnv)$dev )

  results <- tryCatch({
    tc <- textConnection("textfromconnection",open="w")
    sink(file=tc)
    guiSource(textcommands)
    sink()
    close(tc)
    textfromconnection
  }, error=function(e) {
    sink()
    close(tc)
    cat('ERROR1: ')
    paste(paste(textfromconnection, "\n", collapse=""), '\n', e)},
                      finally=close(textcommands))
  dev.set(.dev.active)
  formattedresults <- paste(results,"\n",sep="",collapse="")
#  cat(formattedresults)
  escapeBrackets <- function(x) gsub("(\\{|\\})", "\\\\\\1", x)
  .Tcl(paste("set RpadTclResults {", escapeBrackets(formattedresults), "}", sep=""))
}


"Rpad" <-
function(file = "", defaultfile = "LocalDefault.Rpad", port = 8079) {
    startRpadServer(defaultfile, port)
    browseURL(paste("http://127.0.0.1:", port, "/", file, sep = ""))
}

"startRpadServer" <-
function(defaultfile = "LocalDefault.Rpad", port = 8079) {
    require("tcltk")
    # This is the main function that starts the server
	# This function implements a basic http server on 'port'
 	# The server is written in Tcl.
    # This way it is not blocking the R command-line!

    if (!require("tcltk")) stop("package tcltk required for the local Rpad http server")
    assign("RpadLocal", TRUE, envir = .RpadEnv)
    assign("RpadDir",   ".",  envir = .RpadEnv)
    assign("RpadPort",  port, envir = .RpadEnv)
    graphoptions(type = "Rpng")
    tclfile <- file.path(.find.package(package = "Rpad"), "tcl", "mini1.1.tcl")
    htmlroot <- file.path(.find.package(package = "Rpad"), "basehtml")
    tcl("source", tclfile)
    tcl("Httpd_Server", htmlroot, port, defaultfile)
    unlink(dir(pattern="Rpad_plot.*\\.png")) # delete the Rpad graphics files in the dir
    unlink(dir(pattern="Rpad_plot.*\\.eps"))
    if(interactive() && .Device == "null device") x11() # turn on the interactive plotting device so as not to confuse the command-line user if they later plot
    dev <- dev.cur() 
    newgraph()
    dev.set(dev) # switch back to the existing device to not confuse the user
    return(TRUE)
}

"stopRpadServer" <-
function() {
    require("tcltk")
    assign("RpadLocal",    FALSE, envir = .RpadEnv)
    unlink(dir(pattern="Rpad_plot.*\\.png")) # delete the Rpad graphics files in the dir
    unlink(dir(pattern="Rpad_plot.*\\.eps"))
    .Tcl("close $Httpd(listen)")
    .Tcl("unset Httpd")
}

"restartRpadServer" <-
function() {
  stopRpadServer()
  startRpadServer()
}


