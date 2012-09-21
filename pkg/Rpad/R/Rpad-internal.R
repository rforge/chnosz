.RpadEnv <- new.env()

".onLoad" <-
function(lib, pkg) {
    require("graphics")
    require("utils")
    isR2HMTLAvailable <- length(.find.package("R2HTML", quiet = TRUE)) != 0
    if (isR2HMTLAvailable) { 
      options(R2HTML.sortableDF = TRUE)
      options(R2HTML.format.digits = 3)
      options(R2HTML.format.nsmall = 0)
      options(R2HTML.format.big.mark = "")
      options(R2HTML.format.big.interval = 3)
      options(R2HTML.format.decimal.mark = Sys.localeconv()[["decimal_point"]])
      .HTML.file <<- ""
    }
    # The following uses the environment variable DOCUMENT_ROOT with apache to find
	# the directory of the R process. Change may be required for another server.
	if (Sys.getenv("DOCUMENT_ROOT") != "") { # works for Apache 1.3 linux & win
		RpadDir <- gsub(Sys.getenv("DOCUMENT_ROOT"), "",
                        getwd(),  ignore.case = TRUE) # strip off the document root
	} else if (Sys.getenv("SCRIPT_NAME") != "") { # for Apache 2.0
        RpadDir <- paste(gsub("R_process.pl", "", Sys.getenv("SCRIPT_NAME"), ignore.case = TRUE),
                         gsub(".*/", "", getwd()),
                         sep="")
	} else if (Sys.getenv("PATH_INFO") != "") { # for microsoft IIS
        RpadDir <- paste(gsub("R_process.pl", "", Sys.getenv("PATH_INFO"), ignore.case = TRUE),
                         gsub(".*/", "", getwd()),
                         sep="")
	} else {
		.rootdir = ifelse(.Platform$OS.type == "windows", "C:/www", "/var/www")
		RpadDir <- gsub(.rootdir, "", getwd(),  ignore.case = TRUE)
	}

	options(R.output.format = "text") # do we need or want this anymore? If we do, why don't we put it in .RpadEnv?

    
	assign(".RpadGraphOptions",
           list(type = "pngalpha", extension = "png",
                res = 120, width = 0, height = 0, deviceUsesPixels = TRUE, pointsize = 10,
                sublines = 0, toplines = .6, ratio = 4/3, leftlines = 0, lwd = 0.6),
           envir = .RpadEnv)
    assign("RpadLocal", FALSE, envir = .RpadEnv)
    assign("RpadDir", RpadDir, envir = .RpadEnv)
    assign("Rpad.plot.counter",  0, envir = .RpadEnv)
}

".onUnload" <-
function(libpath) {
 	if (interactive()) stopRpadServer()
}

".packageName" <- "Rpad"
