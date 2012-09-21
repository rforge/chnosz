



# Rpad utility functions.

"guiSource" <-
function (file, out.form = getOption("R.output.format"), local = FALSE,
          echo = verbose, print.eval = TRUE,
          verbose = getOption("verbose"), prompt.echo = getOption("prompt"),
          max.deparse.length = 150, chdir = FALSE)
{
    eval.with.vis <- function(expr, envir = parent.frame(), enclos = if (is.list(envir) ||
        is.pairlist(envir))
        parent.frame()) .Internal(eval.with.vis(expr, envir,
        enclos))
    envir <- if (local)
        parent.frame()
    else .GlobalEnv
    if (!missing(echo)) {
        if (!is.logical(echo))
            stop("echo must be logical")
        if (!echo && verbose) {
            warning(paste("verbose is TRUE, echo not; ... coercing",
                sQuote("echo <- TRUE")))
            echo <- TRUE
        }
    }
    if (verbose) {
        cat(sQuote("envir"), "chosen:")
        print(envir)
    }
    Ne <- length(exprs <- parse(n = -1, file = file))
    if (verbose)
        cat("--> parsed", Ne, "expressions; now eval(.)ing them:\n")
    if (Ne == 0)
        return(invisible())
    if (chdir && (path <- dirname(file)) != ".") {
        owd <- getwd()
        on.exit(setwd(owd))
        setwd(path)
    }
    if (echo) {
        sd <- "\""
        nos <- "[^\"]*"
        oddsd <- paste("^", nos, sd, "(", nos, sd, nos, sd, ")*",
            nos, "$", sep = "")
    }
    for (i in 1:Ne) {
        if (verbose)
            cat("\n>>>> eval(expression_nr.", i, ")\n\t\t =================\n")
        ei <- exprs[i]
        if (echo) {
            dep <- substr(paste(deparse(ei), collapse = "\n"),
                12, 1e+06)
            nd <- nchar(dep) - 1
            do.trunc <- nd > max.deparse.length
            dep <- substr(dep, 1, if (do.trunc)
                max.deparse.length
            else nd)
            cat("\n", prompt.echo, dep, if (do.trunc)
                paste(if (length(grep(sd, dep)) && length(grep(oddsd,
                  dep)))
                  " ...\" ..."
                else " ....", "[TRUNCATED] "), "\n", sep = "")
        }
        yy <- eval.with.vis(ei, envir)
        i.symbol <- mode(ei[[1]]) == "name"
        if (!i.symbol) {
            curr.fun <- ei[[1]][[1]]
            if (verbose) {
                cat("curr.fun:")
                str(curr.fun)
            }
        }
        if (verbose >= 2) {
            cat(".... mode(ei[[1]])=", mode(ei[[1]]), "; paste(curr.fun)=")
            str(paste(curr.fun))
        }
        if ( yy$visible ) {
          printoutput = capture.output(print(yy$value)) # always print, even if not shown for side effects
          if (out.form == "html" && exists("HTML"))
              HTML(yy$value)
          else if (out.form != "none")
              cat(paste(printoutput,collapse="\n"),"\n")
        }
        if (verbose)
          cat(" .. after ", sQuote(deparse(ei)), "\n", sep = "")
    }
    invisible(yy)
}


"RpadURL" <- function(filename = "") {
  # returns the URL for the given filename
  #   "./filename" for the local version
  #   "/Rpad/server/dd????????/filename" for the server version
  # use this to output HTML links for the user
  paste(get("RpadDir", envir = .RpadEnv), "/", filename, sep="")
}

"RpadBaseURL" <- function(filename = "") {
  # returns the base URL
  #   "filename" for the local version
  #   "/Rpad/filename" for the server version
  # use this to read in data files or save data files somewhere permanent
  if ( RpadIsLocal() )
    filename
  else
    paste("../../", filename, sep="")
}

"RpadBaseFile" <- function(filename = "") {
  # returns the file name relative to the base R directory
  #   "filename" for the local version
  #   "../../filename" for the server version
  # use this to read in data files or save data files somewhere permanent
  if ( RpadIsLocal() )
    paste("./", filename, sep="")
  else
    paste("../../", filename, sep="")
}

"RpadIsLocal" <- function() 
  get("RpadLocal", envir = .RpadEnv)    

