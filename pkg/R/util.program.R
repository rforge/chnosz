# CHNOSZ/util.program.R
# various programming-related functions

caller.name <- function(n=2) {
  # returns the name of the calling function n frames up
  # (n=2: the caller of the function that calls this one)
  # or character() if called interactively
  if(sys.nframe() < n) name <- character()
  else {
    sc <- sys.call(-n)[[1]]
    name <- try(as.character(sc),silent=TRUE)
    # also return character() if the value from sys.call is
    # the function itself (why does this sometimes happen,
    # e.g. when called from affinity()?)
    if(class(name)=="try-error") name <- character()
  }
  return(name)
}

palply <- function(X, FUN, ...) {
  # a wrapper function to run parLapply if length(X) > 100 
  # and package 'parallel' is available, otherwise run lapply
  if(length(X) > 100 & "parallel" %in% (.packages())) {
    # tell the user the number of iterations and the name of the calling function
    cn <- caller.name()
    ###if(length(cn)==0) ntext <- "" else ntext <- paste("(called by ", cn, ")", sep="")
    ###cat("palply: ", ntext, " ", length(X), " parLapply calculations ...", sep="")
    # the actual calculations - modified from ?parLapply
    ## Use option mc.cores to choose an appropriate cluster size.
    # or detectCores if that is NULL
    # don't load methods package
    cl <- makeCluster(getOption("mc.cores", detectCores()), methods=FALSE)
    out <- parLapply(cl, X, FUN, ...)
    stopCluster(cl)
    ###cat(" done!\n")
  } else out <- lapply(X,FUN,...)
  return(out)
}
