# CHNOSZ/zzz.R
# this has the .onAttach function for package startup message and data initialization

.onAttach <- function(libname,pkgname) {
  # version figuring adapted from package mgcv
  pkghelp <- library(help=CHNOSZ)$info[[1]]
  # things are different for older versions of R
  if(length(pkghelp)==1) pkghelp <- library(help=CHNOSZ)$info[[2]]
  version <- pkghelp[pmatch("Version:", pkghelp)]
  um <- strsplit(version, " ")[[1]]
  version <- um[nchar(um)>0][2]
  date <- pkghelp[pmatch("Date:", pkghelp)]
  um <- strsplit(date, " ")[[1]]
  date <- um[nchar(um)>0][2]
  # identify the program and version
  packageStartupMessage(paste("CHNOSZ version ", version, " (", date, ")", sep=""))
  # message if 'parallel' package is available but not loaded
  if(!"parallel" %in% sessionInfo()$basePkgs)
    if(length(find.package("parallel", quiet=TRUE)) > 0) 
      packageStartupMessage("Suggested package 'parallel' available but not loaded")
  # load the 'thermo' data object
  data(thermo)

  ## load data files in user's directory
  #if(file.exists('protein.csv')) {
  #  tr <- try(rbind(read.csv('protein.csv',as.is=TRUE),thermo$protein),silent=TRUE)
  #  if(identical(class(tr),'try-error')) cat("thermo: protein.csv in current directory is not compatible with thermo$protein data table.\n")
  #  else add.protein("protein.csv")
  #}
  #if(file.exists('obigt.csv')) {
  #  tr <- try(rbind(read.csv('obigt.csv',as.is=TRUE),thermo$obigt),silent=TRUE)
  #  if(identical(class(tr),'try-error')) cat("thermo: obigt.csv in current directory is not compatible with thermo$obigt data table.\n")
  #  else add.obigt("obigt.csv")
  #}
}
