package.contents <- function (pkg, lib.loc = NULL)
{
    file <- system.file("CONTENTS", package = pkg, lib.loc = lib.loc)
#    if (file == "") {
#        warning(paste("Cannot find CONTENTS file of package",
#            pkg))
#        return(NA)
#    }
    read.dcf(file = file, fields = c("Entry", "Aliases", "Keywords", "Description"))
#    read.dcf(file = file)
}

makeIndexData <- function(lib.loc = rev(.libPaths()) ) {
  res <- NULL
  firstpg <- c("base","utils", "stats", "graphics", "methods") 
  for (i in firstpg) 
    res <- rbind(res,package.contents(i))
  for (lib in lib.loc) {
    pg <- .packages(all.available=TRUE, lib.loc=lib)
    pg <- pg[!(pg %in% firstpg)]
    for (i in pg)
      try(res <- rbind(res,package.contents(i)))
  }
  res
}

writeIndexData <- function(res) { # then must edit after
  t <- file("helpIndex.txt","w")
  cat(file=t,"var Entry =\n")
  dput(res[,"Entry"],file=t)
  cat(file=t,"var Keywords =\n")
  dput(res[,"Keywords"],file=t)
  cat(file=t,"var Aliases =\n")
  dput(res[,"Aliases"],file=t)
  cat(file=t,"var Description =\n")
  dput(res[,"Description"],file=t)
  close(t)
#  system('sed s/^c/Array/ helpIndex.txt > js/helpIndex.js')
  st = readLines("helpIndex.txt")
  st = gsub("^c\\(", "Array\\(", st)
  writeLines(text=st, "js/helpIndex.js")
}

id=makeIndexData()
writeIndexData(id)
