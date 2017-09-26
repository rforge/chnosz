# CHNOSZ/yeast.aa.R
# get amino acid compositions of proteins from Saccharomyces cerevisiae

yeast.aa <- function(protein=NULL) {
  # return the composition of one or more proteins from S. cerevisiae (Sce)
  # extracted from get.protein 20120519
  datapath <- paste("extdata/protein/Sce.csv.xz", sep="")
  datafile <- system.file(datapath, package="CHNOSZ")
  if(datafile=="") stop(paste("missing", datapath))
  mydata <- read.csv(datafile, as.is=TRUE)
  # if protein is not supplied, just give some information about the datafile
  if(is.null(protein)) {
    message("yeast.aa: ", datapath, " has data for ", nrow(mydata), " proteins")
    return(invisible())
  }
  # which columns to search for matches
  searchcols <- c("ORF", "SGDID", "GENE")
  # which columns have the amino acids, in the order of thermo$protein 
  iaa <- match(toupper(aminoacids(3)), toupper(colnames(mydata)))
  # iterate over a list
  waslist <- TRUE
  out <- list()
  if(!is.list(protein)) {
    waslist <- FALSE
    protein <- list(protein)
  }
  for(i in 1:length(protein)) {
    # find the matches
    imatch <- rep(NA, length(protein[[i]]))
    for(cname in searchcols) {
      icol <- match(cname, colnames(mydata))
      if(is.na(icol)) next
      iimatch <- match(protein[[i]], mydata[, icol])
      imatch[!is.na(iimatch)] <- iimatch[!is.na(iimatch)]
    }
    # report and remember the unsuccessful matches
    if(all(is.na(imatch))) stop("no proteins found!")
    inotmatch <- which(is.na(imatch)) 
    if(length(inotmatch) > 0) {
      if(length(inotmatch)==1) verb <- " was" else verb <- " were"
      message("yeast.aa: ", paste(protein[[i]][inotmatch], collapse=" "), verb, " not matched")
    }
    aa <- data.frame(mydata[imatch, iaa])
    # add the identifying columns
    ref <- mydata$SGDID[imatch]
    abbrv <- mydata$GENE[imatch]
    chains <- rep(1, length(protein[[i]]))
    chains[inotmatch] <- NA
    org <- rep("Sce", length(protein[[i]]))
    precols <- data.frame(protein[[i]], organism=org, ref, abbrv, chains, stringsAsFactors=FALSE)
    colnames(precols)[1] <- "protein"
    colnames(aa) <- aminoacids(3)
    aa <- cbind(precols, aa)
    out <- c(out, list(aa))
  }
  # done!
  if(!waslist) return(out[[1]])
  else return(out)
}


