# CHNOSZ/more.aa.R
# get amino acid compositions of proteins from 
# model organisms really exciting!
# (Eco.csv or Sce.csv)

more.aa <- function(protein=NULL, organism) {
  # return the composition of one or more proteins from
  # a "model organism", E. coli (Eco) or S. cerevisiae (Sce)
  # extracted from get.protein 20120519
  datapath <- paste("extdata/protein/", organism, ".csv.xz", sep="")
  datafile <- system.file(datapath, package="CHNOSZ")
  if(datafile=="") stop(paste("missing", datapath))
  mydata <- read.csv(datafile)
  # if protein is not supplied, just give some information about the datafile
  if(is.null(protein)) {
    msgout("more.aa: ", datapath, " has data for ", nrow(mydata), " proteins\n")
    return(invisible())
  }
  if(organism=="Sce") {
    # which columns to search for matches
    searchcols <- c("OLN", "OLN")
    # which columns have the amino acids in the order of thermo$protein 
    iaa <- c(1,5,4,7,14,8,9,10,12,11,13,3,15,6,2,16,17,20,18,19) + 2
  } else if(organism=="Eco") {
    # which columns to search for matches
    searchcols <- c("protein", "abbrv")
    # which columns have the amino acids in the order of thermo$protein 
    iaa <- 1:20 + 5
  }
  # find the matches
  icols <- match(searchcols, colnames(mydata))
  imatch <- match(protein, mydata[, icols[1]])
  imatch2 <- match(protein, mydata[, icols[2]])
  # use not-NA matches for "abbrv" in Eco.csv
  imatch[!is.na(imatch2)] <- imatch2[!is.na(imatch2)]
  # report and remember the unsuccessful matches
  if(all(is.na(imatch))) stop("no proteins found!")
  inotmatch <- which(is.na(imatch)) 
  if(length(inotmatch) > 0) {
    if(length(inotmatch)==1) verb <- " was" else verb <- " were"
    msgout("more.aa: ", paste(protein[inotmatch], collapse=" "), verb, " not matched\n")
  }
  aa <- data.frame(mydata[imatch, iaa])
  # add the identifying columns
  organism <- rep(organism, length(protein))
  abbrv <- rep(NA, length(protein))
  ref <- rep(NA, length(protein))
  chains <- rep(1, length(protein))
  chains[inotmatch] <- NA
  precols <- data.frame(protein, organism, ref, abbrv, chains, stringsAsFactors=FALSE)
  colnames(aa) <- aminoacids(3)
  aa <- cbind(precols, aa)
  # done!
  return(aa)
}


