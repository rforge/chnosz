# CHNOSZ/species.R
# define species of interest 

# to retrieve the coefficients of reactions to form the species from the basis species
species.basis <- function(species) {
  # current basis matrix
  bmat <- basis.matrix()
  tbmat <- t(bmat)
  # what are the elements?
  belem <- rownames(tbmat)
  # get the species makeup into a matrix
  mkp <- as.matrix(sapply(makeup(species, count.zero=TRUE), c))
  # the positions of the species elements in the basis elements
  ielem <- match(rownames(mkp), belem)
  # the elements of the species must be contained by the basis species
  if(any(is.na(ielem))) stop(paste("element(s) not in the basis:", 
    paste(rownames(mkp)[is.na(ielem)], collapse=" ")))
  # the positions of the basis elements in the species elements
  jelem <- match(belem, rownames(mkp))
  # keep track of which ones are NA's; 
  # index them as one here but turn them to zero later
  ina <- is.na(jelem)
  jelem[ina] <- 1
  # now put the species matrix into the same order as the basis
  mkp <- mkp[jelem, , drop=FALSE]
  # fill zeros for any basis element not in the species
  mkp[ina, ] <- 0
  # solve for the basis coefficients and transpose
  nbasis <- t(apply(mkp, 2, function(x) solve(tbmat, x)))
  # very small numbers are probably a floating point artifact
  # can cause problems in situations where zeros are needed
  # (manifests as issue in longex("phosphate"), where which.balance()
  #  identifies H2O as conserved component)
  out <- zapsmall(nbasis)
  # add names of species and basis species
  colnames(out) <- colnames(tbmat)
  # add names of species only if it was a character argument
  if(all(is.character(species))) rownames(out) <- species
  return(out)
} 

species <- function(species=NULL, state=NULL, delete=FALSE, index.return=FALSE) {
# 20080925 changed default to quiet=TRUE 
# 20101003 changed default to quiet=FALSE
# 20120128 remove 'quiet' argument (messages can be hidden with suppressMessages())
# 20120523 return thermo$species instead of rownumbers therein, and remove message showing thermo$species
  missingstate <- missing(state)
  state <- state.args(state)

  # we can't deal with NA species
  if(identical(species, NA)) {
    cn <- caller.name()
    if(length(cn) > 0) ctext <- paste("(calling function was ", cn, ")", sep="") else ctext <- ""
    stop(paste("'species' is NA", ctext))
  }

  # delete the entire species definition or only selected species
  if(delete) {
    # remember the old species definition
    oldspecies <- thermo$species
    # delete the entire definition if requested
    if(is.null(species)) {
      thermo$species <<- NULL
      return(oldspecies)
    }
    # from here we're trying to delete already defined species
    if(is.null(oldspecies)) stop("nonexistent species definition")
    # match species to delete by name and species number
    isp <- rep(NA, length(species))
    ispname <- match(species, thermo$species$name)
    ispnum <- match(species, 1:nrow(oldspecies))
    isp[!is.na(ispname)] <- ispname[!is.na(ispname)]
    isp[!is.na(ispnum)] <- ispnum[!is.na(ispnum)]
    # filter out non-matching species
    ina <- is.na(isp)
    if(any(ina)) warning(paste("species:",
      paste(species[ina], collapse=" "), "not present, so can not be deleted"))
    isp <- isp[!ina]
    # go on to delete this/these species
    if(length(isp) > 0) {
      thermo$species <<- thermo$species[-isp,]
      if(nrow(thermo$species)==0) thermo$species <<- NULL
      else rownames(thermo$species) <<- 1:nrow(thermo$species)
    }
    return(thermo$species)
  } 

  # if no species or states are given, just return the species list
  if(is.null(species) & is.null(state)) return(thermo$species)
  # if no species are given use all of them if available
  if(is.null(species) & !is.null(thermo$species)) species <- 1:nrow(thermo$species)

  if(length(species) > length(state) & !is.null(state)) state <- rep(state,length.out=length(species)) else 
  if(length(state) > length(species) & !is.null(species)) species <- rep(species,length.out=length(state))

  # if they don't look like states (aq,gas,cr) or activities (numeric), 
  # use them as a suffix for species name (e.g., a protein-organism)
  if( length(which(state %in% unique(as.character(thermo$obigt$state)))) < 
    length(state) & !can.be.numeric(state[1]) & !can.be.numeric(species[1]) ) {
      for(i in 1:length(state)) species[i] <- paste(species[i],'_',state[i],sep='')
      state <- rep(thermo$opt$state,length.out=length(state))
  }

  # append/change species entries
  if(is.null(thermo$basis)) stop('basis species are not defined')
  if(is.character(species[1])) {
    # character first argument, species in thermo$obigt
    # but only give states if they are numeric
    is <- NULL
    if(!can.be.numeric(state[[1]])) is <- state
    ispecies <- suppressMessages(info(species, is))
    # check if we got all the species
    ina <- is.na(ispecies)
    # info() returns a list if any of the species had multiple approximate matches
    # we don't accept any of those
    if(is.list(ispecies)) ina <- ina | sapply(ispecies,length) > 1
    if(any(ina)) stop(paste("species not available:",paste(species[ina],collapse=" ")))
    if(length(ispecies)==0) return(species())
    was.character <- TRUE
  } else if(is.numeric(species[1])) {
    ispecies <- species
    ispecies <- ispecies[!is.na(ispecies)]
    if(length(ispecies)==0) return(species())
    was.character <- FALSE
  }
  jspecies <- ispecies
  myspecies <- NULL
  if(!is.null(thermo$species)) 
    if(TRUE %in% (ispecies %in% thermo$species$ispecies)) {
      myspecies <- match(ispecies,thermo$species$ispecies)
      myspecies <- thermo$species$ispecies[myspecies[!is.na(myspecies)]]
      ispecies <- ispecies[!ispecies%in%myspecies]
    }
  # only add to an existing dataframe if the indices can't
  # all possibly refer to the rows of the dataframe
  doit <- TRUE
  ## might (not) work well when you want to add e.g. H2O after some others
  if(!is.null(thermo$species)) if(all(jspecies %in% 1:nrow(thermo$species))) 
    if(!was.character) doit <- FALSE
  if(length(ispecies) > 0 & !(is.numeric(ispecies[1]) & is.numeric(state[1]) ) & doit) {
    # the coefficients in reactions to form the species from basis species
    f <- (species.basis(ispecies))
    # the default states and activities
    state <- as.character(thermo$obigt$state[ispecies])
    logact <- numeric()
    for(i in 1:length(state)) {
      if(length(agrep('missing',rownames(f)[i]))>0) la <- NA
      else { if(state[i]=='aq') la <- -3 else la <- 0 }
      logact <- c(logact,la)
    }
    # yes, the species names too
    name <- as.character(thermo$obigt$name[ispecies])
    # add the ispecies values
    t <- data.frame(f,ispecies=ispecies,logact=logact,state=state,name=name,stringsAsFactors=FALSE)
    # nasty for R, but "H2PO4-" looks better than "H2PO4."
    colnames(t)[1:nrow(thermo$basis)] <- rownames(thermo$basis)
    if(is.null(thermo$species)) thermo$species <<- t else 
      thermo$species <<- rbind(thermo$species,t)
    rownames(thermo$species) <<- seq(1:nrow(thermo$species))
  }
  #  update activities or states
  if(!is.null(state)) {
    state <- rep(state,length.out=length(jspecies))
    # if is looks like species aren't set yet, try to do so
    if(is.null(thermo$species)) { 
      species(jspecies)
    } else {
      mj <- jspecies[!jspecies %in% thermo$species$ispecies]
      if(!can.be.numeric(species[1])) species(mj)
    }
    # we bet that the number of rows is smaller
    # than the indices of whatever species we have
    if(can.be.numeric(species[1]) & max(jspecies) <= nrow(thermo$species))
      jspecies <- thermo$species$ispecies[jspecies]
    mj <- match(jspecies,thermo$species$ispecies)
    if(can.be.numeric(state[1])) {
      if(NA %in% mj[1]) warning(paste('can\'t update activity of species',
        c2s(which(is.na(mj))),' requested'),call.=FALSE)
      thermo$species$logact[mj] <<- state
    } else {
      mj <- match(jspecies,thermo$species$ispecies)
      state <- rep(state,length.out=length(mj))
      name <- thermo$species$name[mj]
      # try to check that the states actually exist
      for(k in 1:length(mj)) {
        doit <- TRUE
        if(NA %in% mj[k]) doit <- FALSE
        myform <- thermo$obigt$formula[thermo$species$ispecies[mj[k]]]
        #iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] | thermo$obigt$formula==myform)
        # 20080925 don't match formula -- two proteins might have the
        # same formula (e.g. YLR367W and YJL190C)
        #iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]])
        # 20091112 do match formula if it's not a protein -- be able to 
        # change "carbon dioxide(g)" to "CO2(aq)"
        if(length(grep("_",thermo$species$name[mj[k]])) > 0)  
          iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]])
        else {
          iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] & thermo$obigt$state==state[k])
          if(length(iobigt)==0)
            iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] | thermo$obigt$formula==myform)
        }
        if(!state[k] %in% thermo$obigt$state[iobigt]) 
          doit <- FALSE
        if(!doit) warning(paste('can\'t update state of species ',
          mj[k],' to ',state[k],'.\n',sep=''),call.=FALSE)
        else {
          ii <- match(state[k],thermo$obigt$state[iobigt])
          thermo$species$state[mj[k]] <<- state[k]
          thermo$species$name[mj[k]] <<- thermo$obigt$name[iobigt[ii]]
          thermo$species$ispecies[mj[k]] <<- as.numeric(rownames(thermo$obigt)[iobigt[ii]])
        }
      }
    }
  } else {
    # this message turns out to be kinda distracting. what should be here?
    #if(!is.null(myspecies)) 
      #cat(paste('species: keeping ',c2s(thermo$obigt$name[myspecies],sep=', '),'.\n',sep='')) 
  }
  # return the new species definition or index(es) of affected species
  if(index.return) return(match(jspecies, thermo$species$ispecies))
  else return(thermo$species)
}

