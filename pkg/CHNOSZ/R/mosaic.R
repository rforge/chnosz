# CHNOSZ/mosaic.R
# calculate affinities with changing basis species
# 20141220 jmd

# function to calculate affinities with mosaic of basis species
mosaic <- function(bases, blend=FALSE, ...) {
  # the arguments for affinity()
  myargs <- list(...)
  # are the swapped basis species on the plot?
  # (the first one should be present in the starting basis set)
  iswap <- match(bases[1], names(myargs))
  # the log activity of the starting basis species
  logact.swap <- basis()$logact[match(bases[1], row.names(basis()))]
  # a list where we'll keep the affinity calculations
  affs <- list()
  for(i in seq_along(bases)) {
    # set up argument list: name of swapped-in basis species
    if(!is.na(iswap)) names(myargs)[iswap] <- bases[i]
    # calculate affinities
    affs[[i]] <- do.call(affinity, myargs)
    # change the basis species; restore the original at the end of the loop
    if(i < length(bases)) {
      swap.basis(bases[i], bases[i+1]) 
      basis(bases[i+1], logact.swap)
    } else {
      swap.basis(bases[i], bases[1])
      basis(bases[1], logact.swap)
      names <- row.names(basis())
    }
  }
  # calculate affinities of formation of basis species
  ispecies <- species()$ispecies
  species.logact <- species()$logact
  species(delete=TRUE)
  species(bases)
  A.basis <- do.call(affinity, myargs)
  # restore original species with original activities
  species(delete=TRUE)
  species(ispecies, species.logact)
  # affinities calculated using the first basis species
  A.species <- affs[[1]]
  if(blend) {
    # calculate affinities using relative abundances of basis species
    e <- equilibrate(A.basis)
    for(j in seq_along(affs)) {
      for(i in seq_along(A.species$values)) {
        # start with zero affinity
        if(j==1) A.species$values[[i]][] <- 0
        # add affinity scaled by relative abundance of this basis species
        A.species$values[[i]] <- A.species$values[[i]] + affs[[j]]$values[[i]] * 10^e$loga.equil[[j]]
      }
    }
  } else {
    # use affinities from the single predominant basis species
    d <- diagram(A.basis, plot.it=FALSE)
    # merge affinities using the second, third, ... basis species
    for(j in tail(seq_along(affs), -1)) {
      is.predominant <- d$predominant==j
      for(i in seq_along(A.species$values)) {
        A.species$values[[i]][is.predominant] <- affs[[j]]$values[[i]][is.predominant]
      }
    }
  }
  # return the affinities for the species and basis species
  return(list(A.species=A.species, A.basis=A.basis))
}
