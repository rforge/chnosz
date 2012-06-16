# CHNOSZ/equil.R
# functions to calculation logarithm of activity
# of species in (metastable) equilibrium

equil.boltzmann <- function(Astar, n.balance, loga.balance) {
  # 20090217 new "abundance" function
  # return equilibrium logarithms of activity of species
  # works using Boltzmann distribution
  # A/At = e^(Astar/n.balance) / sum(e^(Astar/n.balance))
  # A/At = e^(Astar/n.balance) / sum(e^(Astar/n.balance))
  # where A is activity of the ith residue and
  # At is total activity of residues
  # advantages over equil.reaction
  # 2) loops over species only - much faster
  # 3) no root finding - those games might fail at times
  # disadvantage:
  # 1) only works for per-residue reactions
  # 2) can create NaN logacts if the Astars are huge/small

  # initialize output object
  A <- Astar
  # remember the dimensions of elements of Astar (could be vector,matrix)
  Astardim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # first loop: make vectors
  A <- palply(1:length(A),function(i) as.vector(A[[i]]))
  # second loop: get the exponentiated Astars (numerators)
  # need to convert /2.303RT to /RT
  #A[[i]] <- exp(log(10)*Astar[[i]]/n.balance[i])/n.balance[i]
  A <- palply(1:length(A),function(i) exp(log(10)*Astar[[i]]/n.balance[i]))
  # third loop: accumulate the denominator
  # initialize variable to hold the sum
  At <- A[[1]]; At[] <- 0
  for(i in 1:length(A)) At <- At + A[[i]]*n.balance[i]
  # fourth loop: calculate log abundances and replace the dimensions
  A <- palply(1:length(A),function(i) loga.balance + log10(A[[i]]/At))
  # fifth loop: replace dimensions
  for(i in 1:length(A)) dim(A[[i]]) <- Astardim
  # add names and we're done!
  names(A) <- Anames
  return(A)
}

equil.reaction <- function(Astar, n.balance, loga.balance) {
  # to turn the affinities/RT (A) of formation reactions into 
  # logactivities of species (logact(things)) at metastable equilibrium
  # 20090217 extracted from diagram and renamed to abundance.old
  # 20080217 idea / 20120128 cleaned-up strategy
  # for any reaction stuff = thing,
  #   A = logK - logQ 
  #     = logK - logact(thing) + logact(stuff)
  # given Astar = A + logact(thing),
  # given Abar = A / n.balance,
  #   logact(thing) = Astar - Abar * n.balance  [2]
  # where n.balance is the number of the balanced quanitity
  # (conserved component) in each species
  # equilibrium values of logact(thing) satifsy:
  # 1) Abar is equal for all species
  # 2) log10( sum of (10^logact(thing) * n.balance) ) = loga.balance  [1]
  #
  # because of the logarithms, we can't solve the equations directly
  # instead, use uniroot() to compute Abar satisfying [1]
  # remember the dimensions and names
  Adim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # make a matrix out of the list of Astar
  Astar <- list2array(lapply(Astar, c))
  # that produces the same result (other than colnames) and is much faster than
  #Astar <- sapply(Astar, c)  
  # also, latter has NULL nrow for length(Astar[[x]])==1
  # some function definitions
  # to calculate log of activity of balanced quantity from logact(thing) of all species [1]
  logafun <- function(logact) log10(sum(10^logact * n.balance))
  # to calculate logact(thing) from Abar for the ith condition [2]
  logactfun <- function(Abar, i) Astar[i,] - Abar * n.balance
  # to calculate difference between logafun and loga.balance for the ith condition
  logadiff <- function(Abar, i) loga.balance - logafun(logactfun(Abar, i))
  # to calculate a range of Abar that gives negative and positive values of logadiff for the ith condition
  Abarrange <- function(i) {
    # starting guess of Abar (min/max) from range of Astar / n.balance
    Abar.range <- range(Astar[i, ] / n.balance)
    # the range of logadiff
    logadiff.min <- logadiff(Abar.range[1], i)
    logadiff.max <- logadiff(Abar.range[2], i)
    # we're out of luck if they're both infinite
    if(is.infinite(logadiff.min) & is.infinite(logadiff.max))
      stop("FIXME: there are no initial guesses for Abar that give
        finite values of the differences in logarithm of activity
        of the conserved component")
    # if one of them is infinite we might have a chance
    if(is.infinite(logadiff.min)) {
      # decrease the Abar range by increasing the minimum
      # the 0.99 is necessary at least for the residue=FALSE plot in protactiv.Rnw; 0.9 isn't enough!
      Abar.range[1] <- Abar.range[1] + 0.99 * diff(Abar.range)
      logadiff.min <- logadiff(Abar.range[1], i)
      if(is.infinite(logadiff.min)) stop("FIXME: the second initial guess for Abar.min failed")
    }
    if(is.infinite(logadiff.max)) {
      # decrease the Abar range by decreasing the maximum
      Abar.range[2] <- Abar.range[2] - 0.99 * diff(Abar.range)
      logadiff.max <- logadiff(Abar.range[2], i)
      if(is.infinite(logadiff.max)) stop("FIXME: the second initial guess for Abar.max failed")
    }
    # the change of logadiff with Abar
    # it's a weighted mean of the n.balance
    dlogadiff.dAbar <- (logadiff.max - logadiff.min) / diff(Abar.range)
    # change Abar to center logadiff (min/max) on zero
    logadiff.mean <- mean(c(logadiff.min, logadiff.max))
    Abar.range <- Abar.range - logadiff.mean / dlogadiff.dAbar
    # one iteration is enough for the examples in the package
    # but there might be a case where the range of logadiff doesn't cross zero
    # (e.g. for the carboxylic acid example in revisit.Rd)
    logadiff.min <- logadiff(Abar.range[1], i)
    logadiff.max <- logadiff(Abar.range[2], i)
    if(logadiff.min > 0 | logadiff.max < 0) {
      # do again what we did above
      dlogadiff.dAbar <- (logadiff.max - logadiff.min) / diff(Abar.range)
      logadiff.mean <- mean(c(logadiff.min, logadiff.max))
      Abar.range <- Abar.range - logadiff.mean / dlogadiff.dAbar
      logadiff.min <- logadiff(Abar.range[1], i)
      logadiff.max <- logadiff(Abar.range[2], i)
      if(logadiff.min > 0 | logadiff.max < 0) {
        stop("FIXME: make this function (Abarrange() in equil.reaction()) 
          iterate again to find a range of Abar such that the differences in 
          logarithm of activity of the conserved component cross zero")
      }
    }
    return(Abar.range)
  }
  # to calculate an equilibrium Abar for the ith condition
  Abarfun <- function(i) {
    # get limits of Abar where logadiff brackets zero
    Abar.range <- Abarrange(i)
    # now for the real thing: uniroot!
    Abar <- uniroot(logadiff, interval=Abar.range, i=i)$root
    return(Abar)
  }
  # calculate the logact(thing) for each condition
  logact <- palply(1:nrow(Astar), function(i) {
    # get the equilibrium Abar for each condition
    Abar <- Abarfun(i)
    return(logactfun(Abar, i))
  }) 
  # restore the dimensions and names
  if(length(Adim)==1) logact <- list2array(logact)
  else logact <- sapply(logact, c)
  logact <- lapply(1:nrow(logact), function(i) {
    thisla <- list(logact[i,])[[1]]
    dim(thisla) <- Adim
    return(thisla)
  })
  names(logact) <- Anames
  # all done!
  return(logact)
}
