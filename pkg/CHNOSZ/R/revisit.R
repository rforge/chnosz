# CHNOSZ/revisit.R
# 20090415 functions related to diversity calculations
# 20100929 merged draw.diversity and revisit

where.extreme <- function(z, target, do.sat=FALSE) {
  if(missing(target)) stop("no target specified")
  # are we interested in a maximum or minimum?
  if(tolower(target) %in% c("sd", "sd.log", "cv", "cv.log", "rmsd", "cvrmsd", "dgxf"))
    myext <- "minimum" else myext <- "maximum"
  # do we care about the sign of the index?
  if(tolower(target) %in% c("sd", "sd.log", "cv", "cv.log", "rmsd", "cvrmsd")) 
    doabs <- TRUE else doabs <- FALSE
  # takes a matrix, returns the x,y coordinates of the extremum
  if(doabs) z <- abs(z)
  if(myext=="minimum") iext <- which.min(z)
  else iext <- which.max(z)
  ret.val <- iext
  # for a matrix, it gets more complicated esp.
  # if there are multiple instances of the extremum
  # (happens with saturated richnesses)
  if(do.sat & length(dim(z))==2) {
    iext <- which(z==z[iext])
    x.out <- y.out <- numeric()
    xres <- ncol(z)
    yres <- nrow(z)
    for(i in 1:length(iext)) {
      # column (x coord)
      x <- ceiling(iext[i]/xres)
      # and row (y coord)
      y <- iext[i] - floor(iext[i]/yres)*yres
      if(y==0) y <- yres  # there's a more eloquent way...
      x.out <- c(x.out,x)
      y.out <- c(y.out,y)
    }
    ret.val <- list(ix=y.out,iy=x.out)
  }
  return(ret.val)
}

extremes <- function(z, target) {
  if(missing(target)) stop("no target specified")
  # are we interested in a maximum or minimum?
  if(tolower(target) %in% c("sd", "sd.log", "cv", "cv.log", "rmsd", "cvrmsd", "dgxf")) 
    myext <- "minimum" else myext <- "maximum"
  # do we care about the sign of the index?
  if(tolower(target) %in% c("sd", "sd.log", "cv", "cv.log", "rmsd", "cvrmsd")) 
    doabs <- TRUE else doabs <- FALSE
  # takes a matrix, returns the y as f(x) and x as f(y)
  # trajectories of the extreme
  if(doabs) z <- abs(z)
  y <- x <- numeric()
  xres <- ncol(z)
  yres <- nrow(z)
  if(myext=="minimum") {
    for(i in 1:xres) y <- c(y, which.min(z[i,]))
    for(i in 1:yres) x <- c(x, which.min(z[,i]))
  } else {
    for(i in 1:xres) y <- c(y, which.max(z[i,]))
    for(i in 1:yres) x <- c(x, which.max(z[,i]))
  }
  return(list(x=x, y=y))
}

revisit <- function(d, target="cv", loga.ref=NULL,
  plot.it=NULL, col=par("fg"), yline=2, ylim=NULL, ispecies=NULL, add=FALSE,
  cex=par("cex"), lwd=par("lwd"), mar=NULL, side=1:4, xlim=NULL, labcex=0.6,
  pch=1, legend="", legend.x=NULL, lpch=NULL, main=NULL, lograt.ref=NULL, plot.ext=TRUE, DGxf.swap12=FALSE) {
  # calculate and plot diversity indices of relative abundances
  # 20090316 jmd
  # d can be the output from diagram (enables plotting)
  # or simply a list of logarithms of activity 
  # (each list entry must have the same dimensions)
  # test if the entries have the same dimensions
  ud <- unique(lapply(1:length(d),function(x) dim(d[[x]])))
  if(length(ud)==1) {
    # d is list of logarithms of activity
    if(missing(plot.it)) plot.it <- FALSE
    if(plot.it) stop("can't make a plot if argument 'd' is not the output from diagram()")
    logact <- d
  } else {
     # d is the output from diagram()
    if(!"logact" %in% names(d)) {
      stop(paste("the list provided in 'd' is not a usable result from diagram()",
        "(for two variables, diagram(..., mam=FALSE) is required)"))
    }
    if(missing(plot.it)) plot.it <- TRUE
    logact <- d$logact
  }
  # check that all needed arguments are present
  target.lower <- tolower(target)
  if(target.lower %in% c("richness", "cvrmsd", "spearman", "pearson", "logact", "dgxf") & is.null(loga.ref))
    stop(paste("for '", target, "' target, loga.ref must be supplied", sep=""))
  # take a subset (or all) of the species
  if(is.null(ispecies)) ispecies <- 1:length(logact)
  logact <- logact[ispecies]
  # number of species
  ns <- length(logact)
  # the dimensions 
  mydim <- dim(logact[[1]])
  nd <- length(mydim)
  if(nd==1) if(mydim==1) nd <- 0
  msgout(paste("revisit: calculating", target, "in", nd, "dimensions\n"))

  ## on to diversity calculations
  # given a list of logarithms of activities of species
  # (as vectors or matrices or higher dimensional arrays) 
  # calculate a diversity index of the same dimensions
  # these targets only depend on the logarithms of activities from diagram() (logact):
  # "shannon" shannon entropy
  # "sd" standard deviation
  # "cv" coefficient of variation
  # "sd.log", "cv.log" SD/CV for the logarithms of activity
  # "qqr" correlation coefficient on q-q plot (i.e., normality test)
  # these targets also depend on reference logarithms of activities (loga.ref):
  # "richness" species richness
  # "cvrmsd" coefficient of variation of rmsd
  # "spearman" spearman correlation coefficient
  # "pearson" pearson correlation coefficient
  # "logact" maximize the activity of a species
  # "DGxf" minimize the Gibbs energy of transformation
  
  # vectorize the entries in the logact list
  for(i in 1:ns) {
    logact[[i]] <- as.vector(logact[[i]])
    # convert infinite values to NA
    logact[[i]][is.infinite(logact[[i]])] <- NA
  }
  # make a place for the results
  H <- logact[[1]]
  H[] <- 0
  # ratio-specific calculations
  if(!is.null(lograt.ref)) {
    if(!target.lower %in% c("rmsd", "cvrmsd", "spearman", "pearson"))
      stop(paste("target",target,"not available when comparing activity ratios"))
    else(msgout("revisit: calculating activity ratios\n"))
    # instead of logact we use calculated log activity ratio
    logact <- lograt(loga.ref, logact)
    loga.ref <- lograt.ref
  }
  # target-specific calculations

  if(target.lower %in% c("sd", "cv")) {
    # build a matrix; rows are species
    myad <- t(list2array(logact))
    # remove logarithms
    myad <- 10^myad
    # get the standard deviations
    H <- palply(1:ncol(myad), function(i) sd(myad[,i]))
    H <- as.numeric(H)
    # for coefficient of variation, divide by the mean
    if(target.lower=="cv") H <- H / colMeans(myad)

  } else if(target.lower %in% c("sd.log", "cv.log")) {
    # build a matrix; rows are species
    myad <- t(list2array(logact))
    # get the standard deviations
    H <- palply(1:ncol(myad), function(i) sd(myad[,i]))
    H <- as.numeric(H)
    # for coefficient of variation, divide by the mean
    if(target.lower=="cv.log") H <- H / colMeans(myad)

  } else if(target.lower=="richness") {
    # given a list of logarithms of activities of species
    # (as matrices or vectors) calculate the richness 
    # make a place for the results
    # loop over species
    for(j in 1:length(logact)) {
      isthere <- logact[[j]] > loga.ref
      H[isthere] <- H[isthere] + 1
    }

  } else if(target.lower=="shannon") {
    # for this calculation we need to loop twice;
    # first to convert logact to act and get acttotal
    act <- logact
    for(i in 1:ns) {
      # exponentiate the logarithmic values
      act[[i]] <- 10^logact[[i]]
      if(i==1) acttotal <- act[[i]] else acttotal <- acttotal + act[[i]]
    }
    # now do the calculation
    for(i in 1:ns) {
      dH <- -act[[i]]/acttotal*log(act[[i]]/acttotal)
      if(!any(is.na(dH))) H <- H + dH
      else(warning(paste("revisit: skipping species", i, "which gives NA")))
    }

  } else if(target.lower=="qqr") {
    # normality test using correlation coefficient of a q-q plot 20100901
    actarr <- list2array(logact)
    qqrfun <- function(i, actarr) {
      y <- actarr[i,]
      # this is to catch errors from qqr (qqnorm)
      out <- try(qqr(y), silent=TRUE)
      if(class(out)=="try-error") out <- NA
      return(out)
    }
    H <- as.numeric(palply(1:length(H), qqrfun, actarr))

  } else if(target.lower=="spearman") {
    # spearman rank correlation coefficient 20100912
    actarr <- list2array(logact)
    spearfun <- function(i, actarr) {
      y <- actarr[i,]
      out <- spearman(y, loga.ref)
      return(out)
    }
    H <- as.numeric(palply(1:length(H), spearfun, actarr))

  } else if(target.lower=="pearson") {
    # pearson correlation coefficient
    actarr <- list2array(logact)
    pearfun <- function(i, actarr) {
      y <- actarr[i, ]
      out <- cor(y, loga.ref)
      return(out)
    }
    H <- as.numeric(palply(1:length(H), pearfun, actarr))

  } else if(target.lower=="rmsd") {
    # root mean squared deviation
    actarr <- list2array(logact)
    rmsdfun <- function(i, actarr) {
      y <- actarr[i, ]
      out <- rmsd(loga.ref, y)
      return(out)
    }
    H <- as.numeric(palply(1:length(H), rmsdfun, actarr))

  } else if(target.lower=="cvrmsd") {
    # coefficient of variation of the root mean squared deviation
    actarr <- list2array(logact)
    cvrmsdfun <- function(i, actarr) {
      y <- actarr[i, ]
      out <- cvrmsd(loga.ref, y)
      return(out)
    }
    H <- as.numeric(palply(1:length(H), cvrmsdfun, actarr))

  } else if(target.lower=="logact") {
    # where the activity of a species is maximal
    H <- logact[[loga.ref]]

  } else if(target.lower=="dgxf") {
    # Gibbs energy of transformation to the observed assemblage 
    actarr <- list2array(logact)
    # select species, vectorize, then put the Astar values into an array
    Astar <- d$Astar[ispecies]
    for(i in 1:ns) Astar[[i]] <- as.vector(Astar[[i]])
    Astararr <- list2array(Astar)
    Gfun <- function(i, actarr, Astararr) {
      loga.equil <- actarr[i, ]
      Astar <- Astararr[i, ]
      # direction of transformation:
      # swap12=FALSE: loga.equil(Astar) --> loga.ref
      # swap12=TRUE:  loga.ref(Astar) --> loga.equil
      if(DGxf.swap12) out <- -DGxf(loga.ref, loga.equil, Astar)
      else out <- DGxf(loga.equil, loga.ref, Astar)
      return(out)
    }
    H <- as.numeric(palply(1:length(H), Gfun, actarr, Astararr))
  
  } else stop(paste("specified target '", target, "' not available", sep=""))
  # replace dims
  dim(H) <- mydim
  ## now on to plotting + assembling return values
  # get information about the x-axis
  if(plot.it & nd > 0 & nd < 3) {
    xname <- d$xname
    yname <- d$yname
    xres <- d$xlim[3]
    xrange <- d$xlim[1:2]
    # special operations for pH
    if(xname=="H+") {
      xname <- "pH"
      xrange <- -xrange 
    }
    # the x-values
    xs <- seq(xrange[1],xrange[2],length.out=xres)
  }
  # make plots and return values
  if(nd==0) {
    # a 0-D plot
    if(plot.it) {
      plotted <- FALSE
      if(target.lower=="qqr") {
        actarr <- list2array(logact)
        qqnorm(actarr,col=col,pch=pch,main="")
        qqline(actarr)
        plotted <- TRUE
      } else if(target.lower %in% c("rmsd","cvrmsd","spearman","pearson")) {
        # plot the points
        ylab <- "loga.calc"
        xlab <- "loga.ref"
        if(!is.null(lograt.ref)) {
          ylab <- "lograt.calc"
          xlab <- "lograt.ref"
        }
        plot(loga.ref,actarr,xlab=xlab,ylab=ylab,pch=pch,col=col)
        # add a 1:1 line
        lines(range(loga.ref),range(loga.ref),col="grey")
        # add a lowess line
        ls <- loess.smooth(loga.ref,actarr)
        lines(ls$x,ls$y,col="red")
        plotted <- TRUE
      }
      if(plotted) {
        # add a title
        if(missing(main)) main <- paste(target,"=",round(H,3)) 
        title(main=main)
        if(!is.null(legend.x)) {
          if(is.null(lpch)) lpch <- unique(pch)
          legend(legend.x,legend=legend,pch=lpch)
        }
      }
    }
    ret.val <- list(H=H)

  } else if(nd==1) {
    # locate the extremum
    ix <- where.extreme(H,target)
    extval <- H[ix]
    # a 1-D plot
    if(plot.it) {
      if(is.null(ylim)) {
        if(target.lower=="richness") {
            ylim <- 0
            if(max(H) > ylim) ylim <- max(H) + 1
            ylim <- c(0,ylim)
        }
        else ylim <- extendrange(H,f=0.075)
      }
      if(is.null(xlim)) xlim <- xrange
      # format the target name if it's DGxf
      if(target.lower=="dgxf") ylab <- expr.property("DGxf/2.303RT")
      else ylab <- target
      if(!add) thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label(xname),ylab=ylab,yline=yline,
        cex=cex,lwd=lwd,mar=mar,side=side)
      # plot the values
      lines(xs,as.vector(H),col=col)
      x <- xs[ix]
      if(plot.ext) abline(v=x,lty=2)
      ret.val <- list(H=H,ix=ix,x=x,extval=extval)
    } else ret.val <- list(H=H,ix=ix,extval=extval)

  } else if(nd==2) {
    # a 2-D plot
    iext <- where.extreme(H,target.lower,do.sat=TRUE)
    # what is the extreme value
    ix <- iext$ix
    iy <- iext$iy
    extval <- H[ix,iy]
    ret.val <- list(H=H,ix=ix,iy=iy,extval=extval) 
    if(plot.it) {
      yres <- d$ylim[3]
      yrange <- d$ylim[1:2]
      #if(yname=="T") yrange <- outvert(yrange,"K")
      if(yname=="H+") {
        yname <- "pH"
        yrange <- -yrange 
      }
      ys <- seq(yrange[1],yrange[2],length.out=yres)
      # start the plot
      if(is.null(xlim)) xlim <- xrange
      if(is.null(ylim)) ylim <- yrange
      if(!add) thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label(xname),ylab=axis.label(yname),
        yline=yline,side=side,cex=cex,mar=mar)
      contour(xs,ys,H,add=TRUE,labcex=labcex)
      # plot the location(s) of the extremum
      points(xs[iext$ix],ys[iext$iy],pch=8,cex=2)
      # show trajectories of the extrema
      iexts <- extremes(H,target.lower)
      # take out large jumps
      yext <- ys[iexts$y]
      yext.1 <- c(yext[2:length(yext)],yext[length(yext)])
      yext.2 <- c(yext[1],yext[1:length(yext)-1])
      yext[abs(yext.1-yext)/abs(diff(range(ys))) > 0.1] <- NA
      yext[abs(yext.2-yext)/abs(diff(range(ys))) > 0.1] <- NA
      lines(xs,yext,lty=3,col="blue")
      xext <- xs[iexts$x]
      xext.1 <- c(xext[2:length(xext)],xext[length(xext)])
      xext.2 <- c(xext[1],xext[1:length(xext)-1])
      xext[abs(xext.1-xext)/abs(diff(range(xs))) > 0.1] <- NA
      xext[abs(xext.2-xext)/abs(diff(range(xs))) > 0.1] <- NA
      lines(xext,ys,lty=3,col="seagreen")
      ret.val <- list(H=H,ix=ix,iy=iy,
        x=xs[ix],y=ys[iy],extval=extval)
    }

  } else {
    # we don't make plots for more than two dimensions
    # just return the values
    ret.val <- list(H=H)
  }
  # return the results
  return(invisible(ret.val))
}
