# CHNOSZ/water.R
# calculate thermodynamic and electrostatic properties of H2O
# 20061016 jmd

water.AW90 <- function(T=298.15,rho=1000,P=0.1) {
  # Equations for the dielectric constant of water
  # from Archer and Wang, 1990
  # T in K
  # rho in kg m-3
  # p in MPa

  # Table 2
  b <- c(-4.044525E-2, 103.6180   , 75.32165   ,
         -23.23778   ,-3.548184   ,-1246.311   ,
         263307.7    ,-6.928953E-1,-204.4473)
  alpha <- 18.1458392E-30 # m^3
  #alpha <- 14.7E-30
  mu <- 6.1375776E-30 # C m
  N.A <- 6.0221367E23 # mol-1
  k <- 1.380658E-23 # Boltzmann constant, J K-1
  M <- 0.0180153 # kg mol-1
  rho.0 <- 1000 # kg m-3
  # Equation 1
  epsilon.0 <- 8.8541878E-12 # permittivity of vacuum, C^2 J-1 m-1
  epsfun.lhs <- function(e) (e-1)*(2*e+1)/(9*e)
  epsfun.rhs <- function(T,V.m) N.A*(alpha+mufun()/(3*epsilon.0*k*T))/(3*V.m)
  epsfun <- function(e,T,V.m) epsfun.lhs(e) - epsfun.rhs(T,V.m)
  mufun <- function() gfun()*mu^2
  gfun <- function() rhofun()*rho/rho.0 + 1
  # Equation 3
  rhofun <- function() b[1]*P*T^-1 + b[2]*T^-0.5 + b[3]*(T-215)^-1 +
    b[4]*(T-215)^-0.5 + b[5]*(T-215)^-0.25 +
    exp(b[6]*T^-1 + b[7]*T^-2 + b[8]*P*T^-1 + b[9]*P*T^-2)
  epsilon <- function(T,rho) {
    tu <- try(uniroot(epsfun,c(1E-1,1E3),T=T,V.m=M/rho)$root,TRUE)
    if(!is.numeric(tu)) {
      warning('water.AW90: no root for density at ',T,' K and ',rho,' kg m-3.',call.=FALSE,immediate.=TRUE)
      tu <- NA
    }
    return(tu)
  }
  # get things the right length
  our.T <- T; our.rho <- rho; our.P <- P
  t <- numeric()
  for(i in 1:length(our.T)) {
    T <- our.T[i]
    rho <- our.rho[i]
    P <- our.P[i]
    t <- c(t,epsilon(T,rho))
  }
  return(t)
}

water <- function(property = NULL,T = thermo$opt$Tr, P = 'Psat') {
  # calculate the properties of liquid H2O as a function of T and P
  # T in Kelvin, P in bar
  if(is.null(property)) stop('property was NULL')
  # this tells us to do the calculations using code taken from SUPCRT
  do.supcrt <- length(agrep(tolower(thermo$opt$water),'supcrt9',max.distance=0.3)) > 0
  eargs <- eos.args('water',property=property,T=T,P=P)
  property <- eargs$prop; Property <- eargs$Prop
  # working out the arguments
  tpargs <- TP.args(T=T,P=P)
  P <- tpargs$P
  T <- tpargs$T
  # Psat stuff
  psat <- function() {
    p <- numeric()
    if(do.supcrt) {
      p <- water.SUPCRT92('',T=T,rep(0,length(T)),isat=1)
      p[p==0] <- NaN
      return(p)
    } else {
      for(i in 1:length(T)) {
        if(T[i] < 373.124) p <- c(p,0.1)
        else p <- c(p,WP02.auxiliary('P.sigma',T[i]))
      }
      return(convert(p,'bar'))
    }
  }
  # a quick return if property = 'Psat', for use by the TP.args() function
  if(length(property)==1 & property[1]=='psat') return(data.frame(Psat=psat()))
  ### maybe we are using the SUPCRT calculations ###
  if(do.supcrt) {
    names.SUPCRT <- c('Speed','alpha','beta','alpha','beta','diel','ZBorn','YBorn','QBorn','XBorn')
    names.CHNOSZ <- c('w','alpha','beta','E','kT','epsilon','Z','Y','Q','X')
    Property.new <- character()
    # convert names to SUPCRT
    for(i in 1:length(Property)) if(Property[i] %in% names.CHNOSZ) 
      Property.new[i] <- names.SUPCRT[match(Property[i],names.CHNOSZ)]
      else Property.new[i] <- Property[i]
    # deal with compressibility and expansivity 20091203
    iE <- which(Property=="E")
    ikT <- which(Property=="kT")
    iV <- numeric()
    if("kT" %in% Property | "E" %in% Property) iV <- length(Property.new <- c(Property.new,"V"))
    # get the value of the property
    w.out <- water.SUPCRT92(Property.new,T=T,P=P)
    # finish dealing with compressibility and expansivity
    if("E" %in% Property) w.out[,iE] <- w.out$V*w.out$alpha
    if("kT" %in% Property) w.out[,ikT] <- w.out$V*w.out$beta
    if(length(iV) > 0) w.out <- w.out[,-iV,drop=FALSE]
    colnames(w.out) <- Property
    return(w.out)
  } else {
    # here we get properties using IAPWS-95 
    w.out <- water.IAPWS95(property, T, P)
    colnames(w.out) <- Property
    return(w.out)
  }
}


water.SUPCRT92 <- function(property,T=298.15,P=1,isat=0) {
  ### interface to H2O92D.f : FORTRAN subroutine taken from 
  ### SUPCRT92 for calculating the thermodynamic and 
  ### electrostatic properties of H2O. 
  ## we restrict the calculations to liquid water
  ## except for getting Psat (vapor-liquid saturation 
  ## pressure as a function of T>100 C). 20071213 jmd
  # H2O92 doesn't output Born functions N or U
  if('n' %in% tolower(property) | 'uborn' %in% tolower(property))
    stop('I can\'t tell you the Born functions N or U (used in calculating compressibilities and expansibilities of aqueous species).')
  # pressure setting
  if(is.null(P)) P <- rep(0,length(T))
  # values to use here gleaned from H2O92D.f and SUP92D.f
  # it, id, ip, ih, itripl, isat, iopt, useLVS, epseqn, icrit
  if(isat) iopt <- 1 else iopt <- 2  # for Psat(T) (1) or T-P (2)
  specs <- c(2,2,2,5,1,isat,iopt,1,4,0)
  states <- rep(0,4)
  # match up properties with the output
  props <- c('a','g','s','u','h','cv','cp','Speed','alpha',
    'beta','diel','visc','tcond','surten','tdiff','Prndtl',
    'visck','albe','ZBorn','YBorn','QBorn','daldT','XBorn')
  iprop <- seq(1,45,length.out=23)
  # now to the actual calculations
  Tc <- convert(T,'C')
  # initialize the output matrix
  w.out <- matrix(NA,nrow=length(T),ncol=23,byrow=TRUE) 
  err.out <- numeric(length(T))
  rho.out <- numeric(length(T))
  p.out <- numeric(length(T))
  # 20091022 TODO: parallelize this
  for(i in 1:length(T)) {
    states[1] <- Tc[i]
    states[2] <- P[i]
    if(any(is.na(c(Tc[i],P[i])))) {
      # if T or P is NA, all properties are NA
      w <- matrix(rep(NA,23),nrow=1)
      w.out[i,] <- w
      p.out[i] <- NA
      rho.out[i] <- NA
    } else {
      inc <- 0
      h2o <- .Fortran('H2O92',as.integer(specs),as.double(states),
        as.double(rep(0,46)),as.integer(0),PACKAGE='CHNOSZ')
      # errors
      err <- h2o[[4]]
      err.out[i] <- err
      # density
      rho <- h2o[[2]][3]
      rho2 <- h2o[[2]][4]
      if(rho2 > rho) {
        # liquid is denser than vapor
        rho <- rho2 
        # for selecting the liquid properties later
        inc <- 1
      }
      rho.out[i] <- rho
      # most of the properties we're interested in
      w <- t(h2o[[3]][iprop+inc])
      if(err==1) w[1,] <- NA
      # update the ith row of the output matrix
      w.out[i,] <- w
      # Psat
      if(isat | 'psat' %in% tolower(property)) {
        p <- h2o[[2]][2]
        p[p==0] <- NA
        # Psat specifies P=1 below 100 degC
        p[p < 1] <- 1
        p.out[i] <- p
      } else {
        p.out[i] <- P[i]
      }
    }
  }
  # convert output to dataframe
  w.out <- as.data.frame(w.out)
  names(w.out) <- props
  # assemble the properties
  mwH2O <- 18.0152 # SUP92.f
  w.out <- cbind(w.out,V=mwH2O/rho.out,rho=rho.out*1000)
  if(isat | 'psat' %in% tolower(property)) w.out <- cbind(w.out,Psat=p.out)
  # tell the user about any problems
  if(any(err.out==1)) {
    if(length(T) > 1) plural <- "s" else plural <- ""
    nerr <- length(which(err.out==1))
    if(nerr > 1) plural2 <- "s" else plural2 <- ""
    if(isat) msgout(paste("water.SUPCRT92: error",plural2," calculating ",
      nerr," of ",length(T)," point",plural,"; for Psat we need 273.16 < T < 647.067 K\n",sep=""))
    else msgout(paste("water.SUPCRT92: error",plural2," calculating ",nerr,
      " of ",length(T)," point",plural,
      "; T < Tfusion@P, T > 2250 degC, or P > 30kb.\n",sep=""))
      # that last bit is taken from SUP92D.f in the SUPCRT92 distribution
  }
  # if isat is 1, just return the calculated pressures
  if(isat) return(w.out$Psat)
  # return only the selected properties
  icol <- match(tolower(property),tolower(colnames(w.out)))
  return(w.out[,icol,drop=FALSE])
}

water.IAPWS95 <- function(property, T=298.15, P=1, quiet=FALSE) {
  # to get the properties of water via IAPWS-95
  if(!quiet) msgout(paste("water.IAPWS95: calculating", length(T), "values for"))
  M <- 18.015268 # g mol-1
  rho <- function() {
    # return a density in kg m-3
    # corresponding to the given pressure (MPa) and temperature (K)
    pfun <- function(rho,T,P) {
      P <- convert(P,'MPa')
      t <- IAPWS95('p',rho=rho,T=T)[,1] - P
      return(t)
    }
    t <- numeric() 
    for(i in 1:length(T)) {
      if(T[i] < 647.096) {
        rho.lower <- WP02.auxiliary('rho.liquid',T=T[i])-2
        rho.upper <- rho.lower + 400
        if(P[i] < 5000) rho.upper <- rho.lower + 300
        if(P[i] < 1000) rho.upper <- rho.lower + 200
        if(P[i] < 300) rho.upper <- rho.lower + 30
      }
      else { rho.lower <- 0.01; rho.upper <- 1200}
      tu <- try(uniroot(pfun,c(rho.lower,rho.upper),T=T[i],P=P[i])$root,TRUE)
      if(!is.numeric(tu)) {
        warning('water: no root for density between ',round(rho.lower,1),
        ' and ',round(rho.upper,1),' kg m-3 at ',T[i],' K and ',P[i],' bar.',call.=FALSE,immediate.=TRUE)
        tu <- NA
      }
      t <- c(t,tu)
    }
    return(t)
  }
  v <- function() return(M*1000/my.rho)
  p <- function() return(P)
  # Psat stuff
  psat <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      if(T[i] < 373.124) p <- c(p,0.1)
      else p <- c(p,WP02.auxiliary('P.sigma',T[i]))
    }
    return(convert(p,'bar'))
  }
  ## thermodynamic properties
  # convert to SUPCRT reference state
  # at the triple point
  # I2S = SUPCRT - IAPWS ( + entropy in G )
  dH <- -68316.76 - 451.75437
  dS <- 16.7123 - 1.581072
  dG <- -56687.71 + 19.64228 - dS * (T - thermo$opt$Tr)
  # does the reference state used for GHS also go here?
  dU <- -67434.5 - 451.3229
  dA <- -55814.06 + 20.07376 - dS * (T - thermo$opt$Tr)
  # convert IAPWS95() (specific, joule) to (molar, cal) 
  s <- function()
    return(convert(IAPWS95('s',T=T,rho=my.rho)$s*M,'cal')+dS) 
  # u (internal energy) is not here because the letter
  # is used to denote one of the Born functions
  # scratch that! let's put u here and call the other one uborn
  u <- function()
    return(convert(IAPWS95('u',T=T,rho=my.rho)$u*M,'cal')+dU)
  a <- function()
    return(convert(IAPWS95('a',T=T,rho=my.rho)$a*M,'cal')+dA)
  h <- function() 
    return(convert(IAPWS95('h',T=T,rho=my.rho)$h*M,'cal')+dH) 
  g <- function() 
    return(convert(IAPWS95('g',T=T,rho=my.rho)$g*M,'cal')+dG) 
  cv <- function() 
    return(convert(IAPWS95('cv',T=T,rho=my.rho)$cv*M,'cal')) 
  cp <- function() 
    return(convert(IAPWS95('cp',T=T,rho=my.rho)$cp*M,'cal')) 
  w <- function()
    return(IAPWS95('w',T=T,rho=my.rho)$w*100) # to cm/s
  ## electrostatic properties
  epsilon <- function() return(water.AW90(T=T,rho=my.rho,P=convert(P,'MPa')))
  de.dt <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- water.IAPWS95("rho", T=c(t1, t2), P=this.P, quiet=TRUE)[, 1]
      e <- water.AW90(T=c(t1,t2),rho=rho,rep(this.P,2))
      p <- c(p,(e[2]-e[1])/(2*dt))
    }
    return(p)
  }
  de.dp <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- water.IAPWS95("rho", T=this.T, P=c(p1, p2), quiet=TRUE)[, 1]
      e <- water.AW90(P=c(p1,p2),rho=rho,T=rep(this.T,2))
      p <- c(p,(e[2]-e[1])/(2*dp))
    }
    return(p)
  }
  ## Born functions
  q <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- water.IAPWS95("rho", T=rep(this.T, 2), P=c(p1, p2), quiet=TRUE)[, 1]
      e <- water.AW90(T=rep(this.T,2),rho=rho,P=convert(c(p1,p2),'MPa'))
      #p <- c(p,convert(-(1/e[2]-1/e[1])/(2*dp),'cm3bar'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dp))
    }
    return(p)
  }
  n <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- water.IAPWS95("rho", T=rep(this.T, 3), P=c(p1, this.P, p2), quiet=TRUE)[, 1]
      e <- water.AW90(T=rep(this.T,3),rho=rho,P=convert(c(p1,this.P,p2),'MPa'))
      #p <- c(p,convert(convert((-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp,'cm3bar'),'cm3bar'))
      p <- c(p,(-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp)
    }
    return(p)
  }
  y <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- water.IAPWS95("rho", T=c(t1, t2), P=rep(this.P, 2), quiet=TRUE)[, 1]
      e <- water.AW90(T=c(t1,t2),rho=rho,P=convert(rep(this.P,2),'MPa'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dt))
    }
    return(p)
  }
  x <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- water.IAPWS95("rho", T=c(t1, this.T, t2), P=rep(this.P, 3), quiet=TRUE)[, 1]
      e <- water.AW90(T=c(t1,this.T,t2),rho=rho,P=convert(rep(this.P,3),'MPa'))
      p <- c(p,(-(1/e[3]-1/e[2])/dt+(1/e[2]-1/e[1])/dt)/dt)
    }
    return(p)
  }
  uborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; this.T1 <- this.T - dt; this.T2 <- this.T + dt
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho1 <- water.IAPWS95("rho", T=rep(this.T1, 2), P=c(p1, p2), quiet=TRUE)[, 1]
      rho2 <- water.IAPWS95("rho", T=rep(this.T2, 2), P=c(p1, p2), quiet=TRUE)[, 1]
      e1 <- water.AW90(T=rep(this.T1,2),rho=rho1,P=convert(c(p1,p2),'MPa'))
      e2 <- water.AW90(T=rep(this.T2,2),rho=rho2,P=convert(c(p1,p2),'MPa'))
      #p1 <- convert(-(1/e1[2]-1/e1[1])/(2*dp),'cm3bar')
      #p2 <- convert(-(1/e2[2]-1/e2[1])/(2*dp),'cm3bar')
      p1 <- -(1/e1[2]-1/e1[1])/(2*dp)
      p2 <- -(1/e2[2]-1/e2[1])/(2*dp)
      p <- c(p,(p2-p1)/(2*dt))
    }
    return(p)
  }
  ### main loop; init dataframe output and density holders
  w.out <- NULL
  my.rho <- NULL
  # get densities and tell about it
  if(!quiet) msgout(" rho")
  my.rho <- rho() 
  for(i in 1:length(property)) {
    if(property[i] %in% c('e','kt')) {
      # expansivity isn't in the table yet... set it to zero
      warning('water: values of ',property[i],' are NA\n',call.=FALSE)
      inew <- rep(NA,length(T))
    } else {
      if(!quiet) msgout(paste(" ", property[i], sep=""))
      inew <- get(property[i])()
    }
    #if(NA %in% inew) na.h2o <- TRUE
    wnew <- data.frame(inew)
    if(i > 1) w.out <- cbind(w.out,wnew) else w.out <- wnew
  }  
  if(!quiet) msgout("\n")
  return(w.out)
}
