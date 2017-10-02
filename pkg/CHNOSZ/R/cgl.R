# CHNOSZ/cgl.R
# calculate standard thermodynamic properties of non-aqueous species
# 20060729 jmd

cgl <- function(property = NULL, parameters = NULL, T = 298.15, P = 1) {
  # calculate properties of crystalline, liquid (except H2O) and gas species
  Tr <- 298.15
  Pr <- 1
  # the number of T, P conditions
  ncond <- max(c(length(T), length(P)))
  # initialize output list
  out <- list()
  # loop over each species
  for(k in 1:nrow(parameters)) {
    # the parameters for *this* species
    PAR <- parameters[k, ]
    if(PAR$state=="cr_Berman") {
      # use Berman equations (parameters not in thermo$obigt)
      properties <- berman(PAR$name, T=T, P=P, thisinfo=PAR)
      iprop <- match(property, colnames(properties))
      values <- properties[, iprop, drop=FALSE]
    } else {
      # start with NA values
      values <- data.frame(matrix(NA, ncol = length(property), nrow=ncond))
      colnames(values) <- property
      # additional calculations for quartz and coesite
      qtz <- quartz_coesite(PAR, T, P)
      isqtz <- !identical(qtz$V, 0)
      for(i in 1:length(property)) {
        PROP <- property[i]
        # a test for availability of the EoS parameters
        # here we assume that the parameters are in the same columns as in thermo$obigt
        # leave T transition (in 20th column) alone
        hasEOS <- any(!is.na(PAR[, 13:19]))
        # if at least one of the EoS parameters is available, zero out any NA's in the rest
        if(hasEOS) PAR[, 13:19][, is.na(PAR[, 13:19])] <- 0
        # equations for lambda adapted from HOK+98
        if(PROP == "Cp") {
          # use constant Cp if the EoS parameters are not available
          if(!hasEOS) p <- PAR$Cp
          else p <- PAR$a + PAR$b * T + PAR$c * T^-2 + PAR$d * T^-0.5 + PAR$e * T^2 + PAR$f * T^PAR$lambda
        }
        if(PROP == "V") {
          if(isqtz) p <- qtz$V
          else p <- rep(PAR$V, ncond)
        }
        if(PROP %in% c("E", "kT")) {
          p <- rep(NA, ncond)
          warning("cgl: E and/or kT of cr, gas and/or liq species are NA.")
        }
        if(PROP == "G") {
          # use constant Cp if the EoS parameters are not available
          if(!hasEOS) p <- PAR$Cp * (T - Tr - T * log(T/Tr)) else {
            # Gibbs energy integral: the value at Tref plus heat capacity terms
            p <-   PAR$a * (T - Tr - T * log(T/Tr)) - 
                   PAR$b * (T - Tr)^2 / 2 - PAR$c * (1/T + T/Tr^2 - 2/Tr) / 2 -
                   PAR$d * (T^0.5 - 0.5 * T * Tr^-0.5 - 0.5 * Tr^0.5) / -0.25 -
                   PAR$e * (T^3 - 3 * T * Tr^2 + 2 * Tr^3) / 6
          }
          # use additional heat capacity term if it's defined
          if(!is.na(PAR$f) & !is.na(PAR$lambda)) if(PAR$f != 0) {
            if(PAR$lambda == -1) p <- p + PAR$f * (log(T/Tr) - T * (1/Tr - 1/T))
            else p <- p + PAR$f * ( T^(PAR$lambda + 1) - (PAR$lambda + 1) * T * Tr^PAR$lambda + 
              PAR$lambda * Tr^(PAR$lambda + 1) ) / ( PAR$lambda * (PAR$lambda + 1) ) 
          }
          # entropy and volume terms
          if(!is.na(PAR$S)) p <- p - PAR$S * (T - Tr)
          if(isqtz) p <- p + qtz$G
          else if(!is.na(PAR$V)) p <- p + convert(PAR$V * (P - Pr), "calories")
          p <- PAR$G + p
        }
        if(PROP == "H") { 
          # use constant Cp if the EoS parameters are not available
          if(!hasEOS) p <- PAR$Cp * (T - Tr) else {
            p <- PAR$a * (T - Tr) + PAR$b * (T^2 - Tr^2) / 2 +
                 PAR$c * (1/T - 1/Tr) / -1 + PAR$d * (T^0.5 - Tr^0.5) / 0.5 + 
                 PAR$e * (T^3 - Tr^3) / 3 
          }
          if(!is.na(PAR$f) & !is.na(PAR$lambda)) if(PAR$f != 0) {
             if(PAR$lambda == -1) p <- p + PAR$f * log(T/Tr) 
             else p <- p - PAR$f * ( T^(PAR$lambda + 1) - Tr^(PAR$lambda + 1) ) / (PAR$lambda + 1)
          }
          if(isqtz) p <- p + qtz$H
          ## SUPCRT seems to ignore this term? ... 20070802
          #else p <- p + convert(PAR$V*(P-Pr),'calories')
          p <- PAR$H + p
        }
        if(PROP=="S") {
          # use constant Cp if the EoS parameters are not available
          if(!hasEOS) p <- PAR$Cp * log(T/Tr) else {
            p <- PAR$a * log(T / Tr) + PAR$b * (T - Tr) + 
                 PAR$c * (T^-2 - Tr^-2) / -2 + PAR$e * (T^2 - Tr^2) / 2 + 
                 PAR$d * (T^-0.5 - Tr^-0.5) / -0.5
          }
          if(!is.na(PAR$f) & !is.na(PAR$lambda)) if(PAR$f != 0) {
            p <- p + PAR$f * (T^PAR$lambda - Tr^PAR$lambda) / PAR$lambda
          }
          p <- PAR$S + p + qtz$S
        }
        values[, i] <- p
      }
    } # end calculations using parameters from thermo$obigt
    out[[k]] <- values
  } # end loop over species
  return(out)
}

### unexported function ###

# calculate GHS and V corrections for quartz and coesite 20170929
# (these are the only mineral phases which SUPCRT92 applies a variable volume)
quartz_coesite <- function(PAR, T, P) {
  # the corrections are 0 for anything other than quartz and coesite
  if(!PAR$name %in% c("quartz", "coesite")) return(list(G=0, H=0, S=0, V=0))
  ncond <- max(c(length(T), length(P)))
  # Tr, Pr and TtPr (transition temperature at Pr)
  Pr <- 1      # bar
  Tr <- 298.15 # K
  TtPr <- 848  # K
  # constants from SUP92D.f
  aa <- 549.824
  ba <- 0.65995
  ca <- -0.4973e-4
  VPtTta <- 23.348
  VPrTtb <- 23.72
  Stran <- 0.342
  # constants from REAC92D.f
  VPrTra <- 22.688 # VPrTr(a-quartz)
  Vdiff <- 2.047   # VPrTr(a-quartz) - VPrTr(coesite)
  #k <- 38.5       # dPdTtr(a/b-quartz)
  k <- 38.45834    # calculated in CHNOSZ: dPdTtr(info("quartz"))
  # code adapted from REAC92D.f
  qphase <- gsub("cr", "", PAR$state)
  if(qphase == 2) {
    Pstar <- P
    Sstar <- rep(0, ncond)
    V <- rep(VPrTtb, ncond)
  } else {
    Pstar <- Pr + k * (T - TtPr)
    Sstar <- rep(Stran, ncond)
    V <- VPrTra + ca*(P-Pr) + (VPtTta - VPrTra - ca*(P-Pr))*(T-Tr) / (TtPr + (P-Pr)/k - Tr)
  }
  Pstar[T < TtPr] <- Pr
  Sstar[T < TtPr] <- 0
  if(PAR$name == "coesite") {
    VPrTra <- VPrTra - Vdiff
    VPrTtb <- VPrTtb - Vdiff
    V <- V - Vdiff
  }
  GVterm <- convert(1, "calories") * (VPrTra * (P - Pstar) + VPrTtb * (Pstar - Pr) -
    0.5 * ca * (2 * Pr * (P - Pstar) - (P^2 - Pstar^2)) -
    ca * k * (T - Tr) * (P - Pstar) +
    k * (ba + aa * ca * k) * (T - Tr) * log((aa + P/k) / (aa + Pstar/k)))
  SVterm <- convert(1, "calories") * (-k * (ba + aa * ca * k) *
    log((aa + P/k) / (aa + Pstar/k)) + ca * k * (P - Pstar)) - Sstar
  list(G=GVterm, S=SVterm, H=GVterm + T*SVterm, V=V)
}

