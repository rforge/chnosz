# CHNOSZ/cgl.R
# calculate standard thermodynamic properties of non-aqueous species
# 20060729 jmd

cgl <- function(property = NULL, parameters = NULL, T = 298.15, P = 1) {
  # calculate properties of crystalline, liquid (except H2O) and gas species
  # argument handling
  thermo <- get("thermo")
  Tr <- thermo$opt$Tr
  Pr <- thermo$opt$Pr
  eargs <- eos.args("mk", property = property)
  prop <- eargs$prop
  EOS.Prop <- eargs$Prop

  # initialize output list
  out <- list()
  for(k in 1:nrow(parameters)) {
    # loop over each species
    PAR <- parameters[k, ]
    w <- NULL
    for(i in 1:length(prop)) {
      PROP <- prop[i]
      # a test for availability of the EoS parameters
      # here we assume that the parameters are in the same position as in thermo$obigt
      # leave T transition (in 20th column) alone
      hasEOS <- any(!is.na(PAR[, 13:19]))
      # if at least one of the EoS parameters is available, zero out any NA's in the rest
      if(hasEOS) PAR[, 13:19][, is.na(PAR[, 13:19])] <- 0
      # equations for lambda adapted from HOK+98
      if(PROP == "cp") {
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- PAR$Cp
        else p <- PAR$a + PAR$b * T + PAR$c * T^-2 + PAR$d * T^-0.5 + PAR$e * T^2 + PAR$f * T^PAR$lambda
      }
      if(PROP == "v") {
        p <- rep(PAR$V, length(T))
      }
      if(PROP %in% c("e", "kt")) {
        p <- rep(NA, length(T))
        warning("cgl: E and/or kT of cr, gas and/or liq species are NA.")
      }
      if(PROP == "g") {
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- PAR$G + PAR$Cp * (T - Tr - T * log(T/Tr)) else {
          # Gibbs energy integral: the value at Tref plus heat capacity terms
          p <-   PAR$G + PAR$a * (T - Tr - T * log(T/Tr)) - 
                 PAR$b * (T - Tr)^2 / 2 - PAR$c * (1/T + T/Tr^2 - 2/Tr) / 2 -
                 PAR$d * (T^0.5 - 0.5 * T * Tr^-0.5 - 0.5 * Tr^0.5) / -0.25 -
                 PAR$e * (T^3 - 3 * T * Tr^2 + 2 * Tr^3) / 6
        }
        # entropy and volume terms
        if(!is.na(PAR$S)) p <- p - PAR$S * (T - Tr)
        if(!is.na(PAR$V)) p <- p + convert(PAR$V * (P - Pr), "calories")
        # use additional heat capacity term if it's defined
        if(!is.na(PAR$f) & !is.na(PAR$lambda)) if(PAR$f != 0) {
          if(PAR$lambda == -1) p <- p + PAR$f * (log(T/Tr) - T * (1/Tr - 1/T))
          else p <- p + PAR$f * ( T^(PAR$lambda + 1) - (PAR$lambda + 1) * T * Tr^PAR$lambda + 
            PAR$lambda * Tr^(PAR$lambda + 1) ) / ( PAR$lambda * (PAR$lambda + 1) ) 
        }
      }
      if(PROP == "h") { 
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- PAR$Cp * (T - Tr) else {
          p <- PAR$a * (T - Tr) + PAR$b * (T^2 - Tr^2) / 2 +
               PAR$c * (1/T - 1/Tr) / -1 + PAR$d * (T^0.5 - Tr^0.5) / 0.5 + 
               PAR$e * (T^3 - Tr^3) / 3 
               # SUPCRT seems to ignore this term? ... 20070802
               # + convert(PAR$V*(P-Pr),'calories')
        }
        p <- PAR$H + p
        if(!is.na(PAR$f) & !is.na(PAR$lambda)) if(PAR$f != 0) {
           if(PAR$lambda == -1) p <- p + PAR$f * log(T/Tr) 
           else p <- p - PAR$f * ( T^(PAR$lambda + 1) - Tr^(PAR$lambda + 1) ) / (PAR$lambda + 1)
        }
      }
      if(PROP=="s") {
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- PAR$Cp * log(T/Tr) else
        p <- PAR$a * log(T / Tr) + PAR$b * (T - Tr) + 
             PAR$c * (T^-2 - Tr^-2) / -2 + PAR$e * (T^2 - Tr^2) / 2 + 
             PAR$d * (T^-0.5 - Tr^-0.5) / -0.5
        p <- PAR$S + p
        if(!is.na(PAR$f) & !is.na(PAR$lambda)) if(PAR$f != 0) {
          p <- p + PAR$f * (T^PAR$lambda - Tr^PAR$lambda) / PAR$lambda
        }
      }
      wnew <- data.frame(p)
      if(i > 1) w <- cbind(w, wnew) else w <- wnew
    }
  colnames(w) <- EOS.Prop
  out[[k]] <- w
 }
 return(out)
}

