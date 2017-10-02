# CHNOSZ/berman.R 20170930
# calculate thermodynamic properties of minerals using Berman formulation

berman <- function(name, T = 298.15, P = 1) {
  # reference temperature and pressure
  Pr <- 1
  Tr <- 298.15
  file <- system.file("extdata/Berman/Ber88.csv", package="CHNOSZ")
  dat <- read.csv(file, as.is=TRUE)
  # remove the multipliers
  multexp <- c(0, 0, 0, 0,          # Ber88 Table 2
               0, -2, -5, -7,             # Table 3a
               6, 12, 6, 10,              # Table 4
               0, 0, 0, 2, 5, 0,          # Table 3b
               0, 0, 0, -3, -5, 2, 6, -4  # Table 5
               )
  dat[, 2:27] <- t(t(dat[, 2:27]) / 10^multexp)
  # which row has data for this mineral?
  irow <- which(dat$name == name)
  # only the immediately following assign() call is needed for the function to work,
  # but an explicit dummy assignment here is used to avoid "Undefined global functions or variables" in R CMD check
  GfPrTr <- HfPrTr <- SPrTr <- Tmax <- Tmin <- VPrTr <-
    d0 <- d1 <- d2 <- d3 <- d4 <- d5 <- k0 <- k1 <- k2 <- k3 <- v1 <- v2 <- v3 <- v4 <- NA
  # assign values to the variables used below
  for(i in 1:ncol(dat)) assign(colnames(dat)[i], dat[irow, i])
  # check that G in data file is the G of formation from the elements --> Benson-Helgeson convention (DG = DH - T*DS)
  # we get the entropy of the elements using the chemical formula in thermo$obigt
  iname <- info(name, "cr_Berman", check.it=FALSE)
  SPrTr_elements <- convert(entropy(info(iname)$formula), "J")
  GfPrTr_calc <- HfPrTr - Tr * (SPrTr - SPrTr_elements)
  Gdiff <- GfPrTr_calc - GfPrTr
  if(abs(Gdiff) >= 1000) warning(paste0(name, ": GfPrTr(calc) - GfPrTr(table) is too big! == ",
                                        round(GfPrTr_calc - GfPrTr), " J/mol"), call.=FALSE)
  # (the tabulated GfPrTr is unused below)

  ### thermodynamic properties ###
  # calculate Cp and V (Berman, 1988 Eqs. 4 and 5)
  Cp <- k0 + k1 * T^-0.5 + k2 * T^-2 + k3 * T^-3
  P_Pr <- P - Pr
  T_Tr <- T - Tr
  V <- VPrTr * (1 + v1 * P_Pr + v2 * P_Pr^2 + v3 * T_Tr + v4 * T_Tr^2)
  # calculate Ga (Ber88 Eq. 6) --> Berman-Brown convention (DG = DH - T*S)
  Ga <- HfPrTr - T * SPrTr + k0 * ( (T - Tr) - T * (log(T) - log(Tr)) ) +
    2 * k1 * ( (T^0.5 - Tr^0.5) + T*(T^-0.5 - Tr^-0.5) ) -
    k2 * ( (T^-1 - Tr^-1) - T / 2 * (T^-2 - Tr^-2) ) -
    k3 * ( (T^-2 - Tr^-2) / 2 - T / 3 * (T^-3 - Tr^-3) ) +
    VPrTr * ( (v1 / 2 - v2) * (P^2 - Pr^2) + v2 / 3 * (P^3 - Pr^3) +
      (1 - v1 + v2 + v3 * (T - Tr) + v4 * (T - Tr)^2) * (P - Pr) )
  # calculate Ha (symbolically integrated using sympy - expressions not simplified)
  intCp <- T*k0 - Tr*k0 + k2/Tr - k2/T + k3/(2*Tr^2) - k3/(2*T^2) + 2.0*k1*T^0.5 - 2.0*k1*Tr^0.5
  intVminusTdVdT <- -VPrTr + P*(VPrTr + VPrTr*v2 - VPrTr*v1 - Tr*VPrTr*v3 + VPrTr*v4*Tr^2 - VPrTr*v4*T^2) +
    P^2*(VPrTr*v1/2 - VPrTr*v2) + VPrTr*v1/2 - VPrTr*v2/3 + Tr*VPrTr*v3 + VPrTr*v4*T^2 - VPrTr*v4*Tr^2 + VPrTr*v2*P^3/3
  Ha <- HfPrTr + intCp + intVminusTdVdT
  # calculate S (also symbolically integrated)
  intCpoverT <- k0*log(T) - k0*log(Tr) - k3/(3*T^3) + k3/(3*Tr^3) + k2/(2*Tr^2) - k2/(2*T^2) + 2.0*k1*Tr^-0.5 - 2.0*k1*T^-0.5
  intdVdT <- -VPrTr*(v3 + v4*(-2*Tr + 2*T)) + P*VPrTr*(v3 + v4*(-2*Tr + 2*T))
  S <- SPrTr + intCpoverT - intdVdT

  ### disorder thermodynamic properties ###
  if(!is.na(Tmin) & !is.na(Tmax) & any(T > Tmin)) {
    # starting disorder contributions are 0
    Cpds <- Hds <- Sds <- Vds <- Gds <- 0
    # the lower integration limit is Tmin
    iTds <- T > Tmin
    Tds <- T[iTds]
    # the upper integration limit is Tmax
    Tds[Tds > Tmax] <- Tmax
    # Ber88 Eqs. 15, 16, 17, 18, 19
    Cpds[iTds] <- d0 + d1 * Tds^-0.5 + d2 * Tds^-2 + d3 * Tds + d4 * Tds^2
    Hds[iTds] <- d0 * (Tds - Tmin) + 2 * d1 * (Tds^-0.5 - Tmin^-0.5) -
      d2 * (Tds^-1 - Tmin^-1) + d3 * (Tds^2 - Tmin^2) / 2 + d4 * (Tds^3 - Tmin^3) / 3
    Sds[iTds] <- d0 * (log(Tds) - log(Tmin)) - 2 * d1 * (Tds^-0.5 - Tmin^-0.5) -
      d2 * (Tds^-2 - Tmin^-2) / 2 + d3 * (Tds - Tmin) + d4 * (Tds^2 - Tmin^2) / 2
    # we can't do this if d5 == 0 (dolomite and gehlenite)
    if(d5 != 0) Vds <- Hds / d5
    Gds <- Hds - T * Sds + Vds * (P - Pr)
    # Gds above Tmax (Eq. 20)
    ihigh <- T > Tmax
    # note that Gds[ihigh] and Sds[ihigh] on the rhs were both calculated at Tmax (above)
    Gds[ihigh] <- Gds[ihigh] - (T[ihigh] - Tmax) * Sds[ihigh]
    # apply the disorder contributions
    Ga <- Ga + Gds
    Ha <- Ha + Hds
    S <- S + Sds
    V <- V + Vds
    Cp <- Cp + Cpds
  }

  ### (for testing) use G = H - TS to check that integrals for H and S are written correctly
  Ga_fromHminusTS <- Ha - T * S
  if(!all.equal(Ga_fromHminusTS, Ga)) stop("incorrect integrals detected using DG = DH - T*S")

  ### thermodynamic and unit conventions used in SUPCRT ###
  # use entropy of the elements in calculation of G --> Benson-Helgeson convention (DG = DH - T*DS)
  Gf <- Ga + Tr * SPrTr_elements
  # convert J to cal
  G <- convert(Gf, "cal")
  H <- convert(Ha, "cal")
  S <- convert(S, "cal")
  Cp <- convert(Cp, "cal")
  # convert J/bar to cm^3/mol
  V <- V * 10

  data.frame(T=T, P=P, G=G, H=H, S=S, Cp=Cp, V=V)
}
