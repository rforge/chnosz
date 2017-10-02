# demo for the Deep Earth Water (DEW) model in CHNOSZ 20170927

# set up subplots
par(mfrow = c(2, 2), mar=c(3.0, 3.5, 2.5, 1.0), mgp=c(1.7, 0.3, 0), las=1, tcl=0.3, xaxs="i", yaxs="i")

# activate DEW model
oldwat <- water("DEW")

#### plot 1: quartz solubility at high pressure
## after Figure 7D of Sverjensky et al., 2014 
## (Geochim. Cosmochim. Acta, https://doi.org/10.1016/j.gca.2013.12.019)

# load SiO2 and Si2O4 data taken from DEW spreadsheet
iSi <- add.obigt("DEW_aq", c("SiO2", "Si2O4"))
# print the data references to confirm we got the right ones
thermo.refs(iSi)
# set temperature ranges for different pressures
# data.frame is used to make P and T the same length
PT0.5 <- data.frame(P=500, T=seq(200, 550, 10))
PT1.0 <- data.frame(P=1000, T=seq(200, 700, 10))
PT2.0 <- data.frame(P=2000, T=seq(200, 700, 10))
PT5.0 <- data.frame(P=5000, T=seq(200, 850, 10))
PT10.0 <- data.frame(P=10000, T=seq(200, 825, 10))
PT20.0 <- data.frame(P=20000, T=seq(200, 800, 10))
PT <- rbind(PT0.5, PT1.0, PT2.0, PT5.0, PT10.0, PT20.0)
# reaction 1: quartz = SiO2(aq) [equivalent to quartz + 3 H2O = Si(OH)4]
SiO2_logK <- subcrt(c("quartz", "SiO2"), c("cr_Berman", "aq"), c(-1, 1), P=PT$P, T=PT$T)$out$logK
# reaction 2: 2 quartz = Si2O4(aq) [equivalent to 2 quartz + 3 H2O = Si2O(OH)6]
Si2O4_logK <- subcrt(c("quartz", "Si2O4"), c("cr_Berman", "aq"), c(-2, 1), P=PT$P, T=PT$T)$out$logK
# plot the sum of molalities (== activities) for each pressure
plot(c(200, 1000), c(-2.5, 0.5), type="n", xlab=axis.label("T"), ylab="log molality")
for(P in unique(PT$P)) {
  icond <- PT$P == P
  SiO2_logm <- SiO2_logK[icond]
  Si2O4_logm <- Si2O4_logK[icond]
  logm <- log10(10^SiO2_logm + 10^Si2O4_logm)
  lines(PT$T[icond], logm)
  # add text label
  lastT <- tail(PT$T[icond], 1)
  Pkb <- paste(format(P/1000, nsmall=1), "kb")
  text(lastT+25, tail(logm, 1), Pkb, adj=0)
}
t1 <- quote("Solubility of"~alpha*"-quartz")
t2 <- "(after Sverjensky et al., 2014)"
mtitle(as.expression(c(t1, t2)))
# TODO: lines are a little low at highest P and P ...
# does the Berman, 1988 quartz data increase high-PT solubilities?

#### plot 2: correlations between non-solvation volume and HKF a1 parameter
## after Figures 12B and 12C of Sverjensky et al., 2014
## (Geochim. Cosmochim. Acta, https://doi.org/10.1016/j.gca.2013.12.019)
# load the fitted parameters for species as used by SHA14
# TODO: also use their Ca+2??
# NOTE: don't load NaCl, NH4+, or HS- here because the DEW spreadsheet lists a1 from the correlation
add.obigt("DEW", c("CO3-2", "BO2-", "MgCl+", "SiO2", "HCO3-", "Si2O4"))
# set up the plot
V0nlab <- expression(Delta * italic(V) * degree[n]~~(cm^3~mol^-1))
a1lab <- expression(italic(a)[1]%*%10~~(cal~mol~bar^-1))
plot(c(-25, 50), c(-4, 12), type="n", xlab=V0nlab, ylab=a1lab)
# a function to get the HKF parameters, calculate nonsolvation volume, plot points, labels, error bars, and correlation lines
plotfun <- function(species, col, pch, cex, dy, error, xlim, corrfun) {
  # get HKF parameters
  par <- info(info(species))
  a1 <- par$a1 * 10
  # get the nonsolvation volume
  Vn <- unlist(hkf("V", par, contrib="n")$aq)
  points(Vn, a1, col=col, pch=pch, cex=cex)
  for(i in 1:length(species)) text(Vn[i], a1[i]+dy, expr.species(species[i]))
  arrows(Vn, a1 - error, Vn, a1 + error, length = 0.03, angle = 90, code = 3, col=col)
  lines(xlim, corrfun(xlim), col=col)
}
# monovalent ions: Na+, K+, Cl-, Br-
monofun <- function(Vn) 2.0754 + 0.10871 * Vn
# for easier reading, set y-offset to NA so the labels aren't plotted
plotfun(c("Na+", "K+", "Cl-", "Br-"), "red", 19, 1, NA, 0.5, c(-7, 35), monofun)
# divalent ions: Mg+2, Ca+2, CO3-2, SO4-2
difun <- function(Vn) 3.5321 + 0.23911 * Vn
plotfun(c("Mg+2", "Ca+2", "CO3-2", "SO4-2"), "black", 15, 1, 1.2, 0.7, c(-20, 25), difun)
# complexes and neutral molecules: BO2-, MgCl+, SiO2, NaCl, HCO3-, Si2O4, NH4+, HS-
compfun <- function(Vn) 1.5204 + 0.19421 * Vn
plotfun(c("MgCl+", "SiO2", "NaCl", "HCO3-", "Si2O4"), "blue1", 18, 1.5, 1, 0.5, c(-20, 50), compfun)
# for easier reading, put some labels below the points
plotfun(c("BO2-", "NH4+", "HS-"), "blue1", 18, 1.5, -1.2, 0.5, c(-20, 50), compfun)
# include an empty subscript for better spacing between the lines
t1 <- quote("Correlations between non-solvation"[])
t2 <- quote("volume and HKF "*italic(a)[1]*" parameter")
mtitle(as.expression(c(t1, t2)))

#### plots 3 and 4: aqueous inorganic and organic carbon species at high pressure
## after Figure 1b of Sverjensky et al., 2014
## (Nature Geoscience, https://doi.org/10.1038/NGEO2291)

# define system with loga.species = 0
basis("CHNOS+")
species(c("CO2", "HCO3-", "CO3-2", "acetic acid", "acetate", "CH4"))
species(1:6, 0)

# a function to make the diagrams
dfun <- function(T = 600, P = 50000, res=300) {
  a <- affinity(pH = c(0, 10, res), O2 = c(-24, -12, res), T = T, P = P)
  diagram(a, limit.water = FALSE, fill=tail(topo.colors(7), -1))
  dp <- describe.property(c("     T", "     P"), c(T, P), digits=0)
  legend("bottomleft", legend=dp, bty="n")
}

# first plot: CHNOSZ default database
data(OBIGT)
dfun()
t1 <- quote("CHNOSZ default database"[])
t2 <- quote("(not recommended for high"~italic(P)*")")
mtitle(as.expression(c(t1, t2)))
# second plot: use CO2, HCO3-, CO3-2, and methane data from DEW spreadsheet
add.obigt("DEW_aq", c("CO2", "HCO3-", "CO3-2", "methane"))
dfun()
CO2quote <- quote(list(CO[2], HCO[3]^"-", CO[3]^"-2"))
DEWexpr <- substitute("DEW data for"~x, list(x=CO2quote))
mtitle(as.expression(c(DEWexpr, "and methane")))

# reset the database and previous water computational option
data(OBIGT)
water(oldwat)
