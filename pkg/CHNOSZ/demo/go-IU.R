# CHNOSZ/demo/go-IU.R  20171018
# diagrams using data from the SUPCRTBL compilation
# (BL = Bloomington campus of Indiana University)

## set up plotting area
par(mfrow=c(2, 2))

## start with default database
data(thermo)

###########
### plot 1: boehmite - kaolinite equilibrium
###########
## experimental data from Table 1 of Hemley et al., 1980
# doi:10.2113/gsecongeo.75.2.210
xT <- c(200, 200, 200, 200, 250, 250, 300, 300, 300, 300)
xlogaSiO2 <- -c(2.54, 2.59, 2.65, 2.77, 2.21, 2.32, 1.90, 1.95, 1.94, 1.90)
## set up basis species so that axis.label shows activity of SiO2
basis(c("Al2O3","SiO2", "H2O", "O2"))
T <- 125:350
thermo.plot.new(xlim=range(T), ylim=c(-3.5, -1.5), xlab = axis.label("T"), ylab=axis.label("SiO2"))
points(xT, xlogaSiO2)
basis(delete=TRUE)
## first calculation: CHNOSZ default (SiO2 from SHS89, kaolinite and boehmite from HDNB78)
r1 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T=T, P=1000, exceed.Ttr = TRUE) 
# we need exceed.Ttr = TRUE because the T limit for boehmite is 500 K (Helgeson et al., 1978)
## second calculation: kaolinite from Berman, 1988
Kln_Berman <- info("kaolinite", "cr_Berman")
r2 <- subcrt(c("boehmite", "H2O", "SiO2", Kln_Berman), c(-1, -0.5, -1, 0.5), T=T, P=1000, exceed.Ttr = TRUE) 
## third calculation: boehmite from Hemingway et al., 1991
add.obigt("SUPCRTBL", "boehmite")
r3 <- subcrt(c("boehmite", "H2O", "SiO2", Kln_Berman), c(-1, -0.5, -1, 0.5), T=T, P=1000) 
## fourth calculation: SiO2 from Apps and Spycher, 2004
add.obigt("SUPCRTBL", "SiO2")
r4 <- subcrt(c("boehmite", "H2O", "SiO2", Kln_Berman), c(-1, -0.5, -1, 0.5), T=T, P=1000) 
## log activity of SiO2 is -ve logK
lines(T, -r1$out$logK)
lines(T, -r2$out$logK, lty=2)
lines(T, -r3$out$logK, lty=2, col="red")
lines(T, -r4$out$logK, col="red")
## add labels, legend, and title
text(182.5, -3.17, "SUPCRT92\n(CHNOSZ default)", srt=48, cex=0.7, font=2)
text(147, -3.1, "SUPCRTBL", srt=45.5, cex=0.7, font=2, col="red")
legend("topleft", lty=c(1, 2, 2, 1), col=c("black", "black", "red", "red"), bty="n", cex=0.9,
       legend=c("Kln,Bhm:HDNB78; SiO2:SHS89", "Kln:Ber88", "+ Bhm:HRA91", "+ SiO2:AS04"))
legend("bottomright", pch=1, legend="Hemley et al., 1980", bty="n", cex=0.9)
mtitle(c("Kaolinite - Boehmite", "After Zhu and Lu, 2009 Fig. A1"), cex=0.95)
# doi:10.1016/j.gca.2009.03.015

###########
### plot 2: dawsonite solubility
###########
## experimental data from Benezeth et al., 2007 Table 5
# doi:10.1016/j.gca.2007.07.003
# (averages for each temperature in a single run)
T <- c(100.1, 100.1, 150.1, 100.1, 150.1, 99.8, 99.8, 200.7, 99.8, 50.1, 75.1, 100.3, 150.1)
logK <- -c(14.825, 14.735, 13.625, 14.79, 13.665, 14.725, 14.1775, 12.74, 14.4925, 16.8625, 15.61, 14.51, 13.455)
plot(T, logK, xlim=c(25, 250), ylim=c(-18, -10), xlab=axis.label("T"), ylab=axis.label("logK"))
# this gets us dawsonite and Al(OH)4-
add.obigt("SUPCRTBL")
T <- 0:250
# calculation 1: dawsonite with non-zero Cp
species <- c("dawsonite", "H2O", "Al(OH)4-", "HCO3-", "Na+", "H+")
coeffs <- c(-1, -2, 1, 1, 1, 1)
Daw1 <- subcrt(species, coeffs, T=T)
# calculation 2: dawsonite with 0 Cp
mod.obigt("dawsonite", Cp=0)
Daw2 <- subcrt(species, coeffs, T=T)
## plot the calculated logKs
lines(T, Daw1$out$logK, col="red")
lines(T, Daw2$out$logK, col="red", lty=2)
## add labels, legend, and title
text(182.5, -3.17, "SUPCRT92\n(CHNOSZ default)", srt=43, cex=0.7, font=2)
text(145, -3.1, "SUPCRTBL", srt=41.5, cex=0.7, font=2, col="red")
legend("topleft", lty=1:2, col="red", bty="n", cex=0.9,
       legend=c("Daw Cp != 0", "Daw Cp = 0"))
legend("bottomright", pch=1, legend="Ben\u00e9z\u00e9th et al., 2007", bty="n", cex=0.9)
mtitle(c("Dawsonite - aqueous species", "After Zimmer et al., 2016 Fig. 2"), cex=0.95)
# doi:10.1016/j.cageo.2016.02.013

###########
### plot 3: Eh-pH diagram for As-O-H-S
###########
add.obigt("SUPCRTBL")
#basis(c("Fe", "As", "H2O", "H2S", "H+", "e-"))
#basis(c("Fe", "H2S"), c(-6, -3))
basis(c("As", "H2O", "H2S", "H+", "e-"))
basis(c("H2S"), c(-3))
As_aq <- c("H3AsO4", "H2AsO4-", "HAsO4-2", "AsO4-3", "H3AsO3", "H2AsO3-", "HAsO3-2", "AsO3-3")
AsS_aq <- c("AsS(OH)HS-", "As3S4(HS)2-")
As_cr <- "As"
AsS_cr <- c("realgar,alpha", "realgar,beta", "orpiment", "orpiment,amorphous")
FeAs_cr <- c("arsenopyrite", "scorodite", "ferric arsenate,amorphous")
#species(c(As_aq, AsS_aq, As_cr, AsS_cr, FeAs_cr))
species(c(As_aq, AsS_aq, As_cr, AsS_cr))
species(c(As_aq, AsS_aq), -5)
## a simple diagram, but using only H2S at all pH
#a <- affinity(pH=c(0, 14), Eh=c(-1, 1.5))
#diagram(a)
# the S basis species depends on pH
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
# calculate affinties of formation reactions using the speciated S basis species
res <- 300
# we "blend" the transitions with pH, unlike LZ11's diagram where
# it appears the S-basis species are switched in an on/off fashion
m <- mosaic(bases, pH=c(2, 14, res), Eh=c(-0.6, 0.8, res), blend=TRUE)
# adjust colors and names
fill <- rev(heat.colors(nrow(species())))
fill[11:15] <- "darkgrey"
m$A.species$species$name <- gsub(",alpha", "", m$A.species$species$name)
diagram(m$A.species, fill=fill)
dprop <- describe.property(c("T", "P"), c(25, 1))
legend("topright", legend=dprop, bty="n")
t1 <- quote("As-O-H-S, "~list(sum(S)==10^-3*M, sum(As)==10^-5*M))
t2 <- "After Lu and Zhu, 2011 Fig. 2b"
mtitle(as.expression(c(t1, t2)), cex=0.95)
# doi:10.1007/s12665-010-0652-x

###########
### plot 4: aqueous Al species
###########
add.obigt("SUPCRTBL")
basis(c("Al+3", "F-", "H+", "O2", "H2O"))
Al <- "Al+3"
AlOH <- c("Al(OH)4-", "AlOH+2", "Al(OH)2+", "Al(OH)3")
AlF <- c("AlF+2", "AlF2+", "AlF3", "AlF4-")
AlOHF <- c("Al(OH)2F2-", "Al(OH)2F", "AlOHF2")
species(c(Al, AlOH, AlF, AlOHF), "aq")
a <- affinity(pH=c(0, 10), `F-`=c(-1, -9), T=200)
diagram(a, fill=cm.colors(nrow(species())))
dprop <- describe.property(c("T", "P"), c(200, "Psat"))
legend("topright", legend=dprop, bty="n")
mtitle(c("Aqueous aluminum species",
         "After Tagirov and Schott, 2001 Fig. 4d"), cex=0.95)
# doi:10.1016/S0016-7037(01)00705-0

###########
### clean up: restore thermodynamic database to default
###########
data(thermo)
