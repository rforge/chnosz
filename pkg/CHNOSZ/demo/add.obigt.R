## pe-pH diagram for hydrated iron sulfides,
## goethite and pyrite, after Majzlan et al., 2006
##  Majzlan, J., Navrotsky, A., McClesky, R. B. and Alpers, C. N. (2006) Thermodynamic properties and crystal structure refinement of ferricopiapite, coquimbite, rhomboclase, and Fe2(SO4)3(H2O)5. \emph{Eur. J. Mineral.} \bold{18}, 175--186. \url{http://dx.doi.org/10.1127/0935-1221/2006/0018-0175}
# add some of these species to the database
add.obigt()
basis(c("Fe+2", "SO4-2", "H2O", "H+", "e-"), 
  c(0, log10(3), log10(0.75), 999, 999))
species(c("rhomboclase", "ferricopiapite", "hydronium jarosite",
  "goethite", "melanterite", "pyrite"))
a <- affinity(pH=c(-1, 4, 256), pe=c(-5, 23, 256))
d <- diagram(a, main="Fe-S-O-H, after Majzlan et al., 2006")
# the first four species show up along the top of the diagram
stopifnot(all.equal(unique(t(d$predominant)[256,]), 1:4))
water.lines(yaxis="pe")
text(3, 22, describe.basis(thermo$basis[2:3,], digits=2, oneline=TRUE))
text(3, 21, describe.property(c("T", "P"), c(25, 1), oneline=TRUE))
# reset the database
data(thermo)

## Aqueous Aluminum Species F-/OH-, after Tagirov and Schott, 2001
##  Tagirov, B. and Schott, J. (2001) Aluminum speciation in crustal fluids revisited. \emph{Geochim. Cosmochim. Acta} \bold{65}, 3965--3992. \url{http://dx.doi.org/10.1016/S0016-7037(01)00705-0}
# some of the species are not in the default databse
add.obigt()
# the 999s have no effect on the diagram:
# pH and log_a(F-) are plotting variables
# O2 is not in the formation reactions
# Al+3 is the balanced quantity
basis(c("Al+3", "F-", "H+", "O2", "H2O"), c(rep(999, 4), 0))
species(c("Al+3", "Al(OH)4-", "AlOH+2", "Al(OH)2+", "Al(OH)3",
  "AlF+2", "AlF2+", "AlF3", "AlF4-", "Al(OH)2F2-", "Al(OH)2F",
  "AlOHF2"), "aq")
a <- affinity(pH=c(0, 10), "F-"=c(-1, -9), T=200)
diagram(a, fill="heat")
title(main=paste("Aqueous aluminium species, T=200 C, P=Psat\n",
  "after Tagirov and Schott, 2001"))
# restore thermodynamic database to default
data(thermo)
