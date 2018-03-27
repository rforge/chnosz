\encoding{UTF-8}
\name{berman}
\alias{berman}
\title{Thermodynamic Properties of Minerals}
\description{
Calculate thermodynamic properties of minerals using the equations of Berman (1988).
}

\usage{
  berman(name, T = 298.15, P = 1, thisinfo = NULL, check.G = FALSE,
         calc.transition = TRUE, calc.disorder = TRUE, units = "cal")
}

\arguments{
  \item{name}{character, name of mineral}
  \item{T}{numeric, temperature(s) at which to calculate properties (K)}
  \item{P}{numeric, pressure(s) at which to calculate properties (bar)}
  \item{thisinfo}{dataframe, row for mineral from \code{thermo$obigt}}
  \item{check.G}{logical, check consistency of G, H, and S?}
  \item{calc.transition}{logical, include calculation of polymorphic transition properties?}
  \item{calc.disorder}{logical, include calculation of disordering properties?}
  \item{units}{character, energy units, \samp{cal} or \samp{J}}
}

\details{
This function calculates the thermodynamic properties of minerals at high \P and \T using equations given by Berman (1988).
The standard state thermodynamic properties and parameters for the calculations are stored in files under \code{extdata/Berman}, or can be read from a user-created file (if available) named \samp{berman.csv} in the working directory.
If the function is called with missing \code{name}, the data frame of parameters for all available minerals is returned.

When \code{name} is given, it refers to a mineral that should listed in \code{thermo$obigt} with the state \samp{cr} and the chemical formula, and optinally an abbreviation and references, but all other properties set to NA.
The chemical formula is required for calculating the entropies of the elements in the mineral, which are are used to convert the apparent Gibbs energies from the Berman-Brown convention to the the Benson-Helgeson convention (cf. Anderson, 2005).

If \code{check.G} is TRUE, the tabulated value of DGfTrPr is compared with one calculated from DHfPrTr - T*DSPrTr (DS is the difference between the summed entropies of the elements and the tabulated entropy for the mineral).
A warning is produced if the absolute value of the difference between tabulated and calculated DGfTrPr is greater than 1000 J/mol.

Providing \code{thisinfo} means that the mineral name is not searched in \code{thermo$obigt}, potentially saving some running time.
}

\examples{
\dontshow{data(thermo)}
# other than the formula, the parameters aren't stored in
# thermo$obigt, so this shows NAs
info(info("quartz", "cr"))
# properties of alpha-quartz (aQz) at 298.15 K and 1 bar
berman("quartz")
# Gibbs energies of aQz and coesite at higher T and P
T <- seq(200, 1300, 100)
P <- seq(22870, 31900, length.out=length(T))
G_aQz <- berman("quartz", T=T, P=P)$G
G_Cs <- berman("coesite", T=T, P=P)$G
# that is close to the univariant curve (Ber88 Fig. 4),
# so the difference in G is close to 0
DGrxn <- G_Cs - G_aQz
stopifnot(all(abs(DGrxn) < 100))

# make a P-T diagram for SiO2 minerals (Ber88 Fig. 4)
basis(c("SiO2", "O2"), c("cr", "gas"))
species(c("quartz", "quartz,beta", "coesite"), "cr")
a <- affinity(T=c(200, 1700, 200), P=c(0, 50000, 200))
diagram(a)

## a longer example, comparing diagrams made using the
## Berman and Helgeson datasets, after Sverjensky et al., 1991
res <- 200
# using the Helgeson data
add.obigt("SUPCRT92")
# set up basis species
basis(c("K+", "Al+3", "quartz", "H2O", "O2", "H+"))
# use pH = 0 so that aK+ = aK+/aH+
basis("pH", 0)
# load the species
species(c("K-feldspar", "muscovite", "kaolinite",
          "pyrophyllite", "andalusite"), "cr")
# calculate affinities in aK+ - temperature space
a <- affinity(`K+`=c(0, 5, res), T=c(200, 650, res), P=1000)
# note that we go just past the quartz transition,
# but it has no effect on the diagram
diagram(a, xlab=ratlab("K+"))
# now using the Berman data: this resets the thermodynamic database
# without affecting the basis and species settings
data(OBIGT)
# it might be good to check that we have Berman's quartz
# and not coesite or some other SiO2 phase
info(basis()$ispecies[3])
a <- affinity(`K+`=c(0, 5, res), T=c(200, 650, res), P=1000)
diagram(a, add=TRUE, names="", col="blue", lwd=2)
legend("topleft", lty=c(1, 1, NA), lwd=c(1, 2, 0), col=c("black", "blue", ""),
       legend=c("Helgeson et al., 1978 (unadjusted)",
       "Berman, 1988", "    (adjusted by Sverjensky et al., 1991)"), bty="n")
title(main="Comparison of Helgeson and Berman datasets at 1000 bar")
}

\references{
Anderson, G. M. (2005) \emph{Thermodynamics of Natural Systems}, 2nd ed., Cambridge University Press, 648 p. \url{http://www.worldcat.org/oclc/474880901}

Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \url{https://doi.org/10.1093/petrology/29.2.445}

Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.} \bold{278-A}, 1--229. \url{http://www.worldcat.org/oclc/13594862}

Sverjensky, D. A., Hemley, J. J. and D'Angelo, W. M. (1991) Thermodynamic assessment of hydrothermal alkali feldspar-mica-aluminosilicate equilibria. \emph{Geochim. Cosmochim. Acta} \bold{55}, 989-1004. \url{https://doi.org/10.1016/0016-7037(91)90157-Z}
}

\concept{Thermodynamic calculations}
