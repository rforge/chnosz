\encoding{UTF-8}
\name{berman}
\alias{berman}
\title{Thermodynamic Properties of Minerals}
\description{
Calculate thermodynamic properties of minerals using the equations of Berman (1988).
}

\usage{
  berman(name, T = 298.15, P = 1, thisinfo = NULL, check.G = FALSE,
         calc.transition = TRUE, calc.disorder = FALSE, units = "cal")
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
The \code{name} refers to a mineral that must be listed in \code{thermo$obigt} with the state \samp{cr_Berman}.
This file also holds the chemical formula, which is required for calculating the entropies of the elements in the mineral.
These entropies are used to convert the apparent Gibbs energies from the Berman-Brown convention to the the Benson-Helgeson convention (cf. Anderson, 2005).

Becuase they use a different set of parameters than Helgeson et al., 1978 (see \code{\link{cgl}}), the standard state thermodynamic properties and parameters for the calculations are stored in files under \code{extdata/Berman}.

If \code{check.G} is TRUE, the tabulated value of DGfTrPr is compared with one calculated from DHfPrTr - T*DSPrTr (DS is the difference between the summed entropies of the elements and the tabulated entropy for the mineral).
A warning is produced if the absolute value of the difference between tabulated and calculated DGfTrPr is greater than 1000 J/mol.

Providing \code{thisinfo} means that the mineral name is not searched in \code{thermo$obigt}, potentially saving some running time.

FIXME: the default for \code{calc.disorder} is temporarily FALSE, as the code produces incorrect Gibbs energies.
}

\examples{
\dontshow{data(thermo)}
# other than the formula, the parameters aren't stored in
# thermo$obigt, so this shows NAs
info(info("quartz", "cr_Berman"))
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

# calculate the properties of a "reaction" between
# the Helgeson and Berman versions of quartz
T <- 1000
P <- c(10000, 20000)
subcrt(rep("quartz", 2), c("cr", "cr_Berman"), c(-1, 1), T=T, P=P)

# make a P-T diagram for SiO2 minerals (Ber88 Fig. 4)
basis(c("SiO2", "O2"), c("cr_Berman", "gas"))
species(c("quartz", "quartz,beta", "coesite"), "cr_Berman")
a <- affinity(T=c(200, 1700, 200), P=c(0, 50000, 200))
diagram(a)
}

\references{
Anderson, G. M. (2005) \emph{Thermodynamics of Natural Systems}, 2nd ed., Cambridge University Press, 648 p. \url{http://www.worldcat.org/oclc/474880901}

Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O–K{\s2}O–CaO–MgO–FeO–Fe{\s2}O{\s3}–Al{\s2}O{\s3}–SiO{\s2}–TiO{\s2}–H{\s2}O–CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \url{https://doi.org/10.1093/petrology/29.2.445}

Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.} \bold{278-A}, 1--229. \url{http://www.worldcat.org/oclc/13594862}
}

\concept{Secondary thermodynamic modeling}
