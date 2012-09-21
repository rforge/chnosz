\name{examples}
\alias{examples}
\alias{longex}
\title{Run Examples from the Documentation}
\description{
  Run the examples contained in each of the documentation topics. }

\arguments{
  \item{do.png}{logical, generate PNG files for the plots?}
  \item{which}{character, which example to run.}
}

\usage{
  examples(do.png = FALSE)
  longex(which = c("sources", "NaCl", "copper", "cordierite",
    "phosphate", "nucleobase", "pie", "orp", "findit", 
    "CO2Ac", "nonideal", "TPX"))
}

\details{

  \code{examples} runs all the examples in the documentation for the package. The \code{\link{example}} function is called for each topic with \code{ask} set to \code{FALSE} (so all of the figures are shown without prompting the user). If \code{do.png} is TRUE, the plots in the examples are saved as \code{\link{png}} files having names beginning with the name of each of the help topics.

  \code{longex} is a function that contains code to run various other examples that are not contained in the documentation. The example(s) to run is/are specified by \code{which} (default is to run all of them). The available examples are listed below. See the comments in the source code for more information about each one.

  \tabular{ll}{
    \code{sources} \tab cross-check the reference list with the thermodynamic database \cr
    \code{NaCl} \tab equilibrium constant for aqueous NaCl dissociation (Shock et al., 1992) \cr
    \code{copper} \tab an Eh-pH diagram for the copper-water-glycine system \cr
    \code{cordierite} \tab equilibrium constant of hydrous cordierite dehydration \cr
    \code{phosphate} \tab phosphate speciation with pH, temperature and ionic strength \cr
    \code{nucleobase} \tab relative stabilities of nucleobases and some amino acids \cr
    \code{pie} \tab pie charts comparing relative abundances of organisms (Spear et al., 2005) and model proteins in metastable equilibrium \cr
    \code{orp} \tab oxidation-reduction potential of redox standards as a function of temperature \cr
    \code{findit} \tab detailed example of usage of \code{\link{findit}} using log-normal distribution as an objective \cr
    \code{CO2Ac} \tab activity of CO2 buffered by acetic acid; comparing \code{\link{affinity}(return.buffer=TRUE)} with \code{\link{diagram}(what="CO2")} \cr
    \code{nonideal} \tab activity coefficient of charged species (Alberty, 2003), using the \code{IS} argument of \code{subcrt} \cr
  }

}

\references{
  Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}
  
  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{degrees }C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{http://dx.doi.org/10.1039/FT9928800803}

  Spear, J. R., Walker, J. J., McCollom, T. M. and  Pace, N. R. (2005) Hydrogen and bioenergetics in the Yellowstone geothermal ecosystem. \emph{Proc. Natl. Acad. Sci. U. S. A.} \bold{102}, 2555--2560. \url{http://dx.doi.org/10.1073/pnas.0409574102}
}


\examples{
\dontshow{data(thermo)}
longex(c("orp", "NaCl"))
\dontshow{par(thermo$opar)}
\dontrun{
# use the following to run examples in all help topics
examples()
}
}


\keyword{extra}
