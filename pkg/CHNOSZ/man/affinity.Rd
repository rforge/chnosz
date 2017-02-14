\encoding{UTF-8}
\name{affinity}
\alias{affinity}
\title{Chemical Affinities of Formation Reactions}
\description{
  Calculate the chemical affinities of formation reactions of species. Do it for single values of temperature, pressure, ionic strength and chemical activities of the basis species, or as a function of one or more of these variables. Or, return other properties including standard molal Gibbs energies of basis species and species of interest, and standard molal Gibbs energies, equilibrium constants and activity products of formation reactions. 
}

\usage{
  affinity(..., property=NULL, sout=NULL, exceed.Ttr=FALSE,
    return.buffer=FALSE, balance="PBB", iprotein=NULL, loga.protein=-3)
}

\arguments{
  \item{...}{numeric, zero or more named arguments, used to identify the variables of interest in the calculations}
  \item{property}{character, denoting the property to be calculated. Default is \samp{A}, for chemical affinity of formation reactions of species of interest}
  \item{sout}{list, output from subcrt function}
  \item{exceed.Ttr}{logical, allow subcrt to compute properties for phases beyond their transition temperature?}
  \item{return.buffer}{logical. If \code{TRUE}, and a \code{\link{buffer}} has been associated with one or more basis species in the system, return the values of the activities of the basis species calculated using the buffer (it is not necessary in this case to have defined any species of interest). Default is \code{FALSE}}
  \item{balance}{character. This argument is passed to \code{\link{buffer}} to identify a conserved basis species (or \samp{PBB}) in a chemical activity buffer. Default is \samp{PBB}}
  \item{iprotein}{numeric, indices of proteins in \code{\link{thermo}$protein} for which to calculate properties}
  \item{loga.protein}{numeric, logarithms of activities of proteins identified in \code{iprotein}}
}

\details{
\code{affinity} calculates the chemical affinities of reactions to form the species of interest from the basis species.
The calculation of chemical affinities relies on the current definitions of the \code{\link{basis}} species and \code{\link{species}} of interest.
It is possible to use the results of \code{affinity} to generate equilibrium activity diagrams using \code{\link{diagram}}.

The equation used to calculate chemical affinity \emph{\bold{A}} is \emph{\bold{A}}=\eqn{RT\ln (K/Q)}{RT*ln(K/Q)} (Kondepudi and Prigogine, 1998), where \eqn{K} denotes the equilibrium constant of the reaction and \eqn{Q} stands for the activity product of the species in the reaction.
(The equilibrium constant is related to standard Gibbs energy of reaction, \eqn{\Delta G^{\circ}_r}{DeltaG0r}, by \eqn{\Delta G^{\circ}_r = -2.303RT\log K}{DeltaG0r = -2.303*RT*logK}, where \eqn{R} and \eqn{T} stand for, respectively, the gas constant and temperature).
With the approach of a given reaction to a state of equilibrium, the chemical affinity tends toward zero, or \eqn{K = Q}. 

Valid properties are \samp{A} or NULL for chemical affinity, \samp{logK} or \samp{logQ} for logarithm of equilibrium constant and reaction activity product, or any of the properties available in \code{\link{subcrt}} except for \samp{rho}.
The properties returned are those of the formation reactions of the species of interest from the basis species.
It is also possible to calculate the properties of the species of interest themselves (not their formation reactions) by setting the \code{property} to \samp{G.species}, \samp{Cp.species}, etc.
Except for \samp{A}, the properties of proteins or their reactions calculated in this manner are restricted to nonionized proteins.

Zero, one, or more leading arguments to the function identify which of the chemical activities of basis species, temperature, pressure and/or ionic strength to vary.
The names of each of these arguments may be the formula of any of the basis species of the system, or \samp{T}, \samp{P}, \samp{pe}, \samp{pH}, \samp{Eh}, or \samp{IS} (but names may not be repeated).
To use the names of charged basis species such as \samp{K+} and \samp{SO4-2} as the arguments, they should be enclosed in quotes (see the example for aluminum speciation in \code{\link{diagram}}).
The value of each argument is of the form \code{c(min,max)} or \code{c(min,max,res)} where \code{min} and \code{max} refer to the minimimum and maximum values of variable identified by the name of the argument, and \code{res} denotes the resolution, or number of points along which to do the calculations (this is assigned a default value of 128 if it is missing).
For any arguments that refer to basis species, the numerical values are the logarithms of the activities of that basis species, or logarithms of fugacities if it is a gas.
Unlike the \code{energy} function, the units of \samp{T} and \samp{P} in \code{affinity} are those the user has set using \code{\link{T.units}} and \code{\link{P.units} }(on program start-up these are \eqn{^{\circ}}{°}C and bar, respectively). 

If one or more buffers are assigned to the definition of \code{\link{basis}} species, \code{affinity} calls \code{\link{buffer}} to calculate the logarithms of activities of these basis species from the buffer.

The \code{iprotein} and \code{loga.protein} arguments can be used to compute the chemical affinities of formation reactions of proteins that are not in the current \code{\link{species}} definition.
This approach can be utilized in order to calculate the properties of many proteins in a fraction of the time it would take to calculate them individually.
The appropriate \code{\link{basis}} species still must be defined prior to calling \code{affinity}.
\code{iprotein} contains indices of desired proteins in \code{\link{thermo}$protein}; \code{affinity} adds to the species list the amino acid residues and and terminal H2O group (indicated by \dQuote{RESIDUE} in \code{thermo$protein}) then calculates the properties of the reactions for the residues and terminal group, including ionization effects, and adds them together to get those of the indicated proteins.

In CHNOSZ version 0.9, \code{energy} gained a new argument \samp{transect} which is set to TRUE by \code{energy.args} when the length(s) of the variables is(are) greater than three.
In this mode of operation, instead of performing the calculations on an \eqn{n}{n}-dimensional grid, the affinities are calculated on a transect of changing T, P, and/or chemical activity of basis species.
}

\value{
For \code{affinity}, a list, elements of which are \code{sout} output from \code{\link{subcrt}}, \code{property} name of the calculated property (\samp{A} for chemical affinity), \code{basis} and \code{species} definition of basis species and species of interest in effect at runtime, \code{T} and \code{P} temperature and pressure, in the system units of Kelvin and bar, set to \code{numeric()} (length=0) if either one is a variable, \code{vars} the names of the variables, \code{vals} the values of the variables (a list, one element for each variable), \code{values} the result of the calculation (a list, one element for each species, with names taken from the species index in \code{\link{thermo}$obigt}).
The elements of the lists in \code{vals} and \code{values} are arrays of \eqn{n}{n} dimensions, where \eqn{n}{n} is the number of variables.
The values of chemical affinity of formation reactions of the species are returned in dimensionless units (for use with decimal logarithms, i.e., \emph{\strong{A}}/\eqn{2.303RT}).

Names other than \samp{T} or \samp{P} in \code{vars} generally refer to basis species, and the corresponding \code{vals} are the logarithms of activity or fugacity. However, if one or more of \code{pe}, \code{Eh} or \code{pH} is among the variables of interest, \code{vals} holds the values of the those variables as indicated.
}

\seealso{
Usually, \code{\link{equilibrate}} is the next step in calculations of chemical equilibrium.
The help for \code{\link{buffer}} has some examples of using chemical activity buffers.
The guts of the calculations provided by \code{affinity} involve \code{\link{energy.args}} and \code{\link{energy}}, which normally are not part of the user interaction.
Longer examples are in \code{demo{affinity}}.
}

\examples{
\dontshow{data(thermo)}## set up a system and calculate
## chemical affinities of formation reactions
basis(c("SiO2", "MgO", "H2O", "O2"), c(-5, -5, 0, 999))
species(c("quartz","enstatite","forsterite"))
# chemical affinities (A/2.303RT) at 25 deg C and 1 bar
affinity()
# at higher temperature and pressure
affinity(T=500, P=2000)
# at 25 temperatures and pressures
affinity(T=c(500, 1000, 5), P=c(1000, 5000, 5))
# equilibrium constants of formation reactions
affinity(property="logK")
# standard molal Gibbs energies of species,
# user units (default: cal/mol)
affinity(property="G.species")
# standard molal Gibbs energies of reactions
affinity(property="G")
# a T,P-transect
# (fluid pressure from Helgeson et al., 2009 Fig. 7)
affinity(T=c(25, 110, 115, 215), P=c(11, 335, 500, 1450))
}

\references{
  Helgeson, H. C., Richard, L, McKenzie, W. F., Norton, D. L. and Schmitt, A. (2009) A chemical and thermodynamic model of oil generation in hydrocarbon source rocks. \emph{Geochim. Cosmochim. Acta} \bold{73}, 594--695. \url{http://dx.doi.org/10.1016/j.gca.2008.03.004}

  Kondepudi, D. K. and Prigogine, I. (1998) \emph{Modern Thermodynamics: From Heat Engines to Dissipative Structures}, John Wiley & Sons, New York, 486 p. \url{http://www.worldcat.org/oclc/38055900}
}


\concept{Primary thermodynamic modeling}
