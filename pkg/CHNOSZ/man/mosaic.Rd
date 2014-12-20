\encoding{UTF-8}
\name{mosaic}
\alias{mosaic}
\title{Chemical Affinities with Changing Basis Species}
\description{
Calculate chemical affinities of formation reactions of species using basis species that change with the conditions.
}

\usage{
  mosaic(bases, blend=FALSE, ...)
}

\arguments{
  \item{bases}{character, basis species to include in the calculation}
  \item{blend}{logical, use relative abundances of basis species?}
  \item{...}{additional arguments to be passed to \code{\link{affinity}}}
}

\details{

\code{mosaic} can be used to calculate the reaction affinities when the basis species listed in \code{bases} change in relative abundance over the range of conditions, due to e.g. ionization, complexation or redox reactions.
This is a way to \dQuote{speciate the basis species}.
An example is consideration of the speciation of sulfur (\samp{SO4-2}, \samp{HSO4-}, \samp{HS-} and \samp{H2S}) as a function of oxygen fugacity (or Eh) and pH in calculating the relative stabilities of minerals and aqueous species in the system Fe-S-O-H.
The first species listed in \code{bases} should correspond to one of the basis species currently defined, and the function only supports calculations using basis species that all derive from that species (i.e. all share the same element).
The arguments in \code{...} are passed to \code{affinity} to specify the conditions.

The function calculates the affinities using each basis species in turn, changing them via \code{\link{swap.basis}}.
If \code{blend} is FALSE (the default), the function returns the affinities calculated using the single predominant basis species in \code{bases} at each condition.
Chemical activity or predominance diagrams constructed using this method have been described as \dQuote{mosaic diagrams} in the literature.
If \code{blend} is TRUE, the function combines the affinities of the formation reactions weighted by the relative abundances of the basis species at each condition.
This tends to produce curved boundaries.

}

\value{
A list containing \code{A.species} (affinities of formation of the species with changing basis speceis) and \code{A.basis} (affinities of formation of the basis species in terms of the first basis species), each having same structure as the list returned by \code{\link{affinity}}.
}

\examples{
\dontshow{data(thermo)}# Fe-minerals and aqueous species in Fe-S-O-H system
# speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
# after Garrels and Christ, 1965 Figure 7.20
pH <- c(0, 14, 200)
Eh <- c(-1, 1, 200)
T <- 25
basis(c("FeO", "SO4-2", "H2O", "H+", "e-"))
basis("SO4-2", -6)
species(c("Fe+2", "Fe+3"), -6)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite"))
# the basis species we'll swap through
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
# calculate affinities using the predominant basis species
m1 <- mosaic(bases, pH=pH, Eh=Eh, T=T)
# make a diagram and add water stability lines
diagram(m1$A.species)
water.lines("pH", "Eh", T=convert(T, "K"), col="seagreen", lwd=1.5)
# show lines for Fe(aq) = 10^-4 M
species(c("Fe+2", "Fe+3"), -4)
m2 <- mosaic(bases, pH=pH, Eh=Eh, T=T)
diagram(m2$A.species, add=TRUE, names=NULL, dotted=3)
title(main=paste("Iron oxides and sulfides in water, log(total S) = -6",
  "After Garrels and Christ, 1965", sep="\n"))
# we could overlay the basis species predominance fields
#diagram(m1$A.basis, add=TRUE, col="blue", col.names="blue", dotted=3)
}

\references{
  Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}
}

\keyword{secondary}
