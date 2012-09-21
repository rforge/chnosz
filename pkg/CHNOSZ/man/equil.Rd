\name{equil}
\alias{equil.boltzmann}
\alias{equil.reaction}
\title{Equilibrium Logarithms of Activity of Species}
\description{
Calculate equilibrium chemical activities of species from the affinities of formation of the species at unit activity. 
}

\usage{
  equil.boltzmann(Astar, n.balance, loga.balance)
  equil.reaction(Astar, n.balance, loga.balance)
}

\arguments{
  \item{Astar}{numeric, affinities of formation reactions excluding species contribution}
  \item{n.balance}{numeric, number of moles of conserved component in the formation reactions of the species of interest}
  \item{loga.balance}{numeric, logarithm of total activity of balanced quantity}
}

\details{
  The equilibrium chemical activities of species shown by \code{\link{diagram}} are calculated using either the \code{equil.boltzmann} or \code{equil.reaction} function. The former is used if \code{residue} in \code{diagram} is set to \code{TRUE}, and the latter if \code{residue} is set to \code{FALSE}. 

  The input values are in a list, \code{Astar}, all elements of the list having the same dimensions; they can be vectors, matrices, or higher-dimensionsal arrays. Each list element contains the chemical affinities of the formation reactions of one of the species of interest (in dimensionless base-10 units, i.e. A/2.303RT), calculated at unit activity of the species of interest. The equilibrium activities (in base-10 log units) returned by either function satisfy the constraints that 1) the resulting chemical affinities of the formation reactions of the species are all equal and 2) the total activity of the conserved component, given by summing the products of the activities of the species by the respective coefficients in \code{n.balance}, is equal to (\code{loga.balance}). 

  The first constraint above does not necessarily impose a complete equilibrium, where the affinities of the formation reactions are all equal to zero, but also allows for a metastable equilibrium, where the affinities of the formation reactions are negative (and all equal).

  In \code{equil.reaction} (the algorithm described in Dick, 2008 and the only one available prior to CHNOSZ-0.8), the calculations of relative abundances of species are based on a solving a system of equations representing the two constraints stated above. A close approximation of the solution is found using \code{\link{uniroot}}. Prior to CHNOSZ_0.9-8, the values in the \code{Astar} were used to generate initial guesses of the logarithms of activities of species; values of \code{loga.balance} that were hugely different from these guesses could generate errors from \code{uniroot} such as "f() values at end points not of opposite sign". Now (from version 0.9-8), a more robust method for generating guesses is in place. 

  In \code{equil.boltzmann} (algorithm available beginning with CHNOSZ-0.8), the chemical activities of species are calculated using the Boltzmann distribution. This calculation is faster than the equation-solving approach used above, but is limited to systems where the coefficients of the conserved components (\code{n.balance}) are all unity. Since the solution is found in a forward manner, it is more exact than the result of \code{equil.reaction} for the same problem; in most practical settings their results are probably indistinguishable. This function is only called by \code{diagram} when \code{residue} is set to TRUE, usually for proteins.

  Both functions can make use of parallelization for intensive parts of the calculations (see \code{\link{palply}}).
}

\seealso{
  The \code{logact} component of the return value of \code{\link{diagram}} is calculated using these functions.
}

\value{
  A list with dimensions and length equal to \code{Astar}, giving the \code{\link{log10}} of the equilibrium activities of the species of interest.
}

\examples{
\dontshow{data(thermo)}
## equilibrium in a simple system:
## ionization of carbonic acid
basis("CHNOS+")
species(c("CO2", "HCO3-", "CO3-2"))
# set unit activity of the species (0 = log10(1))
species(1:3, 0)
# calculate Astar (for unit activity)
res <- 100
Astar <- affinity(pH=c(0, 14, res))$values
# the logarithms of activity for a total activity
# of the balanced quantity (C or CO2) equal to 0.001
loga.boltz <- equil.boltzmann(Astar, c(1, 1, 1), 0.001)
# calculated another way
loga.react <- equil.reaction(Astar, c(1, 1, 1), 0.001)
# probably close enough for most purposes
stopifnot(all.equal(loga.boltz, loga.react))
# the first ionization constant (pKa)
ipKa <- which.min(abs(loga.boltz[[1]] - loga.boltz[[2]]))
pKa.equil <- seq(0, 14, length.out=res)[ipKa]
# calculate logK directly
logK <- subcrt(c("CO2","H2O","HCO3-","H+"), c(-1, -1, 1, 1), T=25)$out$logK
# we could decrease tolerance here by increasing res
stopifnot(all.equal(pKa.equil, -logK, tol=1e-2))
}

\references{
  Dick, J. M. (2008) Calculation of the relative metastabilities of proteins using the CHNOSZ software package. \emph{Geochem. Trans.} \bold{9}:10. \url{http://dx.doi.org/10.1186/1467-4866-9-10}
}

\keyword{secondary}
