\name{findit}
\alias{findit}
\alias{plot.findit}
\title{Multidimensionally Optimize Chemical Activities}
\description{
  Use a gridded search to find a combination of one or more of chemical activities of basis species, temperature and/or pressure that maximize or minimize a target function of the metastable equilibrium chemical activities of the species of interest.
}

\usage{
  findit (lims = list(), target = "cv", niter = NULL, iprotein = NULL, 
    plot.it = TRUE, T = 25, P = "Psat", res = NULL, labcex = 0.6, 
    loga.ref = NULL, loga.balance = 0, rat = NULL,
    balance = NULL)
  \method{plot}{findit}(x, which=NULL, mar=c(3.5,5,2,2), xlab="iteration", \dots)
}

\arguments{
  \item{lims}{list, specification of search limits}
  \item{target}{character, target statistic to optimize}
  \item{niter}{numeric, number of iterations}
  \item{res}{numeric, grid resolution (number of points on one edge)}
  \item{iprotein}{numeric, indices of proteins}
  \item{plot.it}{logical, make a plot?}
  \item{T}{numeric, temperature}
  \item{P}{numeric, pressure; or character, "Psat"}
  \item{labcex}{numeric, character expansion for plot labels}
  \item{loga.ref}{numeric, reference logarithms of activity of species}
  \item{loga.balance}{numeric, logarithm of total activity of balanced quantity (passed to \code{\link{diagram}})}
  \item{rat}{numeric, ratio of edge length in successive iterations}
  \item{balance}{character or numeric, balanced quantity (passed to \code{\link{diagram}})}
  \item{x}{list, object of class \code{findit}}
  \item{which}{numeric, which of the parameters to plot}
  \item{mar}{numeric, plot margin specification}
  \item{xlab}{character, x-axis label}
  \item{\dots}{additional arguments passed to \code{plot}}
}

\details{
  \code{findit} implements a gridded optimization to find the minimum or maximum value of \code{target} as a function of one or more of the chemical activities, temperature and/or pressure whose ranges are listed in \code{lims}. Any target available in \code{\link{revisit}} can be optimized. Generally, the system (\code{\link{basis}} species and \code{\link{species}} of interest) must be set up before calling this function. If \code{iprotein} is supplied, indicating a set of proteins to use in the calculation, the definition of the \code{species} is not required. 

  \code{lims} is a list, each element of which contains two values indicating the range of the parameter indicated by the \code{\link{names}} of \code{lims}. The \code{names} should be formula(s) of one or more of the basis speices, \samp{T} and/or \samp{P}. If the latter two are missing, the calculations are performed at isothermal and/or isobaric conditions indicated by \code{T} and \code{P}. 

  It can be used for one, two, or more parameters. If \eqn{nd} is the number of parameters (dimensions), default values of \code{n} and \code{res} come from the following table. These settings are selected mostly for quickness of testing the function in high-dimensional space. Detailed studies of a system might have to use more iterations and/or higher resolutions.

  \tabular{rrrrr}{
    \code{nd} \tab \code{n} \tab \code{res} \tab \code{res^nd} \tab \code{rat} \cr
    1 \tab 4 \tab 128 \tab 128 \tab 0.7 \cr
    2 \tab 6 \tab 64 \tab 4096 \tab 0.7 \cr
    3 \tab 6 \tab 16 \tab 4096 \tab 0.8 \cr
    4 \tab 8 \tab 8 \tab  4096 \tab 0.9 \cr
    5 \tab 12 \tab 6 \tab 7776 \tab 0.9 \cr
    6 \tab 12 \tab 4 \tab 4096 \tab 0.95 \cr
    7 \tab 12 \tab 4 \tab 16384 \tab 0.95 \cr
  }

  The function performs \code{niter} iterations. At first, the limits of the parameters given in \code{lims} define the extent of a \eqn{nd}-dimensional box around the space of interest. The value of \code{target} is calculated at each of the \eqn{res^{nd}}{res^nd} grid points and and optimum value located (see \code{\link{revisit}} and \code{\link{where.extreme}}). In the next iteration the new search box is centered on the location of the optimum value, and the edges are shrunk so their length is \code{rat} * the length in the previous step. If the limits of any of the parameters extend beyond those in \code{lims}, they are pushed in to fit (preserving the difference between them).

  \code{plot.findit} plots the values of the parameters and the target statistic as a function of the number of iterations.

}

\value{
  \code{findit} returns a list having class \code{findit} with elements \code{value} (values of the parameters, and value of the target statistic, at each iteration), \code{lolim} (lower limits of the parameters) and \code{hilim} (upper limits of the parameters).
}

\seealso{
  \code{\link{longex}("findit")} and \code{tests/test-findit.R} for more examples.
}

\examples{\dontshow{data(thermo)}
\donttest{## an inorganic example: sulfur species
basis("CHNOS+")
basis("pH", 5)
species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2", "S3O6-2",
  "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
# to minimize the standard deviations of the loga of the species
target <- "sd"
# in one dimension, it's close to logfO2=-27.8
f1 <- findit(list(O2=c(-50, -15)), target, T=325, P=350, n=3)
title("S.D. of equilibrium log activities of sulfur species")
# in two dimensions, it's close to logfO2=-30.0 pH=9.38
f2 <- findit(list(O2=c(-50, -15), pH=c(0, 14)), target, T=325, P=350, res=16, n=4)
title("S.D. of equilibrium log activities of sulfur species")
}
}

\keyword{primary}
