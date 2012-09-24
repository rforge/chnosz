\name{revisit}
\alias{revisit}
\alias{richness}
\alias{extremes}
\alias{where.extreme}
\title{Summarize Logarithms of Activity of Chemical Species}
\description{
  For logarithms of activities of chemical species, calculate standard deviation, coefficient of variation, Shannon diversity index, or correlation coefficient with a normal distribution. For logarithms of activities of species compared with reference values, calculate Spearman or Pearson correlation coefficients, or coefficient of variation of the root mean squared deviation. For logarithms of activities compared with one reference value, return the logarithm of activity (of a single species, unchanged) or richness (number of species with activities greater than the reference value). For differences of logarithms of activities of species from reference values, compared with reference differences, calculate Spearman or Pearson correlation coefficients, or coefficient of variation of the root mean squared deviation. Report the position and of the minimum or maximum of the values and make line or contour plots showing the values.
}

\usage{
  revisit(d, target = "cv", loga.ref = NULL,
    plot.it = NULL, col = par("fg"), yline = 2, ylim = NULL, 
    ispecies = NULL, add = FALSE, cex = par("cex"), lwd = par("lwd"), 
    mar = NULL, side = 1:4, xlim = NULL, labcex = 0.6, pch = 1, 
    legend = "", legend.x = NULL, lpch = NULL, main = NULL, 
    lograt.ref = NULL, plot.ext = TRUE, DGxf.swap12 = FALSE)
  extremes(z, target)
  where.extreme(z, target, do.sat = FALSE)
}

\arguments{
  \item{d}{list, output from \code{\link{diagram}}, or list of logarithms of activities of species}
  \item{target}{character, what statistic to calculate}
  \item{loga.ref}{numeric, logarithm of activities for comparison statistics}
  \item{plot.it}{logical, make a plot?}
  \item{col}{character, color to use for points or lines}
  \item{yline}{numeric, margin line for y-axis label}
  \item{ylim}{numeric, limits of y axis}
  \item{ispecies}{numeric, which species to consider}
  \item{add}{logical, add to an existing plot?}
  \item{cex}{numeric, character expansion factor}
  \item{lwd}{numeric, line width}
  \item{mar}{numeric, plot margin specifications}
  \item{side}{numeric, which sides of plot to draw axes}
  \item{xlim}{numeric, limits of x axis}
  \item{labcex}{numeric, character expansion factor for species labels}
  \item{pch}{numeric, plotting symbol(s) to use for points}
  \item{legend}{character, text to use for legend}
  \item{legend.x}{character, placement of legend}
  \item{lpch}{numeric, plotting symbol(s) to use in legend}
  \item{main}{character, main title for plot}
  \item{lograt.ref}{numeric, log10 of reference abundance ratios}
  \item{plot.ext}{logical, show the location of the extreme value(s)?}
  \item{z}{numeric, matrix of values}
  \item{do.sat}{logical, identify multiple extreme values? }
  \item{DGxf.swap12}{logical, swap order of initial and final states?}
}

\details{

  Supply the output of \code{\link{diagram}} as the value for \code{d}. The type of calculation is indicated by \code{target}, as summarized in the following table.

  \tabular{llll}{
    target \tab description \tab extremum \tab additional arguments \cr
    \samp{sd} \tab standard deviation \tab min \tab none \cr
    \samp{cv} \tab coefficient of variation \tab min \tab none \cr
    \samp{shannon} \tab Shannon diversity index \tab max \tab none \cr
    \samp{qqr} \tab correlation coefficient on q-q plot (normal distribution) \tab max \tab none \cr
    \samp{cvrmsd} \tab coefficient of variation of RMSD \tab min \tab \code{loga.ref} \cr
    \samp{spearman} \tab Spearman correlation coefficient \tab max \tab \code{loga.ref} \cr
    \samp{pearson} \tab Pearson correlation coefficient \tab max \tab \code{loga.ref} \cr
    \samp{count} \tab count of species above a certain activity \tab max \tab \code{loga.ref} \cr
    \samp{logact} \tab logarithm of activity of a single species \tab max \tab \code{loga.ref} \cr
    \samp{DGxf} \tab Gibbs energy of transformation of the system \tab min \tab {loga.ref} \cr
  }

  \code{sd}, \code{cv}, \code{shannon} and \code{qqr} all operate on the logarithms of activities of species. 

  \code{cvrmsd}, \code{spearman} and \code{pearson} operate on the logarithms of activities of species in \code{d} compared against reference values in \code{loga.ref}, which must have the same length as the number of species. However, if code \code{lograt.ref} is supplied, the ratio of the activities in \code{d} to those in \code{loga.ref} (the difference in the logarithms) is calculated instead; and those values are compared to \code{lograt.ref} using the selected target.

  \code{richness} counts the numbers of species whose logarithms of activities are above the single value of \code{loga.ref}. \code{logact} simply operates on the logarithm of activity of a single species (as supplied in \code{d}); \code{loga.ref} in this case is an integer value referring to the number of the species in the list in \code{d}. 

  \code{DGxf} is the change in Gibbs energy of a system undergoing transformation among species of interest at constant T, P and chemical potentials of basis species. See \code{\link{DGxf}} for the equation that is used. The transformation is from an initial state (the equilibrium values in \code{d}) to a final state (the values given in \code{loga.ref}); setting \code{DGxf.swap12} reverses the order, so that the equilibrium state is the final state. Unlike the other targets, \code{DGxf} also depends on the values of the starved affinity (\code{Astar}) that are returned by \code{\link{diagram}}.

  Toether with the results of the calculation of the \code{target}, the function returns the position (in as many dimensions as the input values) and value of the extremum of the target value -- either the minimum or maximum, as indicated in the table.

  If \code{plot.it} is TRUE, and \code{d} is the output from \code{diagram}, and the number of variables is 1 or 2, the results are plotted --- a line diagram in 1 dimension or a contour plot in 2 dimensions. No plot is made if the logarithms of activity have no dimensions (a single value for each species), except for the case of \code{qqr}, showing a quantile-quantile plot (\code{\link{qqnorm}}). On the plots the location of the extreme value is indicated by a dashed vertical line on a 1-D plot or a point marked by an asterisk on a 2-D plot. On 2-D plots the valleys (or ridges) leading to the location of the extremum are plotted. The ridges or valleys are plotted as dashed lines and are colored green for the \eqn{x}{x} values returned by \code{extremes} and blue for the \eqn{y}{y} values returned by \code{extremes}.

  An alternative source for the \code{d} argument is any list of numeric values (such as logarithms of activities of species); the list elements each correspond to a different observation (such as a single species) and can be vectors, matrices, or higher-dimensional arrays, but they must all have the same dimensions. However, plotting is disabled if \code{d} is not the output from \code{diagram}, since the names of the variables are not available.
 
  The function name was changed to \code{revisit} in CHNOSZ-0.9 because there was already a function named \code{\link[vegan]{diversity}} in the \pkg{vegan} package. While \code{diversity} operates on a matrix with (biological) species on the columns, \code{revisit} operates on a list with (chemical) species as the elements of the list. 

Given a matrix of numeric values in \code{z}, \code{extremes} locates the maximum or minimum values in both dimensions. That is, the \eqn{x}{x} values that are returned are the column numbers where the extreme is found for each row, and the \eqn{y}{y} values that are returned are the row numbers where the extreme is found for each column. The location of the extreme value in a matrix or vector \code{z} is calculated using \code{where.extreme}. Whether the extreme is the minimum or the maximum value depends on the value of \code{target}. For matrices, if \code{do.sat} is TRUE, and the extreme value is repeated, the row and columns numbers for all instances are returned. 

}

\value{
  A list containing at least an element named \samp{H} giving the calculated values for the target. This has the same dimensions as a single element of \code{d$logact} (or \code{d} itself, if \code{d} was not the output from \code{diagram}). For calculations as a function of one or two variables, the output also contains the elements \code{ix} (location of the extremum in the first direction), \code{x} (x-value at the extremum), and \code{extval} (extreme value). For calculations as a function of two variables, the output also contains the elements \code{iy} (location of the extremum in the second direction) and \code{y} (y-value at the extremum).
}

\seealso{
  \code{\link{findit}} calls \code{revisit} iteratively to optimize one or more chemical activities, temperature and/or pressure for one of the target statistics.
}

\examples{
\donttest{
## 'revisit' calculations for amino acids
opar <- par(mfrow=c(2, 2))
basis("CHNOS+")
species(aminoacids(""))
# chemical affinities as a function of logarithm of oxygen fugacity
a <- affinity(O2=c(-85, -60))
# shows the equilibrium abundances of the amino acids
d <- diagram(a)
mtitle(c("20 amino acids", "balanced on CO2"))
# show a legend with input constraints
db <- describe.basis(ibasis=3)
dp <- describe.property("T", 25)
legend("bottomright", c(dp, db))
# default is to plot coefficient of variation
r <- revisit(d)
# show a title with the optimal conditions
mincv <- format(r$extval, digits=3)
t1 <- paste("minimum coeff of variation,", mincv, "at:")
# the logfO2 that minimized the C.V.
basis("O2", r$x)
t2 <- describe.basis(ibasis=5)
mtitle(c(t1, t2))
# chemical affinities as a function of two other variables
a <- affinity(NH3=c(-10, 10, 40), T=c(0, 80, 40))
d <- diagram(a, mam=FALSE)
# show a legend with input constraints
db <- describe.basis(ibasis=5)
legend("bottomright", db)
r <- revisit(d)
# show a title with the optimal conditions
mincv <- format(r$extval, digits=3)
t1 <- paste("minimum coeff of variation,", mincv, "at:")
# the logaNH3 and T that minimized the C.V.
basis("NH3", r$x)
db <- describe.basis(ibasis=3)
dp <- describe.property("T", r$y)
t2 <- substitute(list(dp, db), list(dp=dp[[1]], db=db[[1]]))
mtitle(c(t1, as.expression(t2)))
par(opar)

## calculations for proteins in Pelagibacter ubique
## using grep.file, read.fasta, add.protein
f <- system.file("extdata/fasta/HTCC1062.faa.xz", package="CHNOSZ")
# what proteins to select (set to "" for all proteins)
w <- "ribosomal"
# locate entries whose names contain w
j <- grep.file(f, w)
# get the amino acid compositions of these protein
aa <- read.fasta(f, j)
# add these proteins to CHNOSZ's inventory
ip <- add.protein(aa)
# set up a the chemical system
basis("CHNOS+")
# calculate affinities of formation in logfO2 space
a <- affinity(O2=c(-85, -60), iprotein=ip)
# show the equilibrium activities
opar <- par(mfrow=c(2, 2))
d <- diagram(a, loga.balance=0)
# make a title
expr <- as.expression(substitute(x~y~"proteins in"~
  italic("P. ubique"), list(x=length(j), y=w)))
mtitle(c("Equilibrium activities of", expr))
# show the coefficient of variation
revisit(d, "CV")
mtitle(c("CV of equilibrium activities of", expr))
# calculate affinities in logfO2-logaH2O space
a <- affinity(O2=c(-85,-65), H2O=c(-10,0), iprotein=ip)
# calculate the equilibrium activities
d <- diagram(a, mam=FALSE, loga.balance=0)
# show the coefficient of variation
r <- revisit(d, "CV")
stopifnot(r$ix==37)
stopifnot(r$iy==53)
mtitle(c("CV of equilibrium activities of", expr))
par(opar)
}
}

\keyword{primary}
