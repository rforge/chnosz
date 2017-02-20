\encoding{UTF-8}
\name{protein.info}
\alias{protein.formula}
\alias{protein.length}
\alias{protein.info}
\alias{protein.basis}
\alias{protein.equil}
\title{Summaries of Thermodynamic Properties of Proteins}

\description{
  Calculate chemical formulas, lengths, standard Gibbs energies and net charges, stoichiometric coefficients of basis species in reactions to form proteins (possibly per residue), and show steps in calculation of chemical activities of proteins in metastable equilibrium.
}

\usage{
  protein.info(protein, organism=NULL, residue=FALSE)
  protein.formula(protein, organism = NULL, residue = FALSE)
  protein.length(protein, organism = NULL)
  protein.basis(protein, T = 25, normalize = FALSE)
  protein.equil(protein, T=25, loga.protein = 0, digits = 4)
}

\arguments{
  \item{protein}{character, names of proteins; numeric, species index of proteins; data frame; amino acid composition of proteins}
  \item{organism}{character, names of organisms}
  \item{residue}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{normalize}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{T}{numeric, temperature in \eqn{^{\circ}}{°}C}
  \item{loga.protein}{numeric, decimal logarithms of reference activities of proteins}
  \item{digits}{integer, number of significant digits (see \code{\link{signif}})}
}

\details{
For character \code{protein}, \code{protein.info} returns the rownumber(s) of \code{thermo$protein} that match the protein names.
The names can be supplied in the single \code{protein} argument (with an underscore) or as individual \code{protein}s and \code{organism}s.
Any protein not matched returns an NA and generates a message.

For numeric \code{protein}, \code{protein.info} returns the corresponding row(s) of \code{thermo$protein}.
Set \code{residue} to TRUE to return the per-residue composition (i.e. amino acid composition of the protein divided by total number of residues).

For dataframe \code{protein}, \code{protein.info} returns it unchanged, except for possibly the per-residue calculation.

The following functions accept any specification of protein(s) described above for \code{protein.info}:

\code{protein.formula} returns a stoichiometrix matrix representing the chemical formulas of the proteins that can be pased to e.g. \code{\link{mass}} or \code{\link{ZC}}.
The amino acid compositions are multiplied by the output of \code{\link{group.formulas}} to generate the result. 

\code{protein.length} returns the lengths (number of amino acids) of the proteins.

The following functions also depend on an existing definition of the basis species:

\code{protein.basis} calculates the numbers of the basis species (i.e. opposite of the coefficients in the formation reactions) that can be combined to form the composition of each of the proteins.
The basis species must be present in \code{thermo$basis}, and if \samp{H+} is among the basis species, the ionization states of the proteins are included.
The ionization state of the protein is calculated at the pH defined in \code{thermo$basis} and at the temperature specified by the \code{T} argument.
If \code{normalize} is TRUE, the coefficients on the basis species are divided by the lengths of the proteins. 

  \code{protein.equil} produces a series of messages showing step-by-step a calculation of the chemical activities of proteins in metastable equilibrium. For the first protein, it shows the standard Gibbs energies of the reaction to form the nonionized protein from the basis species and of the ionization reaction of the protein (if \samp{H+} is in the basis), then the standard Gibbs energy/RT of the reaction to form the (possibly ionized) protein per residue. The per-residue values of \samp{logQstar} and \samp{Astar/RT} are also shown for the first protein. Equilibrium calculations are then performed, only if more than one protein is specified. This calculation applies the Boltzmann distribution to the calculation of the equilibrium degrees of formation of the residue equivalents of the proteins, then converts them to activities of proteins taking account of \code{loga.protein} and protein length. If the \code{protein} argument is numeric (indicating rownumbers in \code{thermo$protein}), the values of \samp{Astar/RT} are compared with the output of \code{\link{affinity}}, and those of the equilibrium degrees of formation of the residues and the chemical activities of the proteins with the output of \code{\link{diagram}}. If the values in any of these tests are are not \code{\link{all.equal}} an error is produced indicating a bug. 
}

\examples{\dontshow{data(thermo)}
# search by name in thermo$protein
ip1 <- protein.info("LYSC_CHICK")
ip2 <- protein.info("LYSC", "CHICK")
# these are the same
stopifnot(all.equal(ip1, ip2))
# two organisms with the same protein name
ip3 <- protein.info("MYG", c("HORSE", "PHYCA"))
# their amino acid compositions
protein.info(ip3)
# their thermodynamic properties by group additivity
aa2eos(protein.info(ip3))

# an example of an unrecognized protein name
ip4 <- protein.info("MYGPHYCA")
stopifnot(is.na(ip4))

## example for chicken lysozyme C
# index in thermo$protein
ip <- protein.info("LYSC_CHICK")
# amino acid composition
protein.info(ip)
# length and chemical formula
protein.length(ip)
protein.formula(ip)
# group additivity for thermodynamic properties and HKF equation-of-state
# parameters of non-ionized protein
aa2eos(protein.info(ip))
# calculation of standard thermodynamic properties
# (subcrt uses the species name, not ip)
subcrt("LYSC_CHICK")
# affinity calculation, protein identified by ip
basis("CHNOS+")
affinity(iprotein=ip)
# affinity calculation, protein loaded as a species
species("LYSC_CHICK")
affinity()
# NB: subcrt() only shows the properties of the non-ionized
# protein, but affinity() uses the properties of the ionized
# protein if the basis species have H+

## these are all the same
protein.formula("P53_PIG")
protein.formula(protein.info("P53_PIG"))
protein.formula(protein.info(protein.info("P53_PIG")))

## using protein.formula: average oxidation state of 
## carbon of proteins from different organisms (Dick, 2014)
# get amino acid compositions of microbial proteins 
# generated from the RefSeq database 
file <- system.file("extdata/refseq/protein_refseq.csv.xz", package="CHNOSZ")
ip <- add.protein(read.aa(file))
# only use those organisms with a certain
# number of sequenced bases
ip <- ip[as.numeric(thermo$protein$abbrv[ip]) > 50000]
pf <- protein.formula(thermo$protein[ip, ])
zc <- ZC(pf)
# the organism names we search for
# "" matches all organisms
terms <- c("Natr", "Halo", "Rhodo", "Acido", "Methylo",
  "Chloro", "Nitro", "Desulfo", "Geo", "Methano",
  "Thermo", "Pyro", "Sulfo", "Buchner", "")
tps <- thermo$protein$ref[ip]
plot(0, 0, xlim=c(1, 15), ylim=c(-0.3, -0.05), pch="",
  ylab=expression(italic(Z)[C]),
  xlab="", xaxt="n", mar=c(6, 3, 1, 1))
for(i in 1:length(terms)) {
  it <- grep(terms[i], tps)
  zct <- zc[it]
  points(jitter(rep(i, length(zct))), zct, pch=20)
}
terms[15] <- paste("all", length(ip))
axis(1, 1:15, terms, las=2)
title(main=paste("Average oxidation state of carbon in proteins",
  "by taxID in NCBI RefSeq (after Dick, 2014)", sep="\n"))
}

\references{
  Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \url{http://dx.doi.org/10.1098/rsif.2013.1095}
}

\concept{Protein thermodynamic modeling}
