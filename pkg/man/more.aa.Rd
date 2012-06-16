\name{more.aa}
\alias{more.aa}
\title{Proteins from Model Organisms}
\description{
  Retrieve the amino acid compositions of one or more proteins from \emph{Escherichia coli} or \emph{Saccharomyces cerevisiae}.
}

\usage{
  more.aa(protein = NULL, organism)
}

\arguments{
  \item{protein}{character, name of protein}
  \item{organism}{character, name of organism (\samp{Eco} or \samp{Sgd})}
}

\details{
  \code{more.aa} retrieves the amino acid composition(s) of the proteins for either \emph{Escherichia coli} or \emph{Saccharomyces cerevisiae} (for \code{organism} equal to \samp{ECO} or \samp{SGD}, respectively). The calculation depends on the data files \code{\link{extdata}/protein/Eco.csv.xz} and \code{\link{extdata}/protein/Sce.csv.xz}, which contain the amino acid compositions of proteins in these organisms. The \code{protein} can one or more Ordered Locus Names (OLN) or Open Reading Frame (ORF) names that are found in these files. The output data frame contains rows with NA compositions for names that are not matched.
}

\value{
  A data frame containing the amino acid composition(s) of the specified protein(s) in the format of \code{\link{thermo}$protein}.
}

\seealso{
  There are examples of using this function in conjunction with lists of protein names from \code{\link{yeastgfp}} and \code{\link{read.expr}}.
}


\examples{
\dontshow{data(thermo)}
# the first 13 names in UniProt for "aminotransferase ecoli"
at.ecoli <- c("BIOA", "ARNB", "SERC", "AAT", "TYRB", "ARGD", 
  "ILVE", "ALAA", "ALAC", "YBDL", "AVTA", "GLMS", "PUUE")
# get the amino acid compositions
# note that ALAA and ALAC are not matched
at.aa <- more.aa(at.ecoli, "Eco")
# what are their formulas and oxidation states
protein.formula(at.aa)
ZC(protein.formula(at.aa))
}

\keyword{protein}
