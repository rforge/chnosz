\encoding{UTF-8}
\name{yeast.aa}
\alias{yeast.aa}
\title{Proteins from Model Organisms}
\description{
  Retrieve the amino acid compositions of one or more proteins from \emph{Saccharomyces cerevisiae}.
}

\usage{
  yeast.aa(protein = NULL)
}

\arguments{
  \item{protein}{character, name of protein}
}

\details{
\code{yeast.aa} retrieves the amino acid composition(s) of the indicated proteins in \emph{Saccharomyces cerevisiae}.
The calculation depends on the data file \code{\link{extdata}/protein/Sce.csv.xz}, which contains the amino acid compositions of the proteins.
The \code{protein} argument should be a vector or a list of vectors of one or more SGD IDs, Open Reading Frame (ORF) or gene names that are found in these files.
The output data frame contains rows with NA compositions for names that are not matched.
}

\value{
A data frame, or list of data frames, containing the amino acid composition(s) of the specified protein(s) in the format of \code{\link{thermo}$protein}.
}

\seealso{
\code{\link{extdata}} describes the sources of compositional data for the proteins. Other examples of usage of \code{yeast.aa} are shown for \code{\link{read.expr}}.
}


\examples{
\dontshow{data(thermo)}
# the first few names in UniProt for "aminotransferase yeast"
genes <- c("AATC", "ARO8", "BCA1", "AMPL", "BCA2", "ARO9")
# the corresponding ORF names
ORF <- c("YLR027C", "YGL202W", "YHR208W", "YKL103C", "YJR148W", "YHR137W")
# we only match two of them by gene name, but all by ORF name
aa <- yeast.aa(genes)
aa <- yeast.aa(ORF)
# what are their formulas and average oxidation states of carbon
protein.formula(aa)
ZC(protein.formula(aa))
}

\concept{Protein thermodynamic modeling}
