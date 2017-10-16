\encoding{UTF-8}
\name{examples}
\alias{examples}
\alias{demos}
\title{Run Examples from the Documentation}
\description{
  Run the examples contained in each of the documentation topics. }

\arguments{
  \item{do.png}{logical, generate PNG files for the plots?}
  \item{which}{character, which example to run}
  \item{to.file}{logical, generate image files for the plots?}
}

\usage{
  examples(do.png = FALSE)
  demos(which = c("sources", "protein.equil", "affinity", "NaCl",
    "density", "ORP", "revisit", "findit", "ionize", "buffer",
    "protbuff", "yeastgfp", "mosaic", "copper", "solubility",
    "wjd", "dehydration", "bugstab", "Shh", "activity_ratios",
    "adenine", "DEW", "lambda", "TCA"),
    to.file=FALSE)
}

\details{

\code{examples} runs all the examples in the documentation for the package. 
\code{\link{example}} is called for each topic with \code{ask} set to \code{FALSE} (so all of the figures are shown without prompting the user).
If \code{do.png} is TRUE, the plots in the examples are saved as \code{\link{png}} files having names beginning with the name of each of the help topics.

\code{demos} is a function to run other examples that are provided as \code{\link{demo}s}.
\code{\link{demo}} is called with settings to not echo the source code and to not ask before making each plot.
The demo(s) to run is/are specified by \code{which}; the default is to run them in the order of the list below.
Demos that are displayed on the CHNOSZ website (\url{http://chnosz.net/demos}) are indicated with an asterisk.
See the comments in the source code for more information about each demo.

  \tabular{ll}{
    \code{sources} \tab Cross-check the reference list with the thermodynamic database \cr
    \code{protein.equil} \tab Chemical activities of two proteins in metastable equilibrium (Dick and Shock, 2011) \cr
    \code{affinity} \tab Affinities of metabolic reactions and amino acid synthesis (Amend and Shock, 1998, 2001) \cr
    \code{NaCl} \tab * Equilibrium constant for aqueous NaCl dissociation (Shock et al., 1992) \cr
    \code{density} \tab * Density of \H2O, inverted from IAPWS-95 equations (\code{\link{rho.IAPWS95}}) \cr
    \code{ORP} \tab * Temperature dependence of oxidation-reduction potential for redox standards \cr
    \code{revisit} \tab * Coefficient of variation of metastable equilibrium activities of proteins \cr
    \code{findit} \tab * Minimize the standard deviation of logarithms of activities of sulfur species \cr
    \code{ionize} \tab ionize.aa(): contour plots of net charge and ionization properties of LYSC_CHICK \cr
    \code{buffer} \tab  * Minerals and aqueous species as buffers of hydrogen fugacity (Schulte and Shock, 1995) \cr
    \code{protbuff} \tab Chemical activities buffered by thiol peroxidases or sigma factors \cr
    \code{yeastgfp} \tab * Subcellular locations: log fO2 - log aH2O and log a - log fO2 diagrams (Dick, 2009) \cr
    \code{mosaic} \tab * Eh-pH diagram with two sets of changing basis species (Garrels and Christ, 1965) \cr
    \code{copper} \tab * Another example of \code{\link{mosaic}}: complexation of Cu with glycine (Aksu and Doyle, 2001) \cr
    \code{solubility} \tab * Solubility of calcite (cf. Manning et al., 2013) or \CO2 (cf. Stumm and Morgan, 1996) \cr
    \code{wjd} \tab * Gibbs energy minimization: prebiological atmospheres (Dayhoff et al., 1964) and cell periphery of yeast \cr
    \code{dehydration} \tab * \logK of dehydration reactions; SVG file contains tooltips and links \cr
    \code{bugstab} \tab * Formation potential of microbial proteins in colorectal cancer (Dick, 2016) \cr
    \code{Shh} \tab * Affinities of transcription factors relative to Sonic hedgehog (Dick, 2015) \cr
    \code{activity_ratios} \tab * Mineral stability plots with activity ratios on the axes \cr
    \code{adenine} \tab * HKF parameters regressed from heat capacity and volume of aqueous adenine (Lowe et al., 2017) \cr
    \code{DEW} \tab * Deep Earth Water (DEW) model for high pressures (Sverjensky et al., 2014a and 2014b) \cr
    \code{lambda} \tab * Effects of lambda transition on thermodynamic properties of quartz (Berman, 1988) \cr
    \code{TCA} \tab * Standard Gibbs energies of the tricarboxylic (citric) acid cycle (Canovas and Shock, 2016) \cr
  }

}

\references{
Aksu, S. and Doyle, F. M. (2001) Electrochemistry of copper in aqueous glycine solutions. \emph{J. Electrochem. Soc.} \bold{148}, B51--B57. \url{https://doi.org/10.1149/1.1344532}

Amend, J. P. and Shock, E. L. (1998) Energetics of amino acid synthesis in hydrothermal ecosystems. \emph{Science} \bold{281}, 1659--1662. \url{https://doi.org/10.1126/science.281.5383.1659}

Amend, J. P. and Shock, E. L. (2001) Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. \emph{FEMS Microbiol. Rev.} \bold{25}, 175--243. \url{https://doi.org/10.1016/S0168-6445(00)00062-0}

Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O–K{\s2}O–CaO–MgO–FeO–Fe{\s2}O{\s3}–Al{\s2}O{\s3}–SiO{\s2}–TiO{\s2}–H{\s2}O–CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \url{https://doi.org/10.1093/petrology/29.2.445}

Canovas, P. A., III and Shock, E. L. (2016) Geobiochemistry of metabolism: Standard state thermodynamic properties of the citric acid cycle. \emph{Geochim. Cosmochim. Acta} \bold{195}, 293--322. \url{https://doi.org/10.1016/j.gca.2016.08.028}

Dayhoff, M. O. and Lippincott, E. R. and Eck, R. V. (1964) Thermodynamic Equilibria In Prebiological Atmospheres. \emph{Science} \bold{146}, 1461--1464. \url{https://doi.org/10.1126/science.146.3650.1461}

Dick, J. M. (2009) Calculation of the relative metastabilities of proteins in subcellular compartments of \emph{Saccharomyces cerevisiae}. \emph{BMC Syst. Biol.} \bold{3}:75. \url{https://doi.org/10.1186/1752-0509-3-75}

Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLoS ONE} \bold{6}, e22782. \url{https://doi.org/10.1371/journal.pone.0022782}

Dick, J. M. (2015) Chemical integration of proteins in signaling and development. \emph{bioRxiv}. \url{https://doi.org/10.1101/015826}

Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}:e2238. \url{https://doi.org/10.7717/peerj.2238}

Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}

Lowe, A. R., Cox, J. S. and Tremaine, P. R. (2017) Thermodynamics of aqueous adenine: Standard partial molar volumes and heat capacities of adenine, adeninium chloride, and sodium adeninate from \emph{T} = 278.15 K to 393.15 K. \emph{J. Chem. Thermodyn.} \bold{112}, 129--145. \url{https://doi.org/10.1016/j.jct.2017.04.005}
  
Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Schulte, M. D. and Shock, E. L. (1995) Thermodynamics of Strecker synthesis in hydrothermal systems. \emph{Orig. Life Evol. Biosph.} \bold{25}, 161--173. \url{https://doi.org/10.1007/BF01581580}

Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \degC and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{https://doi.org/10.1039/FT9928800803}

Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}

Sverjensky, D. A., Harrison, B. and Azzolini, D. (2014a) Water in the deep Earth: The dielectric constant and the solubilities of quartz and corundum to 60 kb and 1,200 \degC. \emph{Geochim. Cosmochim. Acta} \bold{129}, 125--145. \url{https://doi.org/10.1016/j.gca.2013.12.019}

Sverjensky, D. A., Stagno, V. and Huang, F. (2014b) Important role for organic carbon in subduction-zone fluids in the deep carbon cycle. \emph{Nat. Geosci.} \bold{7}, 909--913. \url{https://doi.org/10.1038/ngeo2291}
}


\examples{
\dontshow{data(thermo)}
\dontshow{opar <- par(no.readonly=TRUE)}
demos(c("ORP", "NaCl"))
\dontshow{par(opar)}
\dontrun{
# use the following to run examples in all help topics
examples()}
}


\concept{Extra thermodynamic modeling}
