\encoding{UTF-8}
\name{thermo}
\alias{reset}
\alias{obigt}
\alias{thermo}
\alias{CHNOSZ}
\title{Thermodynamic Database and System Settings}
\description{
Run \code{reset()} to reset all of the data used in CHNOSZ to default values.
This includes the computational settings, thermodynamic database, and system settings (chemical species).

The system settings are changed using \code{\link{basis}} and \code{\link{species}}.
To clear the system settings (the default, i.e. no species loaded), run \code{basis("")}; to clear only the formed species, run \code{species(delete = TRUE)}

The thermodynamic database is changed using \code{\link{add.obigt}} and \code{\link{mod.obigt}}.
To restore the default database without altering the species settings, run \code{obigt()}.

The computational settings are changed using \code{\link{water}}, \code{\link{P.units}}, \code{\link{T.units}}, \code{\link{E.units}}, and some other commands (e.g. \code{\link{mod.buffer}}).

All the data are stored in the \code{thermo} data object in an environment named \code{CHNOSZ}.
\code{thermo()} is a convenience function to access or modify parts of this object, in particular some computational settings, for example, \code{thermo("opt$ideal.H" = FALSE)} (see \code{\link{nonideal}}).

The main data files provided with CHNOSZ, as \code{*.csv} files in the \code{extdata/thermo} and \code{extdata/OBIGT} directories of the package, are used to build the \code{thermo} object, which is described below.
}

\usage{
  reset()
  obigt()
  thermo(...)
}

\arguments{
  \item{...}{list, one or more arguments whose names correspond to the component() to modify}
}

\format{

  \itemize{
     
    \item \code{thermo$opt} 
    List of computational settings. Square brackets indicate default values.
    \tabular{lll}{
      \code{cutoff} \tab numeric \tab Cutoff below which values are taken to be zero [1e-10] (see \code{\link{makeup}})\cr
      \code{E.units} \tab character \tab The user's units of energy ([\samp{cal}] or \samp{J}) (see \code{\link{subcrt}})\cr
      \code{T.units} \tab character \tab The user's units of temperature ([\samp{C}] or \samp{K})\cr
      \code{P.units} \tab character \tab The user's units of pressure ([\samp{bar}] or \samp{MPa})\cr
      \code{state} \tab character \tab The default physical state for searching species [\samp{aq}] (see \code{\link{info}})\cr
      \code{water} \tab character \tab Computational option for properties of water ([\samp{SUPCRT}] or \samp{IAPWS}; see \code{\link{water}})\cr
      \code{G.tol} \tab numeric \tab Difference in G above which \code{\link{checkGHS}} produces a message (cal mol\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}}) [100]\cr
      \code{Cp.tol} \tab numeric \tab Difference in Cp above which \code{\link{checkEOS}} produces a message (cal K\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}} mol\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}}) [1]\cr
      \code{V.tol} \tab numeric \tab Difference in V above which \code{\link{checkEOS}} produces a message (cm\ifelse{latex}{\eqn{^{3}}}{\ifelse{html}{\out{<sup>3</sup>}}{^3}} mol\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}}) [1]\cr
      \code{varP} \tab logical \tab Use variable-pressure standard state for gases? [\code{FALSE}] (see \code{\link{subcrt}})\cr
      \code{IAPWS.sat} \tab character \tab State of water for saturation properties [\samp{liquid}] (see \code{\link{util.water}})\cr
      \code{paramin} \tab integer \tab Minimum number of calculations to launch parallel processes [1000] (see \code{\link{palply}}) \cr
      \code{ideal.H} \tab logical \tab Should \code{\link{nonideal}} ignore the proton? [\code{TRUE}] \cr
      \code{ideal.e} \tab logical \tab Should \code{\link{nonideal}} ignore the electron? [\code{TRUE}] \cr
      \code{nonideal} \tab character \tab Option for charged species in \code{\link{nonideal}} [\code{Bdot}] \cr
      \code{Setchenow} \tab character \tab Option for neutral species in \code{\link{nonideal}} [\code{bgamma0}] \cr
      \code{Berman} \tab character \tab User data file for mineral parameters in the Berman equations [\code{NA}] \cr
      \code{maxcores} \tab numeric \tab Maximum number of cores for parallel calculations with \code{\link{palply}} [\code{2}]
}

    \item \code{thermo$element}
  Dataframe containing the thermodynamic properties of elements taken from Cox et al., 1989 and Wagman et al., 1982.
  The standard molal entropy (\eqn{S}(\code{Z})) at 25 \degC and 1 bar for the \dQuote{element} of charge (\code{Z}) was calculated from \eqn{S}(H2,g) + 2\eqn{S}(\code{Z}) =  2\eqn{S}(H+), where the standard molal entropies of H2,g and H+ were taken from Cox et al., 1989.
  The mass of \code{Z} is taken to be zero.
  Accessing this data frame using \code{\link{mass}} or \code{\link{entropy}} will select the first entry found for a given element; i.e., values from Wagman et al., 1982 will only be retrieved if the properties of the element are not found from Cox et al., 1989.
      \tabular{lll}{
      \code{element}  \tab character  \tab Symbol of element\cr
      \code{state}  \tab character \tab Stable state of element at 25 \degC and 1 bar\cr
      \code{source} \tab character \tab Source of data\cr
      \code{mass}  \tab numeric \tab Mass of element (in natural isotopic distribution;\cr
      \tab \tab referenced to a mass of 12 for \eqn{^{12}}{12}C)\cr
      \code{s}   \tab numeric \tab Entropy of the compound of the element in its stable\cr
      \tab \tab state at 25 \degC and 1 bar (cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1})\cr
      \code{n}  \tab numeric \tab Number of atoms of the element in its stable\cr
      \tab \tab compound at 25 \degC and 1 bar
    }

    \item \code{thermo$obigt}

  This dataframe is a thermodynamic database of standard molal thermodynamic properties and equations of state parameters of species.
  Note the following database conventions:
  \itemize{
     \item The combination of \code{name} and \code{state} defines a species in \code{thermo$obigt}. A species can not be duplicated (this is checked when running \code{reset()}).
     \item English names of inorganic gases are used only for the gas state. The dissolved species is named with the chemical formula. Therefore, \code{info("oxygen")} refers to the gas, and \code{info("O2")} refers to the aqueous species.
     \item Properties of most aqueous species (\code{state} = \samp{aq}) are calculated using the revised Helgeson-Kirkham-Flowers (HKF) model (see \code{\link{hkf}}).
     \item Properties of aqueous species with an NA value of \code{Z} (the final column of thermo$obigt) are calculated using the Akinfiev-Diamond model (see \code{\link{AkDi}}).
     \item Properties of most non-aqueous species (liquids, gases, and minerals) are calculated using a heat capacity polynomial expression with up to six terms (see \code{\link{cgl}}).
     \item Properties of minerals with NA values of all heat capacity parameters are calculated using the Berman model (see \code{\link{berman}}).
  } 

  \samp{OrganoBioGeoTherm} is the name of a GUI program to use SUPCRT in Windows, produced in Harold C. Helgeson's Laboratory of Theoretical Geochemistry and Biogeochemistry at the University of California, Berkeley.
  The \acronym{OBIGT} database was originally developed for that program, and was the original basis for the database in CHNOSZ.
  There may be an additional meaning for the acronym: \dQuote{One BIG Table} of thermodynamic data.

  Each entry is referenced to one or two literature sources listed in \code{thermo$refs}.
  Use \code{\link{thermo.refs}} to look up the citation information for the references.
  See the vignette \emph{Thermodynamic data in CHNOSZ} for a complete description of the sources of data.
  The original \acronym{OBIGT} database was influenced by the \acronym{SUPCRT92} (Johnson et al., 1992) and \acronym{slop98.dat} data files (Shock et al., 1998), and the references in those files are included here.

In order to represent thermodynamic data for minerals with phase transitions, the higher-temperature phases of these minerals are represented as phase species that have states denoted by \samp{cr2}, \samp{cr3}, etc.
The standard molar thermodynamic properties at 25 \degC and 1 bar (\eqn{T_r}{Pr} and \eqn{P_r}{Pr}) of the \samp{cr2} phase species of minerals were generated by first calculating those of the \samp{cr} (lowest-T) phase species at the transition temperature (\eqn{T_{tr}}{Ttr}) and 1 bar then taking account of the volume and entropy of transition (the latter can be retrieved by combining the former with the Clausius-Clapeyron equation and values of \eqn{(dP/dT)} of transitions taken from the \acronym{SUPCRT92} data file) to calculate the standard molar entropy of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr}, and taking account of the enthalpy of transition (\eqn{{\Delta}H^{\circ}}{DeltaH0}, taken from the \acronym{SUPCRT92} data file) to calculate the standard molar enthalpy of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr}.
The standard molar properties of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr} and 1 bar calculated in this manner were combined with the equations-of-state parameters of the species to generate values of the standard molar properties at 25 \degC and 1 bar.
This process was repeated as necessary to generate the standard molar properties of phase species represented by \samp{cr3} and \samp{cr4}, referencing at each iteration the previously calculated values of the standard molar properties of the lower-temperature phase species (i.e., \samp{cr2} and \samp{cr3}).
A consequence of tabulating the standard molar thermodynamic properties of the phase species is that the values of \eqn{(dP/dT)} and \eqn{{\Delta}H^{\circ}}{DeltaH0} of phase transitions can be calculated using the equations of state and therefore do not need to be stored in the thermodynamic database.
However, the transition temperatures (\eqn{T_{tr}}{Ttr}) generally can not be assessed by comparing the Gibbs energies of phase species and are tabulated in the database.

The identification of species and their standard molal thermodynamic properties at 25 \degC and 1 bar are located in the first 12 columns of \code{thermo$obigt}:

    \tabular{lll}{
      \code{name}     \tab character \tab Species name\cr
      \code{abbrv}    \tab character \tab Species abbreviation\cr
      \code{formula}  \tab character \tab Species formula\cr
      \code{state}    \tab character \tab Physical state\cr
      \code{ref1}  \tab character \tab Primary source\cr
      \code{ref2}  \tab character \tab Secondary source\cr
      \code{date}     \tab character \tab Date of data entry (formatted as in SUPCRT92)\cr
      \code{G}        \tab numeric   \tab Standard molal Gibbs energy of formation\cr
      \tab \tab from the elements (cal mol\eqn{^{-1}}{^-1})\cr
      \code{H}        \tab numeric   \tab Standard molal enthalpy of formation\cr
      \tab \tab from the elements (cal mol\eqn{^{-1}}{^-1})\cr
      \code{S}        \tab numeric   \tab Standard molal entropy (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{Cp}       \tab numeric   \tab Standard molal isobaric heat capacity (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{V}	      \tab numeric   \tab Standard molal volume (cm\eqn{^3} mol\eqn{^{-1}}{^-1})
    }


The meanings of the remaining columns depend on the model used for a particular species (see database conventions above).
The names of these columns are compounded from those of the parameters in the HKF equations of state and general heat capacity polynomial; for example, column 13 is named \code{a1.a}.
Scaling of the values by orders of magnitude is adopted for some of the parameters, following common usage in the literature.

Columns 13-20 for aqueous species (parameters in the revised HKF equations of state):

    \tabular{lll}{
      \code{a1} \tab numeric \tab \eqn{a_1\times10}{a1 * 10} (cal mol\eqn{^{-1}}{^-1} bar\eqn{^{-1}}{^-1})\cr
      \code{a2} \tab numeric \tab \eqn{a_2\times10^{-2}}{a2 * 10^{-2}} (cal mol\eqn{^{-1}}{^-1})\cr
      \code{a3} \tab numeric \tab \eqn{a_3}{a3} (cal K mol\eqn{^{-1}}{^-1} bar\eqn{^{-1}}{^-1})\cr
      \code{a4} \tab numeric \tab \eqn{a_4\times10^{-4}}{a4 * 10^-4} (cal mol\eqn{^{-1}}{^-1} K)\cr
      \code{c1} \tab numeric \tab \eqn{c_1}{c1} (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{c2} \tab numeric \tab \eqn{c_2\times10^{-4}}{c2 * 10^-4} (cal mol\eqn{^{-1}}{^-1} K)\cr
      \code{omega} \tab numeric \tab \eqn{\omega\times10^{-5}}{omega * 10^-5} (cal mol\eqn{^{-1}}{^-1})\cr
      \code{Z}  \tab numeric \tab Charge
    }

Columns 13-20 for crystalline, gas and liquid species (\eqn{Cp=a+bT+cT^{-2}+dT^{-0.5}+eT^2+fT^{\lambda}}{Cp = a + bT + cT^-2 + dT^-0.5 + eT^2 + fT^lambda}).

    \tabular{lll}{
      \code{a} \tab numeric \tab \eqn{a} (cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1})\cr
      \code{b} \tab numeric \tab \eqn{b\times10^3}{b * 10^3} (cal K\eqn{^{-2}}{^-2} mol\eqn{^{-1}}{^-1})\cr
      \code{c} \tab numeric \tab \eqn{c\times10^{-5}}{c * 10^-5} (cal K mol\eqn{^{-1}}{^-1})\cr
      \code{d} \tab numeric \tab \eqn{d} (cal K\eqn{^{-0.5}}{^-0.5} mol\eqn{^{-1}}{^-1})\cr
      \code{e} \tab numeric \tab \eqn{e\times10^5}{e * 10^5} (cal K\eqn{^{-3}}{^-3} mol\eqn{^{-1}}{^-1})\cr
      \code{f} \tab numeric \tab \eqn{f} (cal K\eqn{^{-\lambda-1}}{-lambda-1} mol\eqn{^{-1}}{^-1})\cr
      \code{lambda} \tab numeric \tab \eqn{\lambda}{lambda} (exponent on the \eqn{f} term)\cr
      \code{T} \tab numeric \tab Temperature of phase transition or upper\cr
      \tab \tab temperature limit of validity of extrapolation (K)
    }

Columns 13-20 for aqueous species using the Akinfiev-Diamond model. Note that the \code{c} column is used to store the \eqn{\xi}{xi} parameter, and that \code{Z} must be NA to activate the code for this model. The remaining columns are not used.

    \tabular{lll}{
      \code{a} \tab numeric \tab \eqn{a} (cm\S{3} g\S{-1}) \cr
      \code{b} \tab numeric \tab \eqn{b} (cm\S{3} K\S{0.5} g\S{-1}) \cr
      \code{c} \tab numeric \tab \eqn{\xi}{xi} \cr
      \code{d} \tab numeric \tab \eqn{XX1} NA \cr
      \code{e} \tab numeric \tab \eqn{XX2} NA \cr
      \code{f} \tab numeric \tab \eqn{XX3} NA \cr
      \code{lambda} \tab numeric \tab \eqn{XX4} NA \cr
      \code{Z} \tab numeric \tab \eqn{Z} NA \cr
    }

    \item \code{thermo$refs}
    Dataframe of references to sources of thermodynamic data.
    \tabular{lll}{
      \code{key} \tab character \tab Source key\cr
      \code{author} \tab character \tab Author(s)\cr
      \code{year} \tab character \tab Year\cr
      \code{citation} \tab character \tab Citation (journal title, volume, and article number or pages; or book or report title)\cr
      \code{note} \tab character \tab Short description of the compounds or species in this data source\cr
      \code{URL} \tab character \tab URL\cr
    }

    \item \code{thermo$buffers}

    Dataframe which contains definitions of buffers of chemical activity. Each named buffer can be composed of one or more species, which may include any species in the thermodynamic database and/or any protein. The calculations provided by \code{\link{buffer}} do not take into account phase transitions of minerals, so individual phase species of such minerals must be specified in the buffers.
    \tabular{lll}{
      \code{name} \tab character \tab Name of buffer\cr
      \code{species} \tab character \tab Name of species\cr
      \code{state} \tab character \tab Physical state of species\cr
      \code{logact} \tab numeric \tab Logarithm of activity (fugacity for gases)
    }

    \item \code{thermo$protein}
    Data frame of amino acid compositions of selected proteins. Most of the compositions were taken from the SWISS-PROT/UniProt online database (Boeckmann et al., 2003) and the protein and organism names usually follow the conventions adopted there. In some cases different isoforms of proteins are identified using modifications of the protein names; for example, \samp{MOD5.M} and \code{MOD5.N} proteins of \samp{YEAST} denote the mitochondrial and nuclear isoforms of this protein. See \code{\link{pinfo}} to search this data frame by protein name, and other functions to work with the amino acid compositions.
    \tabular{lll}{
      \code{protein} \tab character \tab Identification of protein\cr
      \code{organism} \tab character \tab Identification of organism\cr
      \code{ref} \tab character \tab Reference key for source of compositional data\cr
      \code{abbrv} \tab character \tab Abbreviation or other ID for protein\cr
      \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
      \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
    }

    \item \code{thermo$groups}
    This is a dataframe with 22 columns for the amino acid sidechain, backbone and protein backbone groups ([Ala]..[Tyr],[AABB],[UPBB]) whose rows correspond to the elements C, H, N, O, S. It is used to quickly calculate the chemical formulas of proteins that are selected using the \code{iprotein} argument in \code{\link{affinity}}.

    \item \code{thermo$basis}
    Initially \code{NULL}, reserved for a dataframe written by \code{\link{basis}} upon definition of the basis species. The number of rows of this dataframe is equal to the number of columns in \dQuote{...} (one for each element).
     \tabular{lll}{
        \code{...} \tab numeric \tab One or more columns of stoichiometric\cr
        \tab \tab coefficients of elements in the basis species\cr
        \code{ispecies} \tab numeric \tab Rownumber of basis species in \code{thermo$obigt}\cr
        \code{logact} \tab numeric \tab Logarithm of activity or fugacity of basis species\cr
        \code{state} \tab character \tab Physical state of basis species\cr
     }

    \item \code{thermo$species}
    Initially \code{NULL}, reserved for a dataframe generated by \code{\link{species}} to define the species of interest. The number of columns in \dQuote{...} is equal to the number of basis species (i.e., rows of \code{thermo$basis}).
    \tabular{lll}{
       \code{...} \tab numeric \tab One or more columns of stoichiometric\cr
       \tab \tab coefficients of basis species in the species of interest\cr
       \code{ispecies} \tab numeric \tab Rownumber of species in \code{thermo$obigt}\cr
       \code{logact} \tab numeric \tab Logarithm of activity or fugacity of species\cr
       \code{state} \tab character \tab Physical state of species\cr
       \code{name} \tab character \tab Name of species\cr
    }

    \item \code{thermo$stoich}
    A precalculated stoichiometric matrix for the default database. This is a matrix, not a data frame, and as such can accept duplicated row names, corresponding to chemical formulas of the species. See \code{\link{retrieve}}, and the first test in \code{testthat/test-retrieve.R} for how to update this.
    \tabular{lll}{
       \code{rownames} \tab character \tab Chemical formulas from \code{thermo$obigt}\cr
       \code{...} \tab numeric \tab Stoichiometry, one column for each element present in any species
    }


  }  % end of itemize with long descriptions

} % end of format

\seealso{
Other data files, including those supporting the examples and vignettes, are documented separately at \code{\link{extdata}}. 
}

\examples{
## where are the data files in CHNOSZ?
system.file("extdata", package="CHNOSZ")
# what files make up OBIGT?
# nb. the .csv.xz files are loaded by default,
# and the .csv files have optional data that can be loaded with add.obigt()
dir(system.file("extdata/OBIGT", package = "CHNOSZ"))

## exploring thermo$obigt
# what physical states there are
unique(thermo()$obigt$state)
# formulas of ten random species
n <- nrow(thermo()$obigt)
thermo()$obigt$formula[runif(10)*n]
}

\references{
Cox, J. D., Wagman, D. D. and Medvedev, V. A., eds. (1989) \emph{CODATA Key Values for Thermodynamics}. Hemisphere Publishing Corporation, New York, 271 p. \url{http://www.worldcat.org/oclc/18559968}

Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\degC. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{https://doi.org/10.1016/0098-3004(92)90029-Q}

Shock, E. L. et al. {1998} \emph{SLOP98.dat} (computer data file). http://geopig.asu.edu/supcrt92_data/slop98.dat, accessed on 2005-11-05; moved to http://geopig.asu.edu/?q=tools.

Wagman, D. D., Evans, W. H., Parker, V. B., Schumm, R. H., Halow, I., Bailey, S. M., Churney, K. L. and Nuttall, R. L. (1982) The NBS tables of chemical thermodynamic properties. Selected values for inorganic and C\eqn{_1}{1} and C\eqn{_2}{2} organic substances in SI units. \emph{J. Phys. Chem. Ref. Data} \bold{11} (supp. 2), 1--392. \url{https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf}
}

\concept{Thermodynamic data}