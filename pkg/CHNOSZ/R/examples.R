# CHNOSZ/examples.R
# run examples from the help files, 
# and a function containing extra examples

examples <- function(do.png=FALSE) {
  # run all the examples in CHNOSZ documentation
  .ptime <- proc.time()
  topics <- c("CHNOSZ-package", "thermo", "sideeffects", "examples",
    "util.args", "util.array", "util.blast", "util.character", 
    "util.data", "util.expression", "util.fasta", "util.formula", "util.matrix", 
    "util.misc", "util.program",
    "util.seq", "util.stat", "util.units", "taxonomy", "info", "protein.info", "hkf", "water", "subcrt",
    "makeup", "basis", "swap.basis", "species", "affinity", "util.affinity", "equil.boltzmann", 
    "diagram", "buffer", "iprotein", "protein", "ionize.aa", "more.aa", "read.expr", "revisit", 
    "findit", "transfer", "anim", "EOSregress", "wjd")
  plot.it <- FALSE
  if(is.character(do.png))
    png(paste(do.png,"%d.png",sep=""),width=700,height=700,pointsize=18)
  else if(do.png) plot.it <- TRUE
  for(i in 1:length(topics)) {
    if(plot.it) png(paste(topics[i],"%d.png",sep=""),width=700,height=700,pointsize=18)
    myargs <- list(topic=topics[i],ask=FALSE)
    do.call(example,myargs)
    if(plot.it) dev.off()
  }
  if(is.character(do.png)) dev.off()
  # at the end we attempt to restore the old par() (active as of the first call of thermo.plot.new)
  par(thermo$opar)
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
}

longex <- function(which=c("sources", "NaCl", "copper", "cordierite", 
  "phosphate", "nucleobase", "pie", "orp", "findit", "CO2Ac", "nonideal", "TPX")) {

  # extra examples for fun
  # character case doesn't matter
  which <- tolower(which)
  # run more than one if wanted
  if(length(which) > 1) {
    for(i in 1:length(which)) out <- longex(which[i])
    return(invisible(out))
  }

  if(which=="sources") {
    ## cross-checking sources
    # the reference sources
    ref.source <- thermo$refs$key
    # sources of elemental data
    element.source <- thermo$element$source
    # sources in the primary thermodynamic database
    os1 <- thermo$obigt$ref1
    os2 <- thermo$obigt$ref2
    # sources also in the supplemental database (OBIGT-2.csv)
    add.obigt()
    os3 <- thermo$obigt$ref1
    os4 <- thermo$obigt$ref2
    data(thermo)
    # all of the thermodynamic data sources - some of them might be NA
    obigt.source <- unique(c(os1,os2,os3,os4))
    obigt.source <- obigt.source[!is.na(obigt.source)]
    # sources of protein compositions
    protein.source <- thermo$protein$ref
    # sources of stress response proteins
    stress.source <- as.character(thermo$stress[2,])
    # if the sources are all accounted for 
    # these all produce character(0)
    print("missing these sources for elemental properties:")
    print(unique(element.source[!(element.source %in% ref.source)]))
    print("missing these sources (1) for thermodynamic properties:")
    print(unique(obigt.source[!(obigt.source %in% ref.source)]))
    print("missing these sources for protein compositions:")
    print(unique(protein.source[!(protein.source %in% ref.source)]))
    print("missing these sources for stress response experiments:")
    print(unique(stress.source[!(stress.source %in% ref.source)]))
    # determine if all the reference sources are cited
    my.source <- c(element.source,obigt.source,protein.source,stress.source)
    # this should produce character(0)
    print("these sources are present but not cited:")
    print(ref.source[!(ref.source %in% my.source)])

  } else if(which=="nacl") {

    ## NaCl dissocation logK f(T,P)
    ## after Shock et al., 1992, Fig. 1
    ## (Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) 
    ##  Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: 
    ##  Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 degrees C and 5 kbar. 
    ##  J. Chem. Soc. Faraday Trans. 88, 803-826. http://dx.doi.org/10.1039/FT9928800803 )
    species <- c("NaCl", "Na+", "Cl-")
    coeffs <- c(-1, 1, 1)
    # start a new plot and show the experimental logK
    thermo.plot.new(xlim=c(0, 1000), ylim=c(-5.5, 1),
      xlab=axis.label("T"), ylab=axis.label("logK"))
    expt <- read.csv(system.file("extdata/cpetc/SOJSH.csv", 
      package="CHNOSZ"), as.is=TRUE)
    points(expt$T,expt$logK, pch=expt$pch)
    # we'll be at 9 distinct pressure conditions, including Psat
    P <- c(list("Psat"), as.list(seq(500, 4000, by=500)))
    # for each of those what's the range of temperature
    T <- list()
    # T > 350 degC at Psat is possibly inappropriate; see "Warning" of subcrt.Rd 
    T[[1]] <- seq(0, 370, 5)
    T[[2]] <- seq(265, 465, 5)
    T[[3]] <- seq(285, 760, 5)
    T[[4]] <- seq(395, 920, 5)
    T[[5]] <- T[[6]] <- T[[7]] <- T[[8]] <- T[[9]] <- seq(400, 1000, 5)
    # calculate and plot the logK
    logK <- numeric()
    for(i in 1:length(T)) {
      s <- subcrt(species, coeffs, T=T[[i]], P=P[[i]])
      lines(s$out$T, s$out$logK)
      # keep the calculated values for each experimental condition
      iexpt <- which(P[[i]]==expt$P)
      Texpt <- expt$T[iexpt]
      logK <- c(logK, splinefun(s$out$T, s$out$logK)(Texpt))
    }
    legend("bottomleft",pch=unique(expt$pch),
      legend=c(unique(expt$source),tail(expt$source,1)))
    title(main=paste("NaCl(aq) = Na+ + Cl-\n",
      "Psat and 500-4000 bar, after Shock et al., 1992"))
    # where do we diverge most from experiment?
    imaxdiff <- which.max(abs(logK - expt$logK))
    stopifnot(all.equal(c("Psat", 347.7),
      as.character(expt[imaxdiff,1:2])))
    # what's our average divergence?
    stopifnot(mean(abs(logK - expt$logK)) < 0.09)

  } else if(which=="copper") {

    ## Eh-pH diagrams for copper-water-glycine
    ## After Fig. 2 of Aksu and Doyle, 2001
    ## (Aksu, S. and Doyle, F. M., 2001. Electrochemistry of copper in aqueous glycine 
    ## solutions. J. Electrochem. Soc., 148, B51-B57. doi:10.1149/1.1344532)
    ##  We need to add some species and change some Gibbs energies.
    # update rows of the database
    i <- info(c("Cu(Gly)+","glycinate","e-","H+"))
    # this was written before mod.obigt() was around... could use that instead
    n <- nrow(thermo$obigt <<- rbind(thermo$obigt,thermo$obigt[rep(i[1],2),]))
    thermo$obigt$name[n-1] <<- "Cu(Gly)2-"
    thermo$obigt$name[n] <<- "HCu(Gly)+2"
    thermo$obigt$formula[n-1] <<- as.chemical.formula(makeup(c(i[1],i[2],i[3]), sum=TRUE))
    thermo$obigt$formula[n] <<- as.chemical.formula(makeup(c(i[1],i[4]), sum=TRUE))
    # In Fig 2b, total log activities of Cu (Cu_T) 
    # and glycine (L_T) are -4 and -1
    basis(c("Cu+2","H2O","H+","e-","glycine","CO2"),c(999,0,999,999,-1,999))
    # solid species
    species(c("copper","cuprite","tenorite"))
    # aqueous species
    species(c("glycinium","glycine","glycinate","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"),-4)
    ispecies <- species()$ispecies
    # update the Gibbs energies using A&D's Table 1 and Table II
    logK <- c(convert(convert(c(0,-146,-129.7,-384.061,-370.647,-314.833,
      49.98,65.49,-183.6,-258.5,-298.2)*1000,"cal"),"logK"),15.64,10.1,2.92) 
    # do it in order so later species take account of prev. species' values
    for(i in 1:length(logK)) {
      G <- convert(logK[i],"G")
      if(i==12) G <- G + thermo$obigt$G[ispecies[8]] + 
        2*thermo$obigt$G[ispecies[6]]
      if(i==13) G <- G + thermo$obigt$G[ispecies[7]] + 
        2*thermo$obigt$G[ispecies[6]]
      if(i==14) G <- G + thermo$obigt$G[ispecies[11]]
      thermo$obigt$G[ispecies[i]] <- G
    }  # done with changing Gibbs free energies!
    # we have to get some leftovers out of there or diagram() gets confused
    species(c("glycinium","glycine","glycinate"),delete=TRUE)
    # make a plot to see if it's working
    ispecies <- ispecies[-(1:6)]
    afun <- function(cu,gly) {
      # from above: our fifth basis species is glycine(-ate,-ium)
      basis(rownames(basis())[5],gly)
      t <- match(ispecies,species()$ispecies)
      species(t,cu)
      affinity(pH=c(0,16),Eh=c(-0.6,1.0))
    }
    diagram(afun(-4,-1))
    title(main=paste("Aqueous Copper + Glycine, 25 deg C, 1 bar",
      "After Aksu and Doyle, 2001 Fig. 2b",sep="\n"))
    # What's missing? Try glycinate not glycine in reactions at ph > ~9.8
    basis(c("Cu+2","H2O","H+","e-","glycinate","CO2"),
      c(999,0,999,999,-2,999))
    species(c("copper","cuprite","tenorite","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"))
    loga_Cu <- -4
    loga_Gly <- -1
    diagram(afun(loga_Cu,loga_Gly),fill=NULL,col="blue",
      names=species()$name,col.names="blue",add=TRUE)
    water.lines()
    # the glycine ionization constants could be calculated using
    # subcrt, here they are taken from A&D Table II
    abline(v=c(2.35,9.778),lty=3)
    # now put glycinium (low pH) in the basis
    basis(c("Cu+2","H2O","H+","e-","glycinium","CO2"),c(999,0,999,999,-2,999))
    species(c("copper","cuprite","tenorite","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"))
    diagram(afun(loga_Cu,loga_Gly),fill=NULL,col="green",
      names=NULL,col.names="green",add=TRUE)
    # let's forget our changes to 'thermo' so that examples
    # below that use glycine will work as expected
#    data(thermo)

  } else if(which=="cordierite") {

    ### 1-D property plot
    ## for hydrous cordierite = cordierite + H2O
    ## after Helgeson et al., 1978
    ## (Summary and critique of the thermodynamic properties of 
    ## rock-forming minerals. Am. J. Sci., 278-A, 1-229)
    basis(c("cordierite,hydrous","Mg+2","SiO2","H2O","O2","H+"))
    species("cordierite")
    # water.SUPCRT92 can only get us up to 5000 bar
    # (lines for 7000 and 10000 bar are in the original diagram)
    P <- c(1,2,3,5)*1000
    col <- rainbow(length(P))
    for(i in 1:length(P)) {
      a <- affinity(property="logK",T=c(20,800),P=P[i])
      diagram(a,add=(i!=1),ylim=c(-4,2),legend.x=NULL,
        col=col[i],main="")
    }
    legend("topright",lty=1,col=col,legend=paste(P,"bar"))
    title(main=paste("hydrous cordierite = cordierite + H2O",
      "After Helgeson et al., 1978",sep="\n"),cex.main=0.9)

  } else if(which=="phosphate") {

    ## speciation of phosphate as a function of ionic strength
    opar <- par(mfrow=c(2, 1))
    basis("CHNOPS+")
    T <- c(25, 100)
    species(c("PO4-3", "HPO4-2", "H2PO4-"))
    d25 <- diagram(affinity(IS=c(0, 0.14), T=T[1]), ylim=c(-3.0, -2.6), legend.x=NULL)
    d100 <- diagram(affinity(IS=c(0, 0.14), T=T[2]), ylim=c(-3.0, -2.6), add=TRUE, col="red")  
    title(main="Non-ideality model for phosphate species")
    dp <- describe.property(c("pH", "T", "T"), c(7, T))
    legend("topright", lty=c(NA, 1, 1), col=c(NA, "black", "red"), legend=dp)
    text(0.07, -2.76, expr.species("HPO4-2"))
    text(0.07, -2.90, expr.species("H2PO4-"))
    # the crossing points of the logarithms of activity at the two temperatures
    # (it's higher IS at higher temperature)
    x25 <- which.min(abs(d25$logact[[3]] - d25$logact[[2]]))
    stopifnot(all.equal(x25, 27))
    x100 <- which.min(abs(d100$logact[[3]] - d100$logact[[2]]))
    stopifnot(all.equal(x100, 45))
    ## phosphate predominance f(IS,pH)
    d <- diagram(affinity(IS=c(0, 0.14), pH=c(6, 13), T=T[1]), fill=NULL)
    diagram(affinity(IS=c(0, 0.14), pH=c(6, 13), T=T[2]), 
      add=TRUE, names=NULL, col="red")
    par(opar)
    # the most stable species are PO4-3, HPO4-2, H2PO4- with decreasing pH,
    # at any ionic strength (in the range used here)
    stopifnot(all.equal(unique(d$out[1:128, 1]), 1:3))
    stopifnot(all.equal(unique(d$out[1:128, 128]), 1:3))

  } else if(which=="nucleobase") {

    ## Nucleobase - Amino Acid Interaction Eh-H2O
    # for this example we try a unique basis definition
    basis(c("CO2","H2O","glutamine","e-","H+"),c(-3,0,-3,0,-7))
    species(c("uracil","cytosine","adenine","guanine",
      "phenylalanine","proline","lysine","glycine"),"aq")
    # this loaded four nucleobases and four related amino acids
    # (coded for by the homocodon triplets)
    # check out the predominance diagrams
    a.1 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0))
    diagram(a.1,fill=NULL)
    # overlay a different temperature
    a.2 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0),T=100)
    diagram(a.2,col="red",add=TRUE,names=NULL)
    # start make a title for the plot
    tb <- thermo$basis   # includes activities of basis species
    # exclude those that are on the axes
    tb <- tb[!((rownames(tb) %in% c("e-","H2O"))),]
    title(main="Nucleobases and amino acids; P=Psat")
    dp <- describe.property(c("T", "T"), c(25, 100))
    db <- describe.basis(tb)
    legend("bottomleft", lty=c(1, 1, NA, NA, NA), col=c("black","red"), legend=c(dp, db))

  } else if(which=="pie") {

    #    This is an attempt to predict some characteristics of the relative
    #    abundances of organisms that are depicted in the pie charts shown by
    #    Spear et al., 2005 (Spear, J. R., Walker, J. J., McCollom, T. M. and
    #    Pace, N. R. Hydrogen and bioenergetics in the Yellowstone geothermal
    #    ecosystem. Proc. Natl. Acad. Sci. U. S. A., 102, 2555-2560, 2005.).
    #    For each type of organism present, we use a single model protein. We
    #    take the reported relative abundances of the four most abundant
    #    organisms to generate an assemblage of proteins that buffers the
    #    activies of CO2, H2O and NH3. Then we calculate relative abundances of
    #    all proteins at different values of oxygen fugacity and H2S activity
    #    and display them on pie charts that can be compared with organismal
    #    abundances. The results show that it is possible to predict conditions
    #    that foster high (many slices in the pie diagrams) or low (few slices)
    #    diversity, even if the make up of the community is not perfectly
    #    replicated.  jmd 2008-2009
    ### Setup
    opar <- par(no.readonly=TRUE)
    # Bacterial species and abundances from Spear et al., 2005 (Fig. 4):
    O1 <- c('Aquificales','Thermotogales','Bacillus','Thermodesulfobacteria')
    O2 <- c('Thermus/Deinococcus','Hydrogenobacter', 'Hydrogenobaculum','Bacteroidetes')
    O3 <- c('a-proteobacterium','d-proteobacterium','g-proteobacterium','OD-1')
    ORGANISMS <- c(O1,O2,O3) 
    # Which species are most abundant in low and high sulfur environments; 
    # the first three of each are aquificales, thermotogales, thermodesulfobacteria:
    lowS.which <- c(1,2,4,3,5)
    lowS.values <- c(75,12,4,7,2)
    highS.which <- c(1,2,4,11,6,9,7,8,10,13)
    highS.values <- c(63,2,10,2,5,2,9,2,3,1)
    # What chemical species we use to represent each organism:
    P1 <- c('ACCA_AQUAE','A8F7V4_THELT','RBL1_THIDA')
    P2 <- c('Q93EV7_9BACT','ACCA_DEIRA','Q05KD0_HYDTH','A7WGI1_9AQUI')
    P3 <- c('Q64VW6_BACFR','ACCA_CAUCR','A1VC70_DESVV','ACCA_PSEAE')
    PROTEINS <- c(P1,P2,P3)
    ### Function definitions
    # To make pie charts for calculated abundances:
    plot.pie.calc <- function(which="high",T=25,main='') {
      # first clean up the buffer definition in case we have
      # been run already
      thermo$buffers <<- thermo$buffers[thermo$buffers$name!='PROTEINS',]
      # we take four predominant proteins from SWM+05
      myprot <- PROTEINS[get(paste(which,"S.which",sep=""))][1:4]
      mypercent <- get(paste(which,"S.values",sep=""))[1:4]
      # use these four proteins to create a buffer
      mybufprot <- paste(myprot,'RESIDUE',sep='.')
      mod.buffer('PROTEINS',mybufprot,'aq',log10(mypercent/100))
      # our species are the residues of all proteins
      species(delete=TRUE)
      # 20120520 the next 3 lines take the place of previous protein.residue call
      aa <- ip2aa(PROTEINS, residue=TRUE)
      aa$organism <- paste(aa$organism, "RESIDUE", sep=".")
      add.protein(aa)
      species(paste(aa$protein, aa$organism, sep="_"))
      # assign the buffer to three basis species
      basis(c('CO2','H2O','NH3'),'PROTEINS')
      # calculate the buffered activities
      a <- affinity(return.buffer=TRUE,balance=1,T=T)
      # make the titles
      sub <- c2s(paste(names(a)[1:3],round(as.numeric(a)[1:3])))
      main <- paste('\n',main,'\n',sub)
      # set the total species activities to those in the buffer
      species(1:nrow(species()),-99)
      species(mybufprot,log10(mypercent/100))
      # get the activities back to numeric
      basis(names(a)[1:3],as.numeric(a)[1:3])
      thermo$basis$logact <<- as.numeric(thermo$basis$logact)
      # colors
      col <- rep('white',99)
      col[match(myprot,PROTEINS)] <- heat.colors(4)
      # calculate the distribution of species
      mylogaH2O <- thermo$basis$logact[rownames(thermo$basis)=='H2O']
      a <- affinity(H2O=c(mylogaH2O,mylogaH2O-1,2),T=T)
      logacts <- equilibrate(a, normalize=TRUE, balance=1)$loga.equil
      # assemble the names and logarithms of activities
      # of species that are above a minimum value
      names <- character()
      values <- numeric()
      cols <- character()
      logactmin <- -2
      for(i in 1:length(logacts)) {
        myvalue <- logacts[[i]][1]
        if(myvalue > logactmin) {
          names <- c(names,ORGANISMS[i])
          values <- c(values,myvalue)
          cols <- c(cols,col[i])
        }
      }
      # remove the logarithms
      values <- 10^values
      # sort them by abundance
      isort <- sort(values,index.return=TRUE,decreasing=TRUE)$ix
      names <- names[isort]
      values <- values[isort]
      cols <- cols[isort]
      # make a pie chart
      pie(values,names,clockwise=FALSE,main=main,col=cols,radius=0.7)
    }
    # To plot pie charts for observed abundances of organisms:
    plot.pie.obs <- function(which="low") {
      # the values from SWM+05
      names <- ORGANISMS[get(paste(which,"S.which",sep=""))]
      values <- get(paste(which,"S.values",sep=""))
      main <- paste("observed at",which,"H2S")
      # colors for the four dominant species
      mycol <- heat.colors(4)
      # colors for the rest
      mycol <- c(mycol,rep('white',length(names)-length(mycol)))
      # sort the values
      isort <- sort(values,index.return=TRUE,decreasing=TRUE)$ix
      values <- values[isort]
      names <- names[isort]
      mycol <- mycol[isort]
      pie(values,names,clockwise=FALSE,main=main,col=mycol,radius=0.7)
    }
    # To plot both types of pie diagrams (showing calculated protein 
    # activities and observed abundances of organisms) at a given temperature:
    plot.pie <- function(T=80) {
      # first deal with the layout
      layout(matrix(1:4,byrow=TRUE,nrow=2))
      opar <- par(mar=c(1,1,1,1))
      # basis definition
      basis(c('CO2','H2O','NH3','H2S','hydrogen','H+'))
      basis('pH',7)
      val <- function(text)
        round(as.numeric(thermo$basis$logact[rownames(thermo$basis)==text]))
      # now to plotting
      # low sulfur and relatively oxidizing
      basis(c('H2','H2S'),c(-9,-9))
      plot.pie.calc("low",T=T,main=paste("calculated for H2",val("H2"),"H2S",val("H2S")))
      label.plot('a')
      plot.pie.obs("low")
      label.plot('b')
      # high sulfur and relatively reducing
      basis(c('H2','H2S'),c(-5,-3))
      plot.pie.calc("high",T=T,main=paste("calculated for H2",val("H2"),"H2S",val("H2S")))
      label.plot('c')
      plot.pie.obs("high")
      label.plot('d')
      par(opar)
    }
    ### Now run it!
    # at 80 degrees C:
    plot.pie(80)
    # reset plot device
    par(opar)

  } else if(which=="orp") {

    # yell2010/orp.R 20100715 jmd
    # calculate the temperature dependence of 
    # potentials vs. SHE of various electrodes (Ag/AgCl)
    # and ORP standards (ZoBell, Light's, (tri)iodide) 
    # CHNOSZ provides functions subcrt() and convert() 
    # used in this example
    #require(CHNOSZ)
    # Bard et al.'s fit to the potential
    # (Bard, Parson, Jordan, Standard Potentials In Aqueous Solution, 1985)
    AgAgCl.Bard <- function(T,high.T=TRUE) {
      # we use the corrected high-T formula from wikipedia
      if(high.T) return(0.23737 - 5.3783e-4 * T - 2.3728e-6 * T^2 - 2.2671e-9 * (T+273))
      else return(0.23695 - 4.8564e-4 * T - 3.4205e-6 * T^2 - 5.869e-9 * (T+273))
    }
    # function to calculate the potential of Ag/AgCl vs. SHE
    # Ag(s) + Cl- = AgCl(s) + e-
    # logK = -pe - logaCl
    # pe = -logK - logmCl - loggamCl
    # ORP = RT/F * (logK - logmCl - loggamCl)
    AgAgCl <- function(T,mKCl=4) {
      # mKCl is the molality of KCl in the electrolyte
      # we take it as a first approximation to be equal to
      # the molality of Cl- (and to the ionic strength)
      logmCl <- log10(mKCl)
      # get the logK for the reaction
      logK <- subcrt(c("Ag","Cl-","AgCl","e-"),c(-1,-1,1,1),c("cr","aq","cr","aq"),T=T)$out$logK
      # get the activity coefficient for Cl-
      loggamCl <- subcrt("Cl-",T=T,IS=mKCl)$out[[1]]$loggam
      # get the pe for the solution
      pe <- -logK - logmCl - loggamCl
      # convert that to Eh
      Eh <- convert(pe,"Eh",T=convert(T,"K"))
      return(Eh)
    }
    ZoBell <- function(T) {
      # doesn't work very well because we ignore the
      # ferricyanide and ferrocyanide complexes
      # Fe+2 = Fe+3 + e-
      # logK = logaFe3 - logaFe2 - pe
      # get the logK for the reaction
      logK <- subcrt(c("Fe+2","Fe+3","e-"),c(-1,1,1),T=T)$out$logK
      # we use the recipe from standard methods (table 2580:II)
      # 1.4080 g K4Fe(CN)6.3H2O -> 0.0033333 mol Fe+2
      # 1.0975 g K3Fe(CN)6      -> 0.0033333 mol Fe+3
      # 7.4555 g KCl            -> 0.1 mol Cl-
      logmFe2 <- logmFe3 <- log10(0.0033333)
      # get the loggam for the iron species
      loggamFe2 <- subcrt("Fe+2",T=T,IS=1)$out[[1]]$loggam
      loggamFe3 <- subcrt("Fe+3",T=T,IS=1)$out[[1]]$loggam
      # get the pe for the solution
      pe <- -logK + logmFe3 + loggamFe3 - logmFe2 - loggamFe2
      # convert to Eh
      Eh <- convert(pe,"Eh",T=convert(T,"K"))
      return(Eh)
    }
    ZoBell.table <- function(T=NULL,which=NULL) {
      # oxidation-reduction potential of ZoBell's solution 
      # from Standard Methods for Water and Wastewater or YSI 
      # (interpolated and/or extrapolated as necessary)
      # standard methods (1997) table 2580:I
      Eh.T.SMW <- 1:30
      Eh.SMW <- c(0.481,0.479,0.476,0.474,0.472,0.47,0.468,0.465,0.463,0.461,
      0.459,0.457,0.454,0.452,0.45,0.448,0.446,0.443,0.441,0.439,0.437,
      0.435,0.432,0.43,0.428,0.426,0.424,0.421,0.419,0.417)
      # from YSI (2005):
      # Measuring ORP on YSI 6-Series Sondes: Tips, Cautions and Limitations
      # NOTE: these values are vs. Ag/AgCl (4 M KCl)
      Eh.T.YSI <- seq(-5,50,by=5)
      Eh.YSI <- c(267.0,260.5,254.0,247.5,241.0,234.5,228.0,221.5,215.0,208.5,202.0,195.5)/1000
      # spline function for each of the tables
      SMW <- splinefun(Eh.T.SMW,Eh.SMW)
      YSI <- splinefun(Eh.T.YSI,Eh.YSI)
      # just one of the tables
      Eh.fun <- get(which)
      Eh.T <- get(paste("Eh.T",which,sep="."))
      if(is.null(T)) T <- Eh.T
      return(data.frame(T=T,Eh=Eh.fun(T)))
    }
    Light <- function(T) {
      # this is going to look something like
      # Fe+2 = Fe+3 + e-
      # logK = logaFe3 - logaFe2 - pe
      # get the logK for the reaction
      logK <- subcrt(c("Fe+2","Fe+3","e-"),c(-1,1,1),T=T)$out$logK
      # we use the recipe from standard methods (table 2580:II)
      # 39.21 g Fe(NH4)2(SO4)2(H2O)6 -> 0.1 mol Fe+2
      # 48.22 g Fe(NH4)(SO4)2(H2O)12 -> 0.1 mol Fe+3
      logmFe2 <- logmFe3 <- log10(0.1)
      # get the loggam for the iron species
      loggamFe2 <- subcrt("Fe+2",T=T,IS=0.2)$out[[1]]$loggam
      loggamFe3 <- subcrt("Fe+3",T=T,IS=0.2)$out[[1]]$loggam
      # get the pe for the solution
      pe <- -logK + logmFe3 + loggamFe3 - logmFe2 - loggamFe2
      # convert to Eh
      Eh <- convert(pe,"Eh",T=convert(T,"K"))
      return(Eh)
    }
    Iodide.table <- function(T=NULL) {
      # oxidation-reduction potential of Thermo's iodide solution
      # from thermo instruction sheet 255218-001 (articlesFile_18739)
      T.Iodide <- seq(0,50,5)
      Eh.Iodide <- c(438,435,431,428,424,420,415,411,406,401,396)/1000
      Iodide <- splinefun(T.Iodide,Eh.Iodide)
      if(is.null(T)) T <- T.Iodide
      return(data.frame(T=T,Eh=Iodide(T)))
    }
    Iodide <- function(T) {
      # this is going to look something like
      # 3I- = I3- + 2e-
      # logK = -2pe + logaI3 - 3logaI
      # get the logK for the reaction
      logK <- subcrt(c("I-","I3-","e-"),c(-3,1,2),T=T)$out$logK
      # could the activities be 0.1 M ... or something else?
      logmI <- log10(2)
      logmI3 <- log10(0.01)
      # get the loggam for the iodine species
      loggamI <- subcrt("I-",T=T,IS=0.2)$out[[1]]$loggam
      loggamI3 <- subcrt("I3-",T=T,IS=0.2)$out[[1]]$loggam
      # get the pe for the solution
      pe <- ( -logK + logmI3 + loggamI3 - 3 * (logmI - loggamI) ) / 2
      # convert to Eh
      Eh <- convert(pe,"Eh",T=convert(T,"K"))
      return(Eh)
    }
    figure <- function() {
      # make some figures
      # the temperatures we're interested in
      # in degrees C
      T <- seq(0,100,5)
      # temperature-Eh diagram for various electrodes
      thermo.plot.new(ylim=c(0,0.8),xlim=c(0,100),
        ylab=axis.label("Eh"),xlab=axis.label("T"))
      # the Ag/AgCl electrode (Bard et al. fit)
      points(T,AgAgCl.Bard(T),pch=0)
      # the Ag/AgCl electrode (equilibrium calculations) 
      lines(T,AgAgCl(T))
      # ZoBell's solution (SMW table 2580)
      SMW <- ZoBell.table(which="SMW")
      points(SMW$T,SMW$Eh,pch=1)
      # ZoBell's solution (YSI tech report table)
      YSI <- ZoBell.table(which="YSI")
      # make these values referenced to SHE instead of Ag/AgCl
      Eh.YSI <- YSI$Eh + AgAgCl(YSI$T)
      points(YSI$T,Eh.YSI,pch=2)
      # Light's solution (equilibrium values)
      lines(T,Light(T))
      # Light's solution (at 25 degrees only)
      points(25,0.475 + 0.200,pch=3)
      # Thermo's I-/I3- solution
      Thermo <- Iodide.table()
      points(Thermo$T,Thermo$Eh,pch=4)
      # calculated I-/I3- values
      lines(T,Iodide(T))
      # add some labels
      text(c(30,30,30,50),c(0.72,0.5,0.35,0.25),
        c("Light","ZoBell","(Tri)Iodide","Ag/AgCl"))
      title(main="Potentials vs SHE")
    }
    # finally, make the plot
    figure()

  } else if(which=="findit") {

    opar <- par(mfrow=c(2,2))
    # an organic example: 
    # find chemical activities where metastable activities of
    # selected proteins in P. ubique have high correlation
    # with a lognormal distribution (i.e., maximize r of q-q plot)
    f <- system.file("extdata/fasta/HTCC1062.faa.xz",package="CHNOSZ")
    # search for three groups of proteins
    myg <- c("ribosomal","nucle","membrane")
    g <- lapply(myg,function(x) grep.file(f,x))
    # note that some proteins match more than one search term
    uug <- unique(unlist(g))
    # read their amino acid compositions from the file
    aa <- read.fasta(f,uug)
    # add these proteins to thermo$protein
    ip <- add.protein(aa)
    # load a predefined set of uncharged basis species
    # (speeds things up as we won't model protein ionization)
    basis("CHNOS")
    # make colors for the diagram
    rgbargs <- lapply(1:3,function(x) as.numeric(uug %in% g[[x]]))
    col <- do.call(rgb,c(rgbargs,list(alpha=0.5)))
    # get point symbols (use 1,2,4 and their sums)
    pch <- colSums(t(list2array(rgbargs)) * c(1,2,4))
    # plot 1: calculated logarithms of chemical activity
    # as a function of logfO2 ... a bundle of curves near logfO2 = -77
    a <- affinity(O2=c(-90,-60),iprotein=ip)
    e <- equilibrate(a, loga.balance=0)
    d <- diagram(e, col=col)
    title(as.expression("Selected proteins in"~italic("Pelagibacter ubique")))
    legend("bottomleft", lty=c(1,1,1), col=unique(col), legend=myg, bg="white")
    db <- describe.basis(ibasis=c(2, 1, 3))
    legend("bottomright", legend=db, bg="white")
    # plot 2: calculate q-q correlation coefficient
    # the lognormal distribution is favored near logfO2 = -73.6
    r <- revisit(d,"qqr")
    title(main="correlation with a normal distribution")
    # plot 3: findit... maximize qqr as a function of O2-H2O-NH3-CO2
    # it shows an optimum at low logaH2O, logaNH3
    f1 <- findit(list(O2=c(-106,-75),H2O=c(-40,-20),CO2=c(-20,10),NH3=c(-15,0)),
      "qqr",iprotein=ip,niter=8)
    title(main="searching 4-D chemical activity space")
    # plot 5: q-q plot at the final loga O2, H2O, CO2, NH3
    # higher correlation coefficient than plot 3
    a <- affinity(iprotein=ip)
    e <- equilibrate(a, loga.balance=0)
    qqr5 <- revisit(e, "qqr",pch=pch)$H
    legend("topleft",pch=c(1,2,4),legend=myg)
    db <- describe.basis(ibasis=c(5, 2, 1, 3))
    legend("bottomright", legend=db)
    # plot 5: trajectory of O2, H2O, CO2, NH3, and the
    # q-q correlation coefficient in the search
    #plot(f1,mar=c(2,5,1,1),mgp=c(4,1,0))
    par(opar)

  } else if(which=="co2ac") {

    # one can solve for the logact of a 
    # basis species using the 'what' argument of diagram
    basis("CHNOS")
    basis("CO2",999)
    species("acetic acid",-3)
    a <- affinity(O2=c(-85,-70,4),T=c(25,100,4))
    # hacking to write a title with formulas and subscripts
    lCO2 <- axis.label("CO2")[[1]]
    lAC <- expr.species("CH3COOH")[[1]]
    main <- substitute(a~~b~~c,list(a=lCO2,b="buffered by",
      c="acetic acid"))
    d <- diagram(a,what="CO2",main=main)
    species(1,-10)
    a <- affinity(O2=c(-85,-70,4),T=c(25,100,4))
    d <- diagram(a,what="CO2",add=TRUE,lty=2)
    # add a legend
    ltext <- c(axis.label("CH3COOH"),-3,-10)
    lty <- c(NA,1,2)
    legend("topright",legend=ltext,lty=lty,bg="white")
    # do return.buffer and diagram(what) give the same results?
    and <- as.numeric(d$logact[[1]])
    basis("CO2","AC")
    mod.buffer("AC",logact=-10)
    a.buffer <- affinity(O2=c(-85,-70,4),T=c(25,100,4),return.buffer=TRUE)
    ana <- as.numeric(unlist(a.buffer[[1]]))
    stopifnot(all.equal(ana,and))

  } else if(which=="nonideal") {

    ### non-ideality calculations -- activity coefficients of 
    ### aqueous species as a function of charge, temperature,
    ### and ionic strength -- after Alberty, 2003 
    ## p. 16 Table 1.3  apparent pKa of acetic acid with
    ## changing ionic strength
    subcrt(c("acetic acid","acetate","H+"),c(-1,1,1),
      IS=c(0,0.1,0.25),T=25,property="logK")
    # note that these *apparent* values of G and logK approach
    # their *standard* counterparts as IS goes to zero.
    ## p. 95: basis and elemental stoichiometries of species 
    ## (a digression here from the nonideality calculations) 
    # note coefficient of O2 and NH3 will be zero for these species
    basis(c("ATP-4","H+","H2O","HPO4-2","O2","NH3"))
    # cf Eq. 5.1-33: (basis composition) 
    species(c("ATP-4","H+","H2O","HPO4-2","ADP-3","HATP-3","HADP-2","H2PO4-"))
    lb <- nrow(basis())
    # cf Eq. 5.1-32: (elemental composition)
    as.matrix(species()[,1:lb]) %*% as.matrix(basis()[,1:lb]) 
    ## p. 273-275: activity coefficient (gamma)
    ## as a function of ionic strength and temperature
    ## (as of 20080304, these do look quantitatively different 
    ## from the plots in Alberty's book.)
    iplotfun <- function(T,col,add=TRUE) {
      IS <- seq(0,0.25,0.0025)
      s <- subcrt(c("H2PO4-","HADP-2","HATP-3","ATP-4"),IS=IS,grid="IS",T=T)
      if(!add) thermo.plot.new(xlim=range(IS),ylim=c(0,1),
        xlab=axis.label("IS"),ylab="gamma")
      for(i in 1:4) lines(IS,10^s$out[[i]]$loggam,col=col)
    }
    iplotfun(0,"blue",add=FALSE)
    iplotfun(25,"black")
    iplotfun(40,"red")
    title(main=paste("activity coefficient (gamma) of -1,-2,-3,-4",
      "charged species at 0, 25, 40 deg C, after Alberty, 2003",
      sep="\n"),cex.main=0.95)

  } else if(which=="tpx") {

    ## Buffer + ionization: Metastablilities of
    ## thiol peroxidases from model bactera
    ## (ECOLI, BACSU mesophile; AQUAE thermophile,
    ## THIDA acidophile, BACHD alkaliphile)
    basis("CHNOS+")
    organisms <- c("ECOLI", "AQUAE", "BACSU", "BACHD", "THIDA")
    species("TPX", organisms)
    # create a buffer with our proteins in it
    mod.buffer("TPX", paste("TPX", organisms, sep="_"))
    # set up the buffered activities
    basis(c("CO2", "H2O", "NH3", "O2"), "TPX")
    a <- affinity(return.buffer=TRUE, T=50)
    basis(c("CO2", "H2O", "NH3", "O2"), as.numeric(a[1:4]))
    a <- affinity(pH=c(0, 14, 200), T=c(25, 70, 200))
    diagram(a, fill=NULL)
    title(main="Thiol peroxidases from bacteria")
    text(0.5, 66, describe.basis(thermo$basis[-6,], oneline=TRUE), adj=0)

  } else stop(paste("example '",which,"' not found",sep=""))

}
