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
