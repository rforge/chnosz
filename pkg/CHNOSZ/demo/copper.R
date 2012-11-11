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
