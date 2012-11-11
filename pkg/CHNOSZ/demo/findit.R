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
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
d <- diagram(e, col=col, names=NULL)
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
  "qqr",iprotein=ip,niter=8,normalize=TRUE)
title(main="searching 4-D chemical activity space")
# plot 5: q-q plot at the final loga O2, H2O, CO2, NH3
# higher correlation coefficient than plot 3
a <- affinity(iprotein=ip)
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
qqr5 <- revisit(e, "qqr",pch=pch)$H
legend("topleft",pch=c(1,2,4),legend=myg)
db <- describe.basis(ibasis=c(5, 2, 1, 3))
legend("bottomright", legend=db)
# plot 5: trajectory of O2, H2O, CO2, NH3, and the
# q-q correlation coefficient in the search
#plot(f1,mar=c(2,5,1,1),mgp=c(4,1,0))
par(opar)
