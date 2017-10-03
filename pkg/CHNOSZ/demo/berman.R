# CHNOSZ/demo/berman.R  20171003
# make some mineral activity diagrams using Berman (1988) and related data

# using the Helgeson data
# set up basis species
basis(c("K+", "Al+3", "quartz", "H2O", "O2", "H+"))
# use pH = 0 so that aK+ = aK+/aH+
basis("pH", 0)
# load the species
species(c("K-feldspar", "muscovite", "kaolinite", "pyrophyllite", "andalusite"), "cr")
# calculate affinities in aK+ - temperature space
a <- affinity(`K+`=c(0, 6), T=c(25, 650), P=1000)
# note that we go just past the quartz transition, but it has no effect on the diagram
diagram(a)

# now using the Berman data
basis("SiO2", "cr_Berman")
# it might be good to check that we have Berman's quartz and not coesite or some other SiO2 phase
info(basis()$ispecies[3])
# remove the Helgeson minerals
species(delete=TRUE)
# load the Berman minerals
species(c("K-feldspar", "muscovite", "kaolinite", "pyrophyllite", "andalusite"), "cr_Berman")
# calculate affinities in aK+ - temperature space
a <- affinity(`K+`=c(0, 6), T=c(25, 650), P=1000)
diagram(a, add=TRUE, names="", col="blue", lwd=2)

legend("topleft", lty=c(1, 1), lwd=c(1, 2), col=c("black", "blue"), legend=c("Helgeson et al., 1978", "Berman, 1988"))
