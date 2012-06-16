context("affinity")

test_that("affinity() in 3D returns values consistent with manual calculation", {
  # our "manual" calculation will be for H2(aq) + 0.5O2(aq) = H2O(l)
  # the equilibrium constants at 25 and 100 degrees C
  # (the logK are tested against literature values in test-subcrt.R)
  logK.25 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T=25)$out$logK
  logK.100 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T=100)$out$logK
  # the value of A/2.303RT at 25 degrees and logaH2=-10, logaO2=-10 and logaH2O=0
  A.2303RT.25.10.10 <- logK.25 - ( (-1)*(-10) + (-0.5)*(-10) )
  # the value of A/2.303RT at 100 degrees and logaH2=-5, logaO2=-10 and logaH2O=0
  A.2303RT.100.5.10 <- logK.100 - ( (-1)*(-5) + (-0.5)*(-10) )
  # set up basis and species
  basis(c("H2", "O2"), "aq")
  species("H2O")
  # we will run affinity() in 3D
  # T = 0, 25, 50, 75, 100, 125 degrees
  # log_a(H2) = -20, -15, -10, -5, 0
  # log_a(O2) = -20, -15, -10, -5, 0
  # first test: the dimensions are correct
  a.logK <- affinity(T=c(0, 125, 6), H2=c(-20, 0, 5), O2=c(-20, 0, 5), property="logK")
  expect_equal(dim(a.logK$values[[1]]), c(6, 5, 5))
  # second and third tests: the logK values used by affinity() are correct
  expect_equal(a.logK$values[[1]][2, 3, 3], logK.25)
  expect_equal(a.logK$values[[1]][5, 4, 3], logK.100)
  # fourth and fifth tests: the A/2.303RT values returned by affinity() are correct
  a.A <- affinity(T=c(0, 125, 6), H2=c(-20, 0, 5), O2=c(-20, 0, 5))
  expect_equal(a.A$values[[1]][2, 3, 3], A.2303RT.25.10.10)
  expect_equal(a.A$values[[1]][5, 4, 3], A.2303RT.100.5.10)
})

test_that("'iprotein' gives consistent results on a transect", {
  # from Dick and Shock, 2011, values of A/2.303RT for the per-residue
  # formation reactions of overall model proteins at five sampling sites
  # at Bison Pool, with different temperature, pH and log_a(H2)
  # these are the maximum values for each site from Table 5 in the paper
  A.2303RT_ref <- c(-18.720, -27.894, -35.276, -36.657, -41.888)
  # the measured temperatures and pHs
  T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
  pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)
  # Eq. 24 of the paper
  H2 <- -11+T*3/40
  basis(c("HCO3-", "H2O", "NH3", "HS-", "H2", "H+"),
    "aq", c(-3, 0, -4, -7, 999, 999))
  sites <- c("N", "S", "R", "Q", "P")
  proteins <- paste("overall_bison", sites, sep="")
  ip <- iprotein(proteins)
  # to reproduce, we need use the "old" parameters for [Met]
  add.obigt()
  a <- affinity(T=T, pH=pH, H2=H2, iprotein=ip)
  # divide A/2.303RT by protein length
  pl <- protein.length(ip)
  A.2303RT <- t(sapply(a$values, c)) / pl
  # find the maximum for each site
  A.2303RT_max <- apply(A.2303RT, 2, max)
  # we're off a bit in the second decimal ... 
  # maybe becuase of rounding of the aa composition?
  expect_equal(A.2303RT_max, A.2303RT_ref, 1e-3)
  # todo: add comparison with results from loading proteins via species()
})

test_that("affinity() for proteins (with/without 'iprotein') returns same value as in previous package versions", {
  # our test case is CSG_HALJP because it has no methionine
  # (aqueous [Met] was updated in 0.9-8)
  # these values were calculated using versions 0.6, 0.8 and 0.9-7
  # (25 degrees C, 1 bar, basis species "CHNOS" or "CHNOS+")
  A.2303RT.nonionized <- -3795.297
  A.2303RT.ionized <- -3075.222
  # first for nonionized protein
  basis("CHNOS")
  # try it with iprotein
  ip <- iprotein("CSG_HALJP")
  expect_equal(affinity(iprotein=ip)$values[[1]][1], A.2303RT.nonionized, 1e-6)
  # then with the protein loaded as a species
  species("CSG_HALJP")
  expect_equal(affinity()$values[[1]][1], A.2303RT.nonionized, 1e-6)
  # now for ionized protein
  basis("CHNOS+")
  expect_equal(affinity(iprotein=ip)$values[[1]][1], A.2303RT.ionized, 1e-6)
  species("CSG_HALJP")
  expect_equal(affinity()$values[[1]][1], A.2303RT.ionized, 1e-6)
})
