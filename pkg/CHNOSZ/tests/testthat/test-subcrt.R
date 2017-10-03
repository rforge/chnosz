context("subcrt")

# delete the basis definition in case there is one
basis(delete=TRUE)

test_that("unbalanced reactions give a warning", {
  expect_warning(subcrt(c("glucose", "ethanol"), c(-1, 3)), "reaction was unbalanced, missing H-6O3")
})

test_that("unbalanced reactions are balanced given sufficient basis species", {
  basis("CHNOS")
  s <- subcrt(c("malic acid", "citric acid"), c(-1, 1))
  expect_equal(s$reaction$coeff, c(-1, 1, -2, -1, 1.5))
  expect_equal(s$reaction$name, c("malic acid", "citric acid", "CO2", "water", "oxygen"))
})

test_that("phase transitions of minerals give expected messages and results", {
  iacanthite <- info("acanthite", "cr2")
  #expect_message(subcrt(iacanthite), "subcrt: some points below transition temperature for acanthite cr2 \\(using NA for G\\)")
  expect_message(subcrt(iacanthite), "subcrt: some points above temperature limit for acanthite cr2 \\(using NA for G\\)")
  expect_equal(subcrt("acanthite")$out$acanthite$state, c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3))
})

test_that("heat capacity of minerals are consistent with literature", {
  # from Helgeson et al., 1978 Table 6
  # here, set P to NA so that water() is not called to get Psat
  expect_equal(subcrt("wollastonite", T=c(200, 500, 800), P=NA)$out$wollastonite$Cp, c(25.4, 28.3, 29.9), tolerance=1e-2)
})

test_that("standard Gibbs energies of reactions involving aqueous species are consistent with the literature", {

  # from Amend and Shock, 2001 [AS01] Table 3
  T <- c(2, 18, 25, 37, 45, 55, 70, 85, 100, 115, 150, 200)
  # standard Gibbs energies in kJ/mol
  E.units("J")

  # H2O(l) = H+ + OH-
  AS01.H2O <- c(78.25, 79.34, 79.89, 80.90, 81.63, 82.59, 84.13, 85.78, 87.55, 89.42, 94.22, 102.21)
  sout.H2O <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T=T)$out
  # tolerances set to lowest order of magnitute to pass
  expect_equal(sout.H2O$G/1000, AS01.H2O, tolerance=1e-4)

  # AS01 Table 4.3 Reaction A1: H2(aq) + 0.5O2(aq) = H2O(l)
  AS01.A1 <- c(-263.94, -263.45, -263.17, -262.62, -262.20, -261.63, -260.67, -259.60, -258.44, -257.18, -253.90, -248.44)
  sout.A1 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T=T)$out
  expect_equal(sout.A1$G/1000, AS01.A1, tolerance=1e-4)

  # AS01 Table 5.1 NO(g) -> NO(aq) 
  DG0.NO.g <- c(91.39, 88.04, 86.57, 84.03, 82.33, 80.20, 76.99, 73.75, 70.50, 67.23, 59.53, 48.38)
  DG0.NO.aq <- c(104.56, 102.87, 102.06, 100.58, 99.53, 98.15, 95.98, 93.67, 91.25, 88.72, 82.46, 72.77)
  sout.NO <- subcrt("NO", c("gas", "aq"), c(-1, 1), T=T)$out
  # higher tolerance, our values for NO(aq) differ slightly from AS01
  expect_equal(sout.NO$G/1000, DG0.NO.aq - DG0.NO.g, tolerance=1e-2)

  # AS01 Table 5.1 Reaction B10: NH3(aq) + 1.5O2(aq) = H+ + NO2- + H2O(l)
  AS01.B10 <- c(-268.85, -268.01, -267.50, -266.46, -265.66, -264.55, -262.66, -260.54, -258.19, -255.63, -248.81, -237.06)
  sout.B10 <- subcrt(c("NH3", "O2", "H+", "NO2-", "H2O"), "aq", c(-1, -1.5, 1, 1, 1), T=T)$out
  expect_equal(sout.B10$G/1000, AS01.B10, tolerance=1e-3)

  # AS01 Table 6.3 Reaction C7: 5S2O3-2 + H2O(l) + 4O2(aq) = 6SO4-2 + 2H+ + 4S(s)
  AS01.C7 <- c(-1695.30, -1686.90, -1682.80, -1675.30, -1670.00, -1663.10, -1652.00, -1640.30, -1628.00, -1615.20, -1583.50, -1533.00)
  s.C7 <- subcrt(c("S2O3-2", "H2O", "O2", "SO4-2", "H+", "S"), c("aq", "liq", "aq", "aq", "aq", "cr"), c(-5, -1, -4, 6, 2, 4), T=T)
  sout.C7 <- s.C7$out
  expect_equal(sout.C7$G/1000, AS01.C7, tolerance=1e-4)
  # we can also check that sulfur has expected phase transitions
  expect_equal(s.C7$state$sulfur, c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3))

  ## AS01 Table 8.3 Reaction E12: 4(2-)propanol(aq) + 3CO2(aq) + 2H2O(l) = 3CH4(aq) + 4lactic acid(aq)
  #DG0.E12 <- c(132.52, 132.26, 132.29, 132.49, 132.74, 133.15, 133.98, 135.04, 136.31, 137.79, 141.97, 149.53)
  #sout.E12 <- subcrt(c("2-propanol", "CO2", "H2O", "CH4", "lactic acid"), c(-4, -3, -2, 3, 4), T=T)$out
  ## this one is on hold until the HKF parameters of 2-propanol can be located
  #expect_equal(sout.E12$G/1000, DG0.E12, 1e-4)

  # return to our favourite units
  E.units("cal")
})

test_that("subzero degree C calculations are possible", {
  ## start with H2O
  s.H2O <- subcrt("H2O", T=c(-20.1, seq(-20, 0)), P=1)$out$water
  # we shouldn't get anything at -20.1 deg C
  expect_equal(s.H2O$G[1], NA_real_)
  # we should get something at -20 deg C
  expect_equal(floor(s.H2O$G[2]), -56001)
  # for historical reasons, an input temperature of 0 was converted to 0.01
  expect_equal(s.H2O$T[22], 0.01)
})

test_that("calculations using IAPWS-95 are possible", {
  oldwat <- water("IAPWS95")
  sb <- subcrt(c("H2O", "Na+"), T=c(-30, -20, 0, 10), P=1)$out
  # the test is not a demanding numerical comparison, more that we got numbers and no error
  expect_that(all(sb$`Na+`$G < sb$water$G), is_true())
  # clean up
  water(oldwat)
})

test_that("calculations for quartz are consistent with SUPCRT92", {
  # output from SUPCRT92 for reaction specified as "1 QUARTZ" run at 1 bar
  # (SUPCRT shows phase transition at 574.850 deg C, and does not give Cp values around the transition)
  S92_1bar <- read.table(header = TRUE, text = "
      T       G       H    S       V
    572	-214482	-209535	24.7	23.3
    573	-214507	-209517	24.7	23.3
    574	-214532	-209499	24.8	23.3
    575	-214557	-209192	25.1	23.7
    576	-214582	-209176	25.1	23.7
    577	-214607	-209159	25.2	23.7
  ")
  CHNOSZ_1bar <- subcrt("quartz", T=seq(572, 577), P=1)$out[[1]]
  expect_equal(CHNOSZ_1bar$G, S92_1bar$G, tolerance = 1e-5)
  expect_equal(CHNOSZ_1bar$H, S92_1bar$H, tolerance = 1e-5)
  expect_equal(CHNOSZ_1bar$S, S92_1bar$S, tolerance = 1e-2)
  expect_equal(CHNOSZ_1bar$V, S92_1bar$V, tolerance = 1e-2)

  # output from SUPCRT92 for reaction specified as "1 QUARTZ" run at 500 bar
  # (SUPCRT shows phase transition at 587.811 deg C)
  S92_500bar <- read.table(header = TRUE, text = "
      T       G       H    S       V
    585	-214523	-209335	24.6	23.3
    586	-214548	-209318	24.7	23.3
    587	-214573	-209301	24.7	23.3
    588	-214602	-208700	25.4	23.7
    589	-214627	-208684	25.4	23.7
    590	-214653	-208668	25.4	23.7
  ")
  CHNOSZ_500bar <- subcrt("quartz", T=seq(585, 590), P=500)$out[[1]]
  expect_equal(CHNOSZ_500bar$G, S92_500bar$G, tolerance = 1e-5)
  expect_equal(CHNOSZ_500bar$H, S92_500bar$H, tolerance = 1e-4)
  expect_equal(CHNOSZ_500bar$S, S92_500bar$S, tolerance = 1e-3)
  expect_equal(CHNOSZ_500bar$V, S92_500bar$V, tolerance = 1e-2)
})

test_that("duplicated species yield correct phase transitions", {
  # If a mineral with phase transitions is in both the basis and species lists,
  # energy()'s call to subcrt() will have duplicated species.
  # This wasn't working (produced NAs at low T) for a long time prior to 20171003.
  s1 <- subcrt("quartz", T=c(100, 1000), P=1000)
  s2 <- subcrt(rep("quartz", 2), T=c(100, 1000), P=1000)
  expect_equal(s1$out[[1]]$logK, s2$out[[1]]$logK)
  expect_equal(s1$out[[1]]$logK, s2$out[[2]]$logK)
  ## another way to test it ...
  #basis(c("quartz", "oxygen"))
  #species("quartz")
  #a <- affinity(T=c(0, 1000, 2), P=1)
  #expect_equal(as.numeric(a$values[[1]]), c(0, 0))
})

# references

# Amend, J. P. and Shock, E. L. (2001) 
#   Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. 
#   FEMS Microbiol. Rev. 25, 175--243. https://doi.org/10.1016/S0168-6445(00)00062-0

# Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) 
#   Summary and critique of the thermodynamic properties of rock-forming minerals. 
#   Am. J. Sci. 278-A, 1--229. http://www.worldcat.org/oclc/13594862

