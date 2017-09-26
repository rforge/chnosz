context("DEW")

# since this is alphabetically the first test,
# we need to load the 'thermo' object (for running tests in R CMD check)
suppressMessages(data(thermo))

test_that("density of water is calculated correctly", {
  pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
  temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
  # density from R translation of DEW macro functions
  RDensity <- calculateDensity(pressure, temperature)
  # density from DEW spreadsheet
  DEWDensity <- c(1.108200, 0.597623, 1.196591, 0.798331, 1.321050, 1.000735, 1.578116, 1.287663)
  expect_equal(RDensity, DEWDensity, tolerance=1e-6)
})

test_that("Gibbs energy of water is calculated correctly", {
  pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
  temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
  # Gibbs energies from R translation of DEW macro functions
  RGibbs <- calculateGibbsOfWater(pressure, temperature)
  # Gibbs energies from DEW spreadsheet
  DEWGibbs <- c(-56019.85419280258, -84262.028821198, -54155.004480575895, -81210.38766217149,
                -50735.122222685815, -76433.07602205424, -41823.26077175943, -65187.48113532527)
  expect_equal(RGibbs, DEWGibbs)
})

test_that("dielectric constant of water is calculated correctly", {
  pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
  temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
  # epsilon from R translation of DEW macro functions
  Repsilon <- calculateEpsilon(calculateDensity(pressure, temperature), temperature)
  # epsilon from DEW spreadsheet
  DEWepsilon <- c(65.63571, 6.10465, 72.40050, 8.97800, 82.16244, 12.13131, 103.12897, 16.97266)
  expect_equal(Repsilon, DEWepsilon, tolerance=1e-7)
})

test_that("Born coefficient Q is calculated correctly", {
  pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
  temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
  # Q from R translation of DEW macro functions
  RQ <- calculateQ(calculateDensity(pressure, temperature), temperature)
  # Q from DEW spreadsheet
  DEWQ <- c(0.32319817, 14.50286092, 0.19453478, 3.12650897,
            0.10918151, 0.87729257,  0.05068788, 0.20640645) / 1e6
  expect_equal(RQ, DEWQ)
})

test_that("g function is calculated correctly", {
  pressure <- c(1000, 1000, 5000, 5000, 10000)
  temperature <- c(100, 1000, 100, 1000, 100)
  # properties of water from DEW implementation in CHNOSZ
  ow <- water("DEW")
  # note that values returned for alpha, daldT, beta are NA
  w <- water(c("rho", "alpha", "daldT", "beta", "Psat"), T=convert(temperature, "K"), P=pressure)
  # g from CHNOSZ functions
  Rg <- gfun(w$rho/1000, temperature, pressure, w$alpha, w$daldT, w$beta)$g
  # g from R translation of DEW macro functions (not used in CHNOSZ)
  DEWg <- calculateG(pressure, temperature, w$rho/1000)
  expect_equal(Rg, DEWg)
  water(ow)
})

test_that("Gibbs energies of species are calculated correctly", {
  P <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
  T <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
  ow <- water("DEW")
  RG_HCl <- subcrt("HCl", P=P, T=T)$out[[1]]$G
  DEWG_HCl <- c(-28784.99, -58496.85, -26520.94, -55276.92, -21928.89, -50337.19, -8014.34, -36746.87)
  # TODO: expect_equal(RG_HCl, DEWG_HCl)
  RG_Cl <- subcrt("Cl-", P=P, T=T)$out[[1]]$G
  DEWG_Cl <- c(-30054.59, -22839.35, -27910.68, -28094.07, -23568.45, -27959.67, -10443.07, -18744.93)
  # TODO: expect_equal(RG_Cl, DEWG_Cl)
  water(ow)
})
