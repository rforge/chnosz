# test-berman.R 20171001
context("berman")

# calculate properties for all available minerals at Tr,Pr
file <- system.file("extdata/Berman/Ber88.csv", package="CHNOSZ")
dat <- read.csv(file, as.is=TRUE)
mineral <- unique(dat$name)
prop_Berman <- NULL

test_that("properties of all minerals are computed without warnings", {
  # running this without error means that:
  # - formulas for the minerals are found in thermo$obigt
  # - there are no warnings for minerals with GfPrTr(calc) >= 1000 J/cal different from GfPrTr(table)
  expect_silent(properties <- lapply(mineral, berman, check.G=TRUE))
  # save the results so we can use them in the next tests
  assign("prop_Berman", properties, inherits=TRUE)
  
  ## - warnings are produced for 7 minerals with GfPrTr(calc) >= 1000 J/cal different from GfPrTr(table)
  #expect_warning(Berman <- lapply(mineral, berman),
  #               "annite|dawsonite|dravite|fluortremolite|greenalite|siderite|Na-Fe-saponite.3W", all=TRUE)
})

# assemble a data frame for Berman properties
prop_Berman <- do.call(rbind, prop_Berman)
# find the mineral data using Helgeson formulation
icr <- suppressMessages(info(mineral, "cr"))
icr1 <- suppressMessages(info(mineral, "cr1"))
# if there's no match for "cr", try to get "cr1" (lowest-temperature phase)
icr[is.na(icr)] <- icr1[is.na(icr)]
# all of these except rutile (Robie et al., 1979) reference Helgeson et al., 1978
# NOTE: with check.it = TRUE (the default), this calculates Cp from the tabulated Maier-Kelley parameters
prop_Helgeson <- suppressMessages(info(icr))

# now we can compare Berman and Helgeson G, H, S, Cp, V
# minerals with missing properties are not matched here
# (i.e. fluorphlogopite, fluortremolite, glaucophane, and pyrope: no G and H in prop_Helgeson data)

test_that("Berman and Helgeson properties have large differences for few minerals", {
  # which minerals differ in DGf by more than 4 kcal/mol?
  idiffG <- which(abs(prop_Berman$G - prop_Helgeson$G) > 4000)
  expect_match(mineral[idiffG],
               "anthophyllite|antigorite|Ca-Al-pyroxene|lawsonite|margarite|merwinite")
  ## we find 9 of them, as follow:
  #expect_match(mineral[idiffG],
  #             "anthophyllite|antigorite|Ca-Al-pyroxene|cordierite,dry|cordierite,hydrous|lawsonite|margarite|merwinite|paragonite")

  # which minerals differ in DHf by more than 4 kcal/mol?
  idiffH <- which(abs(prop_Berman$H - prop_Helgeson$H) > 4000)
  # we get the above, plus phlogopite and clinozoisite:
  expect_match(mineral[idiffH],
               "anthophyllite|antigorite|Ca-Al-pyroxene|lawsonite|margarite|merwinite|phlogopite|clinozoisite")

  # which minerals differ in S by more than 4 cal/K/mol?
  idiffS <- which(abs(prop_Berman$S - prop_Helgeson$S) > 4)
  expect_match(mineral[idiffS], "albite|almandine")
  #expect_match(mineral[idiffS], "albite|almandine|annite|cordierite,hydrous|fluortremolite")

  # which minerals differ in Cp by more than 4 cal/K/mol?
  idiffCp <- which(abs(prop_Berman$Cp - prop_Helgeson$Cp) > 4)
  expect_match(mineral[idiffCp], "antigorite|cristobalite,beta|K-feldspar")
  #expect_match(mineral[idiffCp],
  #             "antigorite|cordierite,hydrous|cristobalite,beta|fluortremolite|glaucophane|greenalite|K-feldspar")

  # which minerals differ in V by more than 1 cm^3/mol?
  idiffV <- which(abs(prop_Berman$V - prop_Helgeson$V) > 1)
  expect_match(mineral[idiffV], "anthophyllite|antigorite|chrysotile|merwinite")
  #expect_match(mineral[idiffV],
  #             "anthophyllite|antigorite|chrysotile|cordierite,hydrous|glaucophane|greenalite|merwinite")
})
