context("revisit")

# initial setup
suppressMessages({
  basis("CHNOS")
  species(c("leucine", "glycine", "glutamic acid"))
  a <- affinity(O2=c(-75, -65))
  d <- diagram(a, plot.it=FALSE)
})

test_that("inconsistent arguments produce an error", {
  expect_error(revisit(d, "logref"), "specified target 'logref' not available") 
  expect_error(revisit(d, "logact"), "for 'logact' target, loga.ref must be supplied") 
})

test_that("target calculations give expected results", {
  r.cv <- revisit(d, "CV", plot.it=FALSE)
  r.sd <- revisit(d, "sd", plot.it=FALSE)
  r.shannon <- revisit(d, "shannon", plot.it=FALSE)
  r.qqr <- revisit(d, "qqr", plot.it=FALSE)
  # the tests will alert us to significant numerical changes
  # but so far haven't been independently verified
  expect_equal(r.cv$extval, 0.30576, tol=1e-5) 
  expect_equal(r.sd$extval, 0.000284694, tol=1e-5) 
  expect_equal(r.shannon$extval, 1.066651, tol=1e-5)
  expect_equal(r.qqr$extval, 0.999783, tol=1e-5)
  # where logarithm of activity of the 3rd species (glutamic acid) maximizes
  r.logact <- revisit(d, "logact", 3, plot.it=FALSE)
  expect_equal(r.logact$ix, 71)
})

test_that("DGT target gives zero at equilibrium and >0 not at equilibrium", {
  # let's use n-alkanes
  basis(c("CH4", "H2"), c("gas", "gas"))
  species(c("methane", "ethane", "propane", "n-butane"), "liq")
  # calculate equilibrium distribution over a range of logaH2
  a <- affinity(H2=c(-10, -5, 101), exceed.Ttr=TRUE)
  d <- diagram(a, plot.it=FALSE)
  # the reference equilibrium distribution at logfH2 = -7.5
  loga.ref <- list2array(d$logact)[51, ]
  # calculate the DGT/RT relative to the reference distribution
  # (we use swap12 so that the reference distribution is the initial one)
  r <- revisit(d, "DGT", loga.ref=loga.ref, DGT.swap12=TRUE, plot.it=FALSE)
  # we should find a minimum of zero at logfH2 = -7.5
  expect_equal(min(r$H), 0)
  expect_equal(a$xvals[which.min(r$H)], -7.5)
  # we can also call the DGT function directly
  # again with reference distribution as the initial one, 
  # but this time the Astar is kept constant (for logfH2 = -7.5)
  Astar <- list2array(d$Astar)[51, ]
  DGT.out <- numeric()
  for(i in 1:length(a$xvals)) {
    loga.equil <- list2array(d$logact)[i, ]
    DGT.out <- c(DGT.out, DGT(loga.ref, loga.equil, Astar))
  }
  # we should find a minimum of zero at logfH2 = -7.5
  expect_equal(min(DGT.out), 0)
  expect_equal(a$xvals[which.min(DGT.out)], -7.5)
})
