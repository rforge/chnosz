context("revisit")

# initial setup
suppressMessages({
  basis("CHNOS")
  species(c("leucine", "glycine", "glutamic acid"))
  a <- affinity(O2=c(-75, -65))
  e <- equilibrate(a)
})

test_that("inconsistent arguments produce an error", {
  expect_error(revisit(e, "logref"), "specified target 'logref' not available") 
  expect_error(revisit(e, "logact"), "for 'logact' target, loga.ref must be supplied") 
})

test_that("target calculations give expected results", {
  r.cv <- revisit(e, "CV", plot.it=FALSE)
  r.sd <- revisit(e, "sd", plot.it=FALSE)
  r.shannon <- revisit(e, "shannon", plot.it=FALSE)
  r.qqr <- revisit(e, "qqr", plot.it=FALSE)
  # the tests will alert us to significant numerical changes
  # but so far haven't been independently verified
  expect_equal(r.cv$extval, 0.30576, tol=1e-5) 
  expect_equal(r.sd$extval, 0.000284694, tol=1e-5) 
  expect_equal(r.shannon$extval, 1.066651, tol=1e-5)
  expect_equal(r.qqr$extval, 0.999783, tol=1e-5)
  # where logarithm of activity of the 3rd species (glutamic acid) maximizes
  r.logact <- revisit(e, "logact", 3, plot.it=FALSE)
  expect_equal(r.logact$ix, 71)
})

test_that("DGxf target gives zero at equilibrium and >0 not at equilibrium", {
  # let's use n-alkanes
  basis(c("CH4", "H2"), c("gas", "gas"))
  species(c("methane", "ethane", "propane", "n-butane"), "liq")
  # calculate equilibrium distribution over a range of logaH2
  a <- affinity(H2=c(-10, -5, 101), exceed.Ttr=TRUE)
  e <- equilibrate(a)
  # the reference equilibrium distribution at logfH2 = -7.5
  loga.ref <- list2array(e$loga.equil)[51, ]
  # calculate the DGxf/RT relative to the reference distribution
  # (we use swap12 so that the reference distribution is the initial one)
  r <- revisit(e, "DGxf", loga.ref=loga.ref, DGxf.swap12=TRUE, plot.it=FALSE)
  # we should find a minimum of zero at logfH2 = -7.5
  expect_equal(min(r$H), 0)
  expect_equal(a$vals[[1]][which.min(r$H)], -7.5)
  # we can also call the DGxf function directly
  # again with reference distribution as the initial one, 
  # but this time the Astar is kept constant (for logfH2 = -7.5)
  Astar <- list2array(e$Astar)[51, ]
  DGxf.out <- numeric()
  for(i in 1:length(a$vals[[1]])) {
    loga.equil <- list2array(e$loga.equil)[i, ]
    DGxf.out <- c(DGxf.out, DGxf(loga.ref, loga.equil, Astar))
  }
  # we should find a minimum of zero at logfH2 = -7.5
  expect_equal(min(DGxf.out), 0)
  expect_equal(a$vals[[1]][which.min(DGxf.out)], -7.5)
})
