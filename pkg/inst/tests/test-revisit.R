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

test_that("target calculations are working", {
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
