context("diagram")

test_that("diagram() handles 2D plots with different x and y resolution", {
  basis("CHNOS")
  species(c("alanine", "glycine", "serine", "methionine"))
  a <- affinity(T=c(0, 200, 6), O2=c(-90, -60, 5))
  # TODO: fix plot.line() function in diagram() so that the plot can be made
  #expect_equal(diagram(a), diagram(a, plot.it=FALSE))
})
