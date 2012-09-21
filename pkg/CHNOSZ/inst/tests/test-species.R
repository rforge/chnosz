context("species")

# clear out any previous basis definition or database alterations
suppressPackageStartupMessages(data(thermo))

test_that("species not contained by basis cause errors", {
  expect_error(species("H2O"), "basis species are not defined")
  expect_error(species.basis("H2O"), "basis species are not defined")
  basis("CHNOS")
  expect_error(species.basis("U"), "element\\(s\\) not in the basis\\: U")
  expect_error(species("fayalite"), "element\\(s\\) not in the basis\\: Fe Si")
})

test_that("for one or more species, species.basis() keeps track of zeroes and puts elements in order of thermo$basis", {
  basis("CHNOS")
  test0 <- count.elements("OHN0")
  test1 <- count.elements("HN0O")
  expect_equal(species.basis(test0), species.basis(test1))
  # we can send multiple species to species.basis() but the argument has to be constructed correctly
  expect_equal(unique(as.numeric(species.basis(makeup(c("C", "CCN"))))), 0)
  expect_equal(species.basis(makeup(c("C", "CCN"), count.zero=TRUE))[2, , drop=FALSE], species.basis(makeup("CCN")))
})

test_that("deleting nonexistent species causes error or warning", {
  expect_error(species("CO2", delete=TRUE), "nonexistent species definition")
  species("H2O")
  expect_warning(species("CO2", delete=TRUE), "not present, so can not be deleted")
  expect_is(species("water", delete=TRUE), "NULL")
})
