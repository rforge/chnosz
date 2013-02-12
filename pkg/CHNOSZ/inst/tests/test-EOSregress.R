context("EOSregress")

test_that("EOSvar stops with unknown variables", {
  expect_error(EOSvar("TX", T=25, P=1), "can't find a variable named TX")
})

