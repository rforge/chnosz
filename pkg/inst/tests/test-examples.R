context("examples")

test_that("longex() gives an error for unknown example", {
  expect_error(longex("Xxx"), "example 'xxx' not found")
})

