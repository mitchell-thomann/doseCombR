context("binomial models, data generation and fit")
library(doseCombR)

set.seed(1)

# monotherapy binomial data
bindat <- gendr_emax()$data
combdat <- gendr_linearcomb()$data

test_that("data generated", {
  expect_equal(sum(bindat$n), 300)
  expect_equal(sum(bindat$resp), 86)
})

# test_that("model code generated",{
#
# })
