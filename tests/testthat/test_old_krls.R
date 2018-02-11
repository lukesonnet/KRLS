# This file is .Rbuildignore'd as it needs the old KRLS package

context("Match KRLS 1.0.0")

n <- 20
X <- cbind(rnorm(n), rbinom(n, 1, 0.5))
y <- rnorm(20)

test_that("", {
  
  # Default sigma (b) is different
  oldout <- KRLS::krls(X = X, y = y, sigma = 2*ncol(X), print.level = 3)
  newout <- KRLS2::inference.krls2(KRLS2::krls(X = X, y = y, printlevel = 3))
  
  expect_equivalent(
    oldout$lambda,
    newout$lambda
  )
  
  # Lambda and coeffs are different because we use eigtrunc
  expect_equal(
    oldout$coeffs,
    newout$coeffs
  )
  
  expect_equal(
    oldout$derivatives,
    newout$derivatives
  )
  
  expect_equivalent(
    oldout$avgderivatives,
    newout$avgderivatives
  )
  
  expect_equivalent(
    oldout$var.avgderivatives,
    newout$var.avgderivatives * c(1, 2) # old variance of FDs had been multiplied by 2
  )
  
})
  