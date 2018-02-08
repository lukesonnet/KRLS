context("Truncation works")

test_that("Truncation returns roughly same results as no truncation", {
  # This largely tests if anything goes really wrong in truncation
  # not that it is a good approximation
  
  krlso <- krls(y = mtcars$am, X = mtcars$mpg, truncate = FALSE)
  
  krlso_trunc <- krls(y = mtcars$am, X = mtcars$mpg, epsilon = 0.0000001)

  expect_true(
    ncol(krlso_trunc$U) < ncol(krlso$U)
  )
  
  expect_true(
    max(
      inference.krls2(krlso)$derivatives - 
        inference.krls2(krlso_trunc)$derivatives
    ) < 1e-6
  )
  
  expect_true(
    max(
      krlso$fitted - 
        krlso_trunc$fitted
    ) < 1e-7
  )

  # Overfitting without truncation with logistic seems to be a serious problem
  # TODO resovlve below
  # krlogo <- krls(y = mtcars$am, X = mtcars$mpg, truncate = FALSE, loss = "logistic", hyperfolds = nrow(mtcars))
  # 
  # krlogo_trunc <- krls(y = mtcars$am, X = mtcars$mpg, epsilon = 0.0000001, loss = "logistic", hyperfolds = nrow(mtcars))
  # 
  # expect_true(
  #   ncol(krlogo_trunc$U) < ncol(krlogo$U)
  # )
  # 
  # expect_true(
  #   max(
  #     inference.krls2(krlogo, vcov = FALSE)$derivatives - 
  #       inference.krls2(krlogo_trunc, vcov = FALSE)$derivatives
  #   ) < 1e-6
  # )
  # 
  # expect_equal(
  #   max(
  #     krlogo$fitted - 
  #       krlogo_trunc$fitted
  #   ) < 1e-7
  # )  
})
  