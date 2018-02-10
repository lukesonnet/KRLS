context("KRLogit gradient/hessian correct")

# Note that this file is ignored by buildignore/CRAN to avoid
# adding numDeriv to `suggests` in DESCRIPTION

test_that("", {
  library(numDeriv)

  set.seed(42)

  # Fit model
  krlog_out <- krls(
    y = mtcars$am,
    X = mtcars[, c("mpg", "wt")],
    loss = "logistic",
    lambda = 0.2,
    epsilon = 0.01,
    returnopt = TRUE
  )
  
  par <- rnorm(ncol(krlog_out$U) + 1) # +1 for beta0
  U <- krlog_out$U
  D <- krlog_out$D
  y <- krlog_out$y
  w <- krlog_out$w
  lambda <- 0.1
  
  # -----
  # Check gradient
  # -----
  
  # numerical deriv
  num_grad <- grad(
    krlogit_fn_trunc,
    x = par,
    U = U,
    D = D,
    y = y,
    w = w,
    lambda = lambda
  )
  
  # our grad fn
  an_grad <- krlogit_gr_trunc(
    par = par,
    U = U,
    D = D,
    y = y,
    w = w,
    lambda = lambda
  )
  
  # looks good!
  expect_equivalent(
    num_grad,
    an_grad
  )
  
  # -----
  # Check hessian
  # -----
  
  # numerical hess
  num_hess <- hessian(
    krlogit_fn_trunc,
    x = par,
    U = U,
    D = D,
    y = y,
    w = w,
    lambda = lambda
  )
  
  # our hess
  an_hess <- solve(krlogit_hess_trunc_inv(
    par = par,
    U = U,
    D = D,
    y = y,
    w = w,
    lambda = lambda
  ))
  
  
  expect_equivalent(
    num_hess,
    an_hess
  )
  
})
