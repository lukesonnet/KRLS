library(testthat)
library(KRLS2)
context("Test different input configurations work\n")

n <- 50
betas <- rnorm(2)
X <- cbind(rbinom(n, 1, 0.4), rnorm(n))
ypure <- X %*% betas 
y <- ypure + 2*sd(ypure)*rnorm(n)
ybin <- rbinom(n, 1, 1/(1+exp(-y)))
weights <- runif(n)
clusters <- sample(1:5, size = n, replace = T)
lambdagrid <- c(0.01, 0.1, 1, 2)

## grid of possibilities
test_grid <- expand.grid(
  loss = c('leastsquares', 'logistic'),
  lambdasearch = c('optim', 'fixed', 'grid'),
  whichkernel = c('gaussian', 'poly2'),
  truncation = c(T, F),
  robust = c('reg', 'rob', 'cl'),
  weighted = c(T, F)
)

for(i in 1:nrow(test_grid)) {
  cl <- NULL
  rob <- F
  w <- NULL
  lambda <- NULL
  epsilon <- NULL
  
  if (test_grid$truncation[i]) {
    epsilon <- 0.01
  }

  if (test_grid$loss[i] == 'leastsquares') {
    ytemp <- y    
  } else {
    ytemp <- ybin
  }
  
  if (test_grid$robust[i] == 'cl') {
    cl <- clusters
  } else if (test_grid$robust[i] == 'rob') {
    rob <- T
  }
  
  if (test_grid$weighted[i]) {
    w <- weights
  }
  
  if (test_grid$lambdasearch[i] == 'fixed') {
    lambda <- 0.9
  } else if (test_grid$lambdasearch[i] == 'grid') {
    lambda <- lambdagrid    
  }

  cat('\n')
  print(test_grid[i,])
  # don't support these, should throw error when fitting KR*
  if ((test_grid$weighted[i] & test_grid$loss[i] == 'leastsquares' & !test_grid$truncation[i])) {
    print('expect fit error')
    expecting_fit_error <- 'weighted' #this means we will expect an error with one of these things in it (regex)
  } else {
    print('expect no fit error')
    expecting_fit_error <- NA #this means we will not expect an error
  }
  
  # don't support these, should throw error when getting pwmfx
  if (test_grid$whichkernel[i] != 'gaussian') {
    print('expect inference error')
    expecting_inf_error <- 'Robust|kernel' #this means we will expect an error with one of these things in it (regex)
  } else {
    print('expect no inference error')
    expecting_inf_error <- NA #this means we will not expect an error
  }
  
  expect_error({
    # Comment this line to print warnings and messages
    capture.output({suppressWarnings({suppressMessages({
      kout <- KRLS2::krls(
        X=X,
        y=ytemp,
        w=w,
        loss=test_grid$loss[i],
        lambda=lambda,
        whichkernel=test_grid$whichkernel[i],
        epsilon = epsilon
      )
    # Also comment out line below to print warnings and messages
    })})})
    },
    expecting_fit_error
    
  )
  
  # only check summary if we expect krls to work under these conditions
  if (is.na(expecting_fit_error)) { # ie. if we are not expecting a fit error
    expect_error({
      # Comment this line to print warnings and messages
      capture.output({suppressWarnings({suppressMessages({
        summary(
          kout,
          robust = rob,
          clusters = cl
        )
      # Also comment out line below to print warnings and messages
      })})})
    },
    expecting_inf_error
    )
  }
}
