
# test code
n <- 2000
Ntru <- 100
p <- 20 # number of covariates
b <- sqrt(p)
Xtru <- matrix(runif(Ntru*p), nrow = Ntru, ncol = p)
coeff <- rnorm(Ntru)
ftru <- function(X) {
  K_X_Xtru <- new_gauss_kern(Xtru, X, b)
  t(K_X_Xtru) %*% coeff
}

sgm <- sqrt(2)*sd(ftru(Xtru))
X <- matrix(runif(n*p), nrow = n, ncol = p)
y_pure  <- ftru(X)
y <- y_pure + rnorm(n, sd = sgm)

nys_result <- krls(X = X, y = y, aprxmethod = "col2")
nys_result$lambda
mean((nys_result$fitted - nys_result$y)^2)

predict_result <- predict.krls2(nys_result, newdata = matrix(runif(100*p), nrow = 100, ncol = p))
mean((predict_result$fit - y[sample(n, 100)])^2)

krls_result <- krls(X = X, y = y)
mean((krls_result$fitted - krls_result$y)^2)


