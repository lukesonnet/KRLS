# find the upper and lower bound of m
rm(list = ls())
# lower bound should be m that cost 30s 
source("select_m.R")

Ntru <- 100
d <- 20 # number of covariates
Xtru <- matrix(runif(Ntru*d), nrow = Ntru, ncol = d)
coeff <- rnorm(Ntru)
ftru <- function(X) {
  K_X_Xtru <- new_gauss_kern(Xtru, X, d)
  t(K_X_Xtru) %*% coeff
}
sgm <- sqrt(2)*sd(ftru(Xtru))


## training data
n <- 10000
m <- 1000
I <- sample(n, m)
X_init <- matrix(runif(n*d), nrow = n, ncol = d)
y_pure  <- ftru(X_init)
y_init <- y_pure + rnorm(n, sd = sgm)
X <- scale(X_init) # do scaling at start and end only
y <- scale(y_init)

R <- new_gauss_kern(X, X[I, ], b=d)
D <- R[I, ]

dt <- system.time( {
  lambda_opt <- select_lambda(R, D, y)
  dh_train <- train_krr(y, R, D, lambda_opt)$dh
})["elapsed"] 
dt
