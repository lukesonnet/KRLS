# compare the time with or without parallel

m <- 100
X <- matrix(runif(n*p), nrow = n, ncol = p)
S <- matrix(rnorm(m*n), nrow = n, ncol = m)/sqrt(m)
b <- p

system.time({
  R <- kern_gauss(X, b) %*% S
})

system.time({
  R <- apply(X, 1, function(x){
    onerow <- new_gauss_kern(matrix(X[1, ], nrow = 1), X, b)
    apply(S, 2, function(x) crossprod(onerow, x))
  })
})