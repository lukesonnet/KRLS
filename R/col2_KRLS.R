#' @export

# KRLS with column sampling 1
# the ordinary column sampling method from literature

col2_KRLS <- function(X, y, b, m,
                      lambdas =  c(10^(seq(-6, 2, length.out = 20))),
                      folds = 5,
                      n_av, control){
  require(expm)
  n <- nrow(X)
  I <- sample(n, m)
  
  # parameters in column sampling version 1
  R <- switch(control$whichkernel,
              gaussian=new_gauss_kern(X, X[I, ], b),
              linear=tcrossprod(X, X[I, ]),
              gausslinear=.5*(new_gauss_kern(X, X[I, ], b)+tcrossprod(X, X[I, ])),
              poly2=(tcrossprod(X, X[I, ])+1)^2,
              poly3=(tcrossprod(X, X[I, ])+1)^3,
              poly4=(tcrossprod(X, X[I, ])+1)^4,
              stop("No valid Kernel specified") )
  D <- diag(m)
  
  # assess the quality of the selected columns
  # calculating the average R square
  K_y <-switch(control$whichkernel,
               gaussian=new_gauss_kern(X, X[sample(n, n_av), ], b),
               linear=tcrossprod(X, X[sample(n, n_av), ]),
               gausslinear=.5*(new_gauss_kern(X, X[sample(n, n_av), ], b)+tcrossprod(X, X[sample(n, n_av), ])),
               poly2=(tcrossprod(X, X[sample(n, n_av), ])+1)^2,
               poly3=(tcrossprod(X, X[sample(n, n_av), ])+1)^3,
               poly4=(tcrossprod(X, X[sample(n, n_av), ])+1)^4,
               stop("No valid Kernel specified") )
  Rsq <- ave_Rsq(R, K_y)
  
  # cross validation
  lables <- sample(folds, n, replace = T)
  cv_MSE <- rep(0, length(lambdas))
  
  for (i in 1:length(lambdas)){
    lambda <- lambdas[i]
    MSE_vali <- rep(0, folds)
    for (v in 1:folds){
      # training indices in cross validation
      cv_idx <- which(lables != v)
      
      # validation MSE
      dh_vali <- train_krr(y[cv_idx], R[cv_idx, ], D, lambda)$dh
      MSE_vali[v] <- test_krr(dh_vali, y[-cv_idx], R[-cv_idx, ])$test_MSE
    }
    
    # cross validation MSE 
    cv_MSE[i] <- mean(MSE_vali)
  }
  
  # select optimal lambda
  lambda_opt <- lambdas[which.min(cv_MSE)]
  train_result <- train_krr(y, R, D, lambda_opt)
  
  return(list(coeffs = train_result$dh,
              fitted = train_result$yfitted,
              I = I,      
              Rsq = Rsq, 
              lambda = lambda_opt))
}
