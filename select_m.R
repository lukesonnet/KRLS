findMSE <- function(i, K_max, y, y_cv, lambda, I_max ,idx, K_cv){
  K_i <- K_max[1:idx[i], 1:idx[i]]
  y_i <- y[I_max[1:idx[i]]]
  ch <- solve(K_i + lambda*diag(nrow(K_i)), y_i)
  yh <- K_cv[, 1:idx[i]] %*% ch 
  MSE <- mean((yh - y_cv)^2)
  return(MSE)
}

# nystrom KRLS, tuning both m and lambda
select_m <- function(X, y, m_0 = 50, m_max= 1000,  
              lambdas = 10^(seq(-6, 2, length.out = 10))){
  n <- nrow(X)
  d <- ncol(X)
  if (m_max >= n){
    m_max  = round(n* 0.8)
  }
  
  # get all the possible index at one time
  I_max <- sample(n, m_max) 
  idx <- unique(c(seq(m_0, m_max, by= m_0),m_max))
  # construct the largest kernal 
  K_max <- kern_gauss(X[I_max, ], b = d)#do not do this! waste time!
  K_cv <- new_gauss_kern(X[-I_max, ], X[I_max, ], b = d)
  y_cv <- y[-I_max]
  
   # find lambda using the first set of the index
  K_1 <- K_max[1:m_0, 1:m_0]
  y_1 <- y[I_max[1:m_0]]
  MSE_1 <- sapply(lambdas, function(x){
    ch <- solve(K_1 + x*diag(m_0), y_1)
    yh <- K_cv[, 1:m_0] %*% ch 
    mean((yh - y_cv)^2)
  })
  lambda <- lambdas[which.min(MSE_1)]
  MSE <- min(MSE_1)
  
  # initialize the first 2 values of the MSE vector
  MSE[2] <- findMSE(2, K_max, y, y_cv, lambda, I_max ,idx, K_cv)
  
  i <- 2
  # select m when the drop in MSE is not so sifnificant
  while ((MSE[i-1] - MSE[i])/ MSE[i-1] > 0.001 & i <length(idx)) {
    i <- i+1
    MSE[i] <- findMSE(i, K_max, y, y_cv, lambda, I_max ,idx, K_cv)
  }
  m_opt <- idx[i-1] # select the optimal m
  I <- I_max[1:m_opt]
  return(I)
}

## select lambda 
select_lambda <- function(R, D, y, folds = 5,
                          lambdas = 10^(seq(-6, 2, length.out = 10))){
  # assign label to each observation
  labels <- sample(folds, nrow(R), replace = T)
  cv_MSE <- sapply(lambdas, function(lambda){
    MSE_vali <- sapply(1:folds, function(v){
      cv_idx <- which(labels != v)
      # validation MSE
      dh_vali <- train_krr(y[cv_idx], R[cv_idx, ], D, lambda)$dh
      yh_vali <- R[-cv_idx, ] %*% dh_vali
      return(mean((yh_vali- y[-cv_idx])^2))
    })
  # average validation MSE
  return(mean(MSE_vali))
  })
  # select optimal lambda
  lambda_opt <- lambdas[which.min(cv_MSE)]
  return(lambda_opt)
}

# select lambda based on ch
# 1/6: still not very sure about this method
# essentially, we would get the whole n*n kernel matrix
select_lambda_2 <- function(R, D, y, X, folds = 5,
                          lambdas = 10^(seq(-6, 2, length.out = 10))){
  # assign label to each observation
  labels <- sample(folds, nrow(R), replace = T)
  cv_MSE <- sapply(lambdas, function(lambda){
    MSE_vali <- sapply(1:folds, function(v){
      cv_idx <- which(labels != v)
      # validation MSE
      R_vali <- R[cv_idx, ]
      dh_vali <- train_krr(y[cv_idx], R[cv_idx, ], D, lambda)$dh
      ch_vali <- R_vali %*% solve(crossprod(R_vali) + diag(length(I)), D %*% dh_vali)
      K_vali <- new_gauss_kern(X[-cv_idx, ], X[cv_idx, ], b=d)
      yh_vali <- K_vali %*% ch_vali
      return(mean((yh_vali- y[-cv_idx])^2))
    })
    # average validation MSE
    return(mean(MSE_vali))
  })
  # select optimal lambda
  lambda_opt <- lambdas[which.min(cv_MSE)]
  return(lambda_opt)
}



