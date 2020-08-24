# nystrom with adaptive sampling in pure R code


# KRLS --------------------------------------------------------------------
createKXY <- function(X, Y){
  p <- ncol(X)  # bandwidth square
  
  n <- nrow(X)
  m <- nrow(Y)
  K <- matrix(1, nrow = n, ncol = m)
  
  for(i in 1:n){
    for (j in 1:m) {
      K[i,j] = exp(-sum((X[i, ]- Y[j, ])**2)/p)
    }
  }
  return(K)
}

krls2 <- function(y, X, lambda = NULL, b = NULL) {
  n <- nrow(X)
  if (is.null(b)){
    b = ncol(X)
  }
  K <-  createKXY(X, X)
  
  if (is.null(lambda)){ #select lambda with cross validation
    lambdas = 10^(seq(-6, 2, length.out = 9))
    # assign label to each observation
    labels <- sample(5, nrow(X), replace = T)
    cv_MSE <- sapply(lambdas, function(lambda){
      MSE_vali <- sapply(1:5, function(v){
        cv_idx <- which(labels != v)
        # validation MSE
        yh_vali <- K[-cv_idx, cv_idx] %*% solve(K[cv_idx, cv_idx] + lambda*diag(length(cv_idx)), y[cv_idx])
        return(mean((yh_vali- y[-cv_idx])^2))
      })
      # average validation MSE
      return(mean(MSE_vali))
    })
    # select optimal lambda
    lambda <- lambdas[which.min(cv_MSE)]
  }
  
  ch <- solve(K + lambda*diag(n), y )
  yh <- K %*% ch
  
  e <- y - yh
  mse <- mean(e^2) 
  #Rsq <- 1-mse/ meansq(y)
  mod <- list(X=X, ch=ch, yh=yh, e=e, mse=mse, lambda = lambda)#, Rsq=Rsq)
  class(mod) <- c("krls2","list")
  mod
}

predict.krls2 <- function(mod,Xnew=NA) {
  if (any(is.na(Xnew))) {
    return(mod$yh)
  } else {
    return( createKXY(Xnew, mod$X) %*% mod$ch )
  }
}


# Adaptive selection ------------------------------------------------------

select_m <- function(y, X, m = 25, ncv = 5, z = 2) {
  n <- nrow(X)
  if (m >= n) {
    return(1:n)
  }
  mse_tst <- matrix(0,ncv,n-m)
  I_rec <- matrix(0,ncv,m)
  # get lambda based on the first batch
  I <- sample(n,m)
  mod <- krls2(y[I], X[I,,drop=F])
  lambda <- mod$lambda
  yh_tst <- predict(mod, X[-I,,drop=F])
  mse_tst[1, ] <- (y[-I] - yh_tst)^2
  I_rec[1, ] <- I
  
  for (t in 2:ncv){
    I <- sample(n,m)
    mod <- krls2(y[I], X[I , , drop=F], lambda)  # fits full krls
    yh_tst <- predict(mod, X[-I,,drop=F])
    mse_tst[t, ] <- (y[-I] - yh_tst)^2
    I_rec[t, ] <- I
  }
  
  pred_mse_vec <- rowMeans(mse_tst)
  i0 <- which.min( pred_mse_vec )
  I <- I_rec[i0,]
  press <- mse_tst[i0,]
  all_idx <- 1:n
  I <- c(I, all_idx[-I][press > mean(press) + z*sd(press)])
  
  I
}


# Nystrom -----------------------------------------------------------------

krls_nys <- function(y, X, I = NULL, lambda = NULL){
  if (is.null(I)){
    I <- select_m(y, X)
  }
  R <- createKXY(X, X[I, , drop = F])
  D <- R[I, ]
  RTR <- crossprod(R)
  Del <-  1e-6*diag(length(I))
  
  if(is.null(lambda)){ # select lambda with cv
    lambdas = 10^(seq(-6, 2, length.out = 9))
    labels <- sample(5, nrow(X), replace = T)
    cv_MSE <- sapply(lambdas, function(lambda){
      MSE_vali <- sapply(1:5, function(v){
        cv_idx <- which(labels != v)
        R_vali <- R[cv_idx, ]
        RTR_vali <- crossprod(R_vali)
        dh_vali <- solve(RTR_vali + lambda*D + Del, t(R_vali)%*%y[cv_idx])
        ch_vali <- R_vali %*% solve(RTR_vali + Del, D %*% dh_vali)
        yh_vali <- createKXY(X[-cv_idx, ], X[cv_idx, ]) %*% ch_vali
        return(mean((yh_vali- y[-cv_idx])^2))
      })
      # average validation MSE
      return(mean(MSE_vali))
    })
    # select optimal lambda
    lambda <- lambdas[which.min(cv_MSE)]
  }
  
  dh <- solve(RTR + lambda*D + Del, t(R)%*%y)
  ch <- R %*% solve(RTR + Del, D %*% dh)
  yh <- R %*% dh
  
  e <- y - yh
  mse <- mean(e^2) 
  #Rsq <- 1-mse/ meansq(y)
  mod <- list(X=X, ch=ch, dh = dh, I =I, yh=yh, e=e, mse=mse, lambda = lambda)#, Rsq=Rsq)
  class(mod) <- c("krls2","list")
  mod
}

predict.krls_nys <- function(mod,Xnew=NA) {
  if (any(is.na(Xnew))) {
    return(mod$yh)
  } else {
    return( createKXY(Xnew, mod$X) %*% mod$ch )
  }
}
