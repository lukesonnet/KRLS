# improve the speed in constructing kernel matrix


# slow, use rcpp version instead
# createKXY <- function(X, Y){
#   p <- ncol(X)  # bandwidth square
#   
#   n <- nrow(X)
#   m <- nrow(Y)
#   K <- matrix(1, nrow = n, ncol = m)
#   
#   for(i in 1:n){
#     for (j in 1:m) {
#       K[i,j] = exp(-sum((X[i, ]- Y[j, ])**2)/p)
#     }
#   }
#   return(K)
# }

my_krls2 <- function(y, X, lambda, b = NULL, I=NA) {
  n <- nrow(X)
  if (is.null(b)){
    b = ncol(X)
  }
  if (any(is.na(I))) {
    # full KRLS
    K <-  kern_gauss(X, b = ncol(X))
    ch <- dh <- solve(K + lambda*diag(n), y )
    yh <- K*ch
  } else {
    R <- kern_gauss_d(X, X[I, ], b= ncol(X))
    D <- R[I, ]
    RTR <- crossprod(R)
    Del <-  1e-6*diag(length(I))
    dh <- solve(RTR + lambda*D + Del, t(R)%*%y)
    ch <- R %*% solve(RTR + Del, D %*% dh)
    yh <- R %*% dh
  }
  
  e <- y - yh
  mse <- meansq(e) 
  Rsq <- 1-mse/ meansq(y)
  mod <- list(X=X, dh=dh, ch=ch, yh=yh, e=e, mse=mse, Rsq=Rsq)
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


# select m
krls_select_m <- function(y, X, lambda, m = 25, ncv = 5, press_pct = 0.25) {
  n <- nrow(X)
  if (m >= n) {
    return(1:n)
  }
  mse_tst <- matrix(0,ncv,n-m)
  I_rec <- matrix(0,ncv,m)
  for (t in 1:ncv){
    I <- sample(n,m)
    # mod <- krls2(y[I], X[I , , drop=F], lambda, I=1:m) # fits full krls
    mod <- krls2(y[I], X[I , , drop=F], lambda, I=NA)  # fits full krls
    yh_tst <- predict(mod, X[-I,,drop=F])
    mse_tst[t, ] <- (y[-I] - yh_tst)^2
    I_rec[t, ] <- I
  }
  pred_mse_vec <- rowMeans(mse_tst)
  i0 <- which.min( pred_mse_vec )
  I <- I_rec[i0,]
  press <- mse_tst[i0,]
  all_idx <- 1:n
  I <- c(I, all_idx[-I][press > press_pct*max(press)])
  I
}

# param_grid <- expand.grid(lambda = 10**( seq(-2,0.5,length.out = 10) ), m = 10)
lam_cv_krls <- function(y, X, lambda_seq, m, tr_pct=0.6, ncv=5, trim=0) {
  n <- nrow(X)
  ntr <- round(tr_pct*n)

  res <- expand.grid(lambda = lambda_seq)
  for (r in 1:nrow(res)) {
    lambda <- res[r, "lambda"]
    mse_tst <- 0
    for (t in 1:ncv){
      tr_idx <- sample(n, ntr)
      mod <- krls2(y[tr_idx], X[tr_idx , , drop=F], lambda, I=sample(ntr,m))
      yh_tst <- predict(mod,X[-tr_idx,,drop=F])
      mse_tst <- mse_tst + meansq(y[-tr_idx] - yh_tst, trim=trim) # compute trimmed mean due to possible bad column batches
    }
    res[r, "mse_tst"] <- mse_tst / ncv
  }
  
  res
}



# cv_krls <- function(y,X, param_grid = expand.grid(lambda = 10**( seq(-2,0.5,length.out = 10) ), m = 10), 
#                     tr_pct=0.8, ncv=5, trim=0) {
#   n <- nrow(X)
#   ntr <- round(tr_pct*n)
#   
#   for (r in 1:nrow(param_grid)) {
#     lambda <- param_grid[r, "lambda"]
#     m <- param_grid[r, "m"]
#     mse_tst <- 0
#     for (t in 1:ncv){
#       tr_idx <- sample(n, ntr)
#       # Xtr <- X[tr_idx , , drop=F]
#       # ytr <- y[tr_idx]
#       
#       mod <- krls2(y[tr_idx], X[tr_idx , , drop=F], lambda, m=m)
#       yh_tst <- predict(mod,X[-tr_idx,,drop=F])
#       mse_tst <- mse_tst + meansq(y[-tr_idx] - yh_tst, trim=trim) # compute trimmed mean due to possible bad column batches
#     }
#     param_grid[r, "mse_tst"] <- mse_tst / ncv
#   }
#   
#   param_grid
# }

meansq <- function(x, trim=0) mean(x^2, trim=trim)


