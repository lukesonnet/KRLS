# this is the very complex version
# new cross validation method to select m 
n <- 50000
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

# the selected m should be in the position (i-1)
findMSE <- function(i, X, y, lambdas, I_max ,idx){
  I <- I_max[1:idx[i]]
  # need to scale before puting them in the model
  X_i <- scale(X[I, ])
  y_i <- scale(y[I])
  K_i <- kern_gauss(X_i, b = d)
  # scale the rest of the data with the training scales
  X_cv <- scale(X[-I_max, ], center = attributes(X_i)$`scaled:center`, scale = attributes(X_i)$`scaled:scale`)  
  y_cv <- scale(y[-I_max, ], center = attributes(y_i)$`scaled:center`, scale = attributes(y_i)$`scaled:scale`)
  K_cv <- new_gauss_kern(X_cv, X_i, b = d)
  # find the optimal lambda with the least MSE (best those m points can do)
  MSE <- min(sapply(lambdas, function(x){
    ch <- solve(K_i + x*diag(length(I)), y_i)
    yh <- (K_cv %*% ch + attributes(y_i)$`scaled:center`) * attributes(y_i)$`scaled:scale`
    sqrt(mean((yh - y[-I_max])^2)) #sqrt(MSE)
    # yh <- K_cv %*% ch 
    # mean((yh - y+_cv)^2)
  }))
  return(MSE)
}

# more complex models
# n <- 1e4
# 
# p <- 5
# fried1 <- sim_friedman(p = 5, ng = 10)

# real data
X <- as.matrix(dat2.train[, -ncol(dat2.train)])
y <- as.matrix(dat2.train[, ncol(dat2.train)])

#tuning <- function(X,y, m_0= 50, m_max = 1000){
  # get all the possible index
  m_0 = 50
  m_max = 1000
  n <- nrow(X)
  d <- ncol(X)
  I_max <- sample(n, m_max)
  idx <- unique(c(seq(m_0, m_max, by= m_0),m_max))
  
  # test different lambdas
  lambdas <- 10^(seq(-6, 2, length.out = 10)) 

  ## Choose m that gives small drop in MSE
  ## TODO: need to be careful about the way we choose the drop though
  MSE <- c()
  # initialize the first 2 values of the MSE vector
  MSE[1] <- findMSE(1, X, y, lambdas, I_max,idx)
  MSE[2] <- findMSE(2, X, y, lambdas, I_max,idx)
  
  i <- 2
  while ((MSE[i-1] - MSE[i])/ MSE[i-1] > 0.001 & i <length(idx)) {
    i <- i+1
    MSE[i] <- findMSE(i, X, y, lambdas, I_max,idx)
  }
  
 
  
#}
  
# the best result is m= m_0 = 50
I <-  I_max[1:idx[4]]
X_scaled <- scale(X)
y_scaled <- scale(y)
R <- new_gauss_kern(X_scaled, X_scaled[I, ], b=d)
D <- R[I, ]
folds <- 5
lables <- sample(folds, n, replace = T)
cv_MSE <- rep(0, length(lambdas))

## TODO: speed up by using sapply
for (i in 1:length(lambdas)){
  lambda <- lambdas[i]
  MSE_vali <- rep(0, folds)
  for (v in 1:folds){
    # training indices in cross validation
    cv_idx <- which(lables != v)
    
    # validation MSE 
    dh_vali <- train_krr(y_scaled[cv_idx], R[cv_idx, ], D, lambda)$dh
    yh_vali <- (R[-cv_idx, ] %*% dh_vali + attributes(y_scaled)$`scaled:center`) * attributes(y_scaled)$`scaled:scale`
    MSE_vali[v] <- sqrt(mean((yh_vali- y[-cv_idx])^2))
  }
  
  # average validation MSE 
  cv_MSE[i] <- mean(MSE_vali)
}

# select optimal lambda
lambda_opt <- lambdas[which.min(cv_MSE)]
dh_train <- train_krr(y_scaled, R, D, lambda_opt)$dh
yh_train <- (R %*% dh_train + attributes(y_scaled)$`scaled:center`) * attributes(y_scaled)$`scaled:scale`
#train MSE
mean((yh_train- y)^2)

## TEST DATA
# X_test <- matrix(runif(100*p), nrow = 100, ncol = p)
# 
# X_test_scaled <- scale(X_test, 
#                 center = attributes(X_scaled)$`scaled:center`,
#                 scale = attributes(X_scaled)$`scaled:scale`)
# y_pure <- ftru(X_test)
# y_test <- y_pure + rnorm(100, sd = sgm)

# real data
# X_test <- as.matrix(winequality.red[-(1:1300), -12])

# y_test <- as.matrix(winequality.red[-(1:1300), 12])

X_test <- as.matrix(dat2.test[, -ncol(dat2.test)])
y_test <- as.matrix(dat2.test[, ncol(dat2.test)])
X_test_scaled <- scale(X_test,
                        center = attributes(X_scaled)$`scaled:center`,
                        scale = attributes(X_scaled)$`scaled:scale`)

yh_test <- ( new_gauss_kern(X_test_scaled, X_scaled[I, ], b=d)%*% dh_train +
               attributes(y_scaled)$`scaled:center`) * attributes(y_scaled)$`scaled:scale`
mean((yh_test- y_test)^2)
