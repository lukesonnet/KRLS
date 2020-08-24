# search for m adaptively
source("select_m.R")


## set the original generating function
Ntru <- 100
d <- 20 # number of covariates
Xtru <- matrix(runif(Ntru*d), nrow = Ntru, ncol = d)
coeff <- rnorm(Ntru)
ftru <- function(X) {
  K_X_Xtru <- new_gauss_kern(Xtru, X, d)
  t(K_X_Xtru) %*% coeff
}
sgm <- sqrt(2)*sd(ftru(Xtru))


## TEST DATA
n_test <-  2000
X_test <- matrix(runif(n_test*d), nrow = n_test, ncol = d)
y_pure <- ftru(X_test)
y_test <- y_pure + rnorm(n_test, sd = sgm)

## training data
n <- 10000
X_init <- matrix(runif(n*d), nrow = n, ncol = d)
y_pure  <- ftru(X_init)
y_init <- y_pure + rnorm(n, sd = sgm)
X <- scale(X_init) # do scaling at start and end only
y <- scale(y_init)
# scale the testing data according to the training
X_test_scaled <- scale(X_test, 
                       center = attributes(X)$`scaled:center`,
                       scale = attributes(X)$`scaled:scale`)
K_test <- new_gauss_kern(X_test_scaled, X, b=d)


# determine m based on computing time
m <- 350
m_0 <- 50 # number of observations for each batch
n_batch <- 5 # number of batches
index_all <- sample(n, m_0*n_batch) # get all the index of batches

# construct the cross validation kernel matrix 
K_cv <- new_gauss_kern(X[-index_all, ], X[index_all, ], b = d)
y_cv <- y[-index_all]
  
# find lambda using the first set of the index
lambdas = 10^(seq(-6, 2, length.out = 10))
index_1 <- index_all[1:m_0]
K_1 <- kern_gauss(X[index_1, ], b = d)
y_1 <- y[index_1]
K_cv <- new_gauss_kern(X[-index_1, ], X[index_1, ], b = d)
y_cv <- y[-index_1]
MSE_1 <- sapply(lambdas, function(x){
  ch <- solve(K_1 + x*diag(m_0), y_1)
  yh <- K_cv[, 1:m_0] %*% ch 
  mean((yh - y_cv)^2)
})
lambda <- lambdas[which.min(MSE_1)]
MSE_1<- min(MSE_1)

MSE_batches <- sapply(2:n_batch, function(i){
  # index for batch i
  sub_index_i <- ((i-1)*m_0+1) : (i*m_0)
  index_i <- index_all[sub_index_i]
  # KRLS on batch
  K_i <- kern_gauss(X[index_i, ], b = d)
  y_i <- y[index_i]
  K_cv <- new_gauss_kern(X[-index_1, ], X[index_1, ], b = d)
  y_cv <- y[-index_1]
  ch <- solve(K_i + lambda*diag(m_0), y_i)
  yh <- K_cv[ , sub_index_i] %*% ch 
  return(mean((yh - y_cv)^2))
})

# combine with the MSE of the first batch
MSE_batches <- c(MSE_1, MSE_batches)
# select the batch with least MSE
i_min <- which.min(MSE_batches)
# save the index for the selected batch
I_0 <- index_all[((i_min-1)*m_0+1) : (i_min*m_0)]
# perform the full KRLS (again)
K_I_0 <- kern_gauss(X[I_0, ], b = d)
ch <- solve(K_I_0 + lambda*diag(m_0), y[I_0])
yh <- K_cv[ , ((i_min-1)*m_0+1) : (i_min*m_0)] %*% ch
MSE_cv <- (yh - y_cv)^2
# sampling probability based on MSE
samp_prob <- MSE_cv/sum(MSE_cv)

# determine the rest of the index
I_1 <- sample((1:n)[-index_all], size = m-m_0, prob = samp_prob)
I <- c(I_0, I_1)

# train the model based on I
R <- new_gauss_kern(X, X[I, ], b=d)
D <- R[I, ]
lambda_opt <- select_lambda(R, D, y)
dh_train <- train_krr(y, R, D, lambda_opt)$dh

## get the test error
# get back to the original scale for the predicted y
y_predicted <- K_test[ , I] %*% dh_train * 
  attributes(y)$`scaled:scale` + attributes(y)$`scaled:center`
# find the testing MSE
result[t,"mse"] <- mean((y_predicted-y_test)^2)

