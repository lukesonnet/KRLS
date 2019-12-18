# For a fixed largish n, see the MSE with different m

rm(list = ls())
library(randomForest)
library(ggplot2)
library(dplyr)
require(RSpectra)
set.seed(123)
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

n <- 100000 # n_train
result <- data.frame(type = factor("nys", levels = c("nys","col")))
ms <- c(seq(10,200,10),seq(250,1000,50))
t <- 1
for (i in 1:4){
  ## training data
  set.seed(i) # change seed to vary training data
  
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
  
  ## nystrom
  I_max <- sample(n, max(ms))
  
  for(m in ms){
    MSE <- c()
    # select the corresponding index and train the model
    I <- I_max[1:m]
    
    ## Nystrom
    R <- new_gauss_kern(X, X[I, ], b=d)
    D <- R[I, ]
    lambda_opt <- select_lambda(R, D, y)
    dh_train <- train_krr(y, R, D, lambda_opt)$dh
    y_predicted <- K_test[ , I] %*% dh_train * 
      attributes(y)$`scaled:scale`+ attributes(y)$`scaled:center`
    result[t, "type"] <- "nys"
    result[t, "mse"] <- mean((y_predicted-y_test)^2)
    result[t, "seed"] <- i
    result[t, "m"] <- m
    t <- t+1
    
    ## column sampling
    R <- new_gauss_kern(X, X[I, ], b=d)
    D <- diag(m)
    lambda_opt <- select_lambda(R, D, y)
    dh_train <- train_krr(y, R, D, lambda_opt)$dh
    y_predicted <- K_test[ , I] %*% dh_train * 
      attributes(y)$`scaled:scale`+ attributes(y)$`scaled:center`
    result[t, "type"] <- "col"
    result[t, "mse"] <- mean((y_predicted-y_test)^2)
    result[t, "seed"] <- i
    result[t, "m"] <- m
    t <- t+1
    
    #return(MSE)
  }
}

result$seed <- factor(result$seed)

result1 <- subset(result, type=="nys")
#ggplot(subset(result2, type != "lm" & type != "rf"), aes(x=n,y=mean_mse,color=type)) + 
ggplot(result1, aes(x=m,y=mse,color=seed)) +    
  geom_line(size=1.1) +
  xlab("m") + ylab("MSE") + #scale_y_log10() +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.8),
        legend.text = element_text(size=14),
        text = element_text(size=20)) 
