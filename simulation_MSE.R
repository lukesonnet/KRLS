# simple version of cross validation
rm(list = ls())
library(randomForest)
library(ggplot2)
library(dplyr)
require(RSpectra)
set.seed(1)
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

## set up saving vectors
# ns <- round(10^seq(3.1, 5.3, 0.3))
ns <- c(200, 500, 1000, 1500, 2000)
t <- 1
Nrep <- 5 #number of simulated training data set
lambdas = 10^(seq(-6, 2, length.out = 10))
result <- data.frame(type=factor("lm", levels=c("lm","rf","nys", "truc", "col", "ske","KRLS")), n=1, mse=0, dt = 0, rep=1)
m_selected <- matrix(nrow = length(ns), ncol = Nrep) # save the selected m for each n
row.names(m_selected) <- paste0("n=", ns)
colnames(m_selected) <- 1:Nrep

## find test_MSE in each methods
for (rep in 1:Nrep) {
  cat('.')
  for (k in 1:length(ns)) {
    n <- ns[k]
    
    ## training data
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
    
    
    ## linear regression
    # put X and y in a data frame
    data_df <- data.frame(X,y)
    # put X and y in a formula
    data_formula <- formula(paste(tail(names(data_df), 1), 
                     paste(names(data_df)[-ncol(data_df)], collapse = " + "), 
                  sep = " ~ "))
    # put the scaled version of X and y in the LM and record the time
    dt <- system.time( m_lm <- lm(data_formula, data = data_df) )["elapsed"]    
    result[t,"type"] <- "lm"
    # get back to the original scale for the predicted y
    y_predicted <- predict.lm(m_lm, newdata = data.frame(X_test_scaled))*
      attributes(y)$`scaled:scale` + attributes(y)$`scaled:center`
    # find the testing MSE
    result[t,"mse"] <- mean((y_predicted-y_test)^2)
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"rep"] <- rep
    t <- t + 1
    
    ## random forest
    dt <- system.time( m_rf <- randomForest(data_formula, data = data_df)  )["elapsed"]
    result[t,"type"] <- "rf"
    # get back to the original scale for the predicted y
    y_predicted <- predict(m_rf, newdata = data.frame(X_test_scaled))*
      attributes(y)$`scaled:scale` + attributes(y)$`scaled:center`
    result[t,"mse"] <- mean((y_predicted-y_test)^2)
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"rep"] <- rep
    t <- t + 1

    ## find optimal indices for nystrom and naive column sampling
    time_m <- system.time(  I <-  select_m(X, y) )["elapsed"]
    m <- length(I)
    m_selected[k, rep] <- m # save the selected m
    
    ## nystrom
    dt <- system.time( {
      R <- new_gauss_kern(X, X[I, ], b=d)
      D <- R[I, ]
      lambda_opt <- select_lambda(R, D, y)
      dh_train <- train_krr(y, R, D, lambda_opt)$dh
    })["elapsed"] + time_m
    result[t,"type"] <- "nys"
    # get back to the original scale for the predicted y
    y_predicted <- K_test[ , I] %*% dh_train * 
      attributes(y)$`scaled:scale`+ attributes(y)$`scaled:center`
    # find the testing MSE
    result[t,"mse"] <- mean((y_predicted-y_test)^2)
    result[t,"dt"] <- dt
    result[t,"n"] <- n

    result[t,"rep"] <- rep
    t <- t + 1

    ## naive column sampling
    dt <- system.time( {
      R <- new_gauss_kern(X, X[I, ], b=d)
      D <- diag(length(I))
      lambda_opt <- select_lambda(R, D, y)
      dh_train <- train_krr(y, R, D, lambda_opt)$dh
    } )["elapsed"] + time_m
    result[t,"type"] <- "col"
    # get back to the original scale for the predicted y
    y_predicted <- K_test[ , I] %*% dh_train * 
      attributes(y)$`scaled:scale` + attributes(y)$`scaled:center`
    # find the testing MSE
    result[t,"mse"] <- mean((y_predicted-y_test)^2)
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"rep"] <- rep
    t <- t + 1
# 
#     ## sketching
#     dt <- system.time({
#       K <- kern_gauss(X, d)
#       S <- matrix(rnorm(m*n), nrow = n, ncol = m)/sqrt(m)
#       R <- K %*% S
#       D <- t(S) %*% K %*% S
#       lambda_opt <- select_lambda(R, D, y)
#       dh_train <- train_krr(y, R, D, lambda_opt)$dh
#     })["elapsed"] + time_m
#     result[t,"type"] <- "ske"
#     # get back to the original scale for the predicted y
#     y_predicted <- K_test %*% S %*% dh_train *
#                        attributes(y)$`scaled:scale`+ attributes(y)$`scaled:center`
#     # find the testing MSE
#     result[t,"mse"] <- mean((y_predicted-y_test)^2)
#     result[t,"dt"] <- dt
#     result[t,"n"] <- n
#     result[t,"rep"] <- rep
#     t <- t + 1
# 
#     ## truncation
#     dt <- system.time({
#       K <- kern_gauss(X, d)
#       out <- eigs(K,m)
#       R <- out$vectors
#       D <- diag(1/out$values)
#       lambda_opt <- select_lambda(R, D, y)
#       dh_train <- train_krr(y, R, D, lambda_opt)$dh
#     })["elapsed"]
#     y_predicted <- K_test %*% R %*% D %*% dh_train *
#       attributes(y)$`scaled:scale`+ attributes(y)$`scaled:center`
#     result[t,"mse"] <- mean((y_predicted-y_test)^2)
#     result[t,"dt"] <- dt
#     result[t,"n"] <- n
#     result[t,"type"] <- "truc"
#     result[t,"rep"] <- rep
#     t <- t + 1
#     
#     ## original KRLS
#     dt <- system.time( m_krls <- krls(X = X_init, y=y_init, b = d))["elapsed"]
#     y_predicted = predict.krls2(object=m_krls, newdata = X_test)$fit
#     result[t,"mse"] <- mean((y_predicted-y_test)^2)
#     result[t,"dt"] <- dt
#     result[t,"n"] <- n
#     result[t,"type"] <- "KRLS"
#     result[t,"rep"] <- rep
#     t <- t + 1
    
  }
}

# plot

## MSE plot
result2 <- result %>% group_by(type,n) %>% summarize(mean_mse = mean(mse), mean_dt = mean(dt))
#ggplot(subset(result2, type != "lm" & type != "rf"), aes(x=n,y=mean_mse,color=type)) + 
ggplot(result2, aes(x=n,y=mean_mse,color=type)) +    
  geom_line(size=1.1) +
  xlab("n(log)") + ylab("MSE(log)") + scale_y_log10() + scale_x_log10()+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.8),
        legend.text = element_text(size=14),
        text = element_text(size=20)) 

## Time plot
ggplot(result2, aes(x=n,y=mean_dt,color=type)) + geom_line(size=1.1) +
  xlab("n") + ylab("Time (log scale)") + scale_y_log10() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.4),
        legend.text = element_text(size=14),
        text = element_text(size=20)) 

