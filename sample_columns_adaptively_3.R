# search for m adaptively version 3
rm(list = ls())
library(dplyr)
library(randomForest)
library(ggplot2)
source("select_m.R")
set.seed(1)

## set the original generating function
Ntru <- 50
d <- 5 # number of covariates
# set all the bandwidth to be d
#Xtru <- matrix(runif(Ntru*d), nrow = Ntru, ncol = d)
Xtru <- scale( matrix(runif(Ntru*d), nrow = Ntru, ncol = d) )
coeff <- rnorm(Ntru)
ftru <- function(X) {
  K_X_Xtru <- new_gauss_kern(Xtru, X, d)
  t(K_X_Xtru) %*% coeff
}
sgm <- sqrt(2)*sd(ftru(Xtru))

## TEST DATA
# n_test <-  2000
# X_test <- matrix(runif(n_test*d), nrow = n_test, ncol = d)
# y_pure <- ftru(X_test)
# y_test <- y_pure + rnorm(n_test, sd = sgm)

#ns <- round(10^seq(3.1, 4.3, 0.3))
ns <- round(seq(150,2000,length.out = 10))
#ns <- c(500, 1000, 1500, 2000)
#ns <- 1000
t <- 1

Nrep <- 5 #number of simulated training data set
lambdas = 10^(seq(-6, 2, length.out = 10)) 
m_selected <- matrix(nrow = length(ns), ncol = Nrep) # save the selected m for each n
row.names(m_selected) <- paste0("n=", ns)
colnames(m_selected) <- 1:Nrep
result <- data.frame() # save all the results in this data frame

for (rep in 1:Nrep) {
  cat('.')
  for (k in 1:length(ns)) {
    n <- ns[k]
    
    ## simulate training data with different sizes for Nrep times
    #X_init <- matrix(runif(n*d), nrow = n, ncol = d)
    X <- scale( matrix(runif(n*d), nrow = n, ncol = d) )
    y_pure  <- ftru(X)
    #y_pure  <- ftru(X_init) # response value without noise
    y <- y_pure + rnorm(n, sd = sgm)
    #X <- scale(X_init) # do scaling at start and end only
    # K_train used for y pure testing only
    K_train <- kern_gauss(X, b=d)
    
    ## scale the testing data according to the training
    # X_test_scaled <- scale(X_test, 
    #                        center = attributes(X)$`scaled:center`,
    #                        scale = attributes(X)$`scaled:scale`)
    #K_test <- new_gauss_kern(X_test_scaled, X, b=d)
    
    # nystrom with adaptive selection
    m_0 <- 25 # number of observations for each batch
    n_batch <- 5 # number of batches
    index_all <- sample(n, m_0*n_batch) # get all the index of batches
    alpha <- 0.25 #  thresholding MSE 
    
    result_batches <- vector(mode="list", length=2)
    names(result_batches) <- c("MSE", "I")
    result_batches[[1]] <- c()
    result_batches[[2]] <- list()
    
    dt <- system.time({for(i in 1:n_batch){
      # # index for batch i
      # sub_index_i <- ((i-1)*m_0+1) : (i*m_0)
      # index_i <- index_all[sub_index_i]
      index_i <- sample(n, m_0)
      # KRLS on batch with the same lambda = 1
      KRLS_result <- krls(X = X[index_i, ], y=y[index_i], b=d) 
      # prediction on the rest of data
      krls.predict <- predict.krls2(object=KRLS_result, newdata = X[-index_i, ]) 
      yh <- krls.predict$fit # predicted value
      rs <- (yh - y[-index_i])^2 # residual squared
      # select observations that give (1-alpha) percentile MSE
      threshold <- max(rs)*alpha
      index_left <- (1:n)[-index_i]
      I <- c(index_i, index_left[rs > threshold])
      result_batches[[1]][i] <- mean(rs)
      result_batches[[2]][[i]] <- I
      #return(result)
    }
      
      # select the batch that gives the smallest MSE
      I <- result_batches[[2]][[which.min(result_batches[[1]])]]
      
      ## train the model based on I and output dh
      R <- new_gauss_kern(X, X[I, ], b=d)
      D <- R[I, ]
      # cv to select lambda
      lambda_opt_1 <- select_lambda(R, D, y)
      #lambda_opt_1 <- lambda_opt_2 <- 1
      dh <- train_krr(y, R, D, lambda_opt_1)$dh
      
      # output on ch
      lambda_opt_2 <- select_lambda_2(R, D, y, X) # cv on ch 
      dh <- train_krr(y, R, D, lambda_opt_2)$dh
      ch <- R %*% solve(crossprod(R) + diag(length(I)), D %*% dh)
    })["elapsed"] 
    
    ## get the test error
    ## out-sample testing
    # y_predicted_1 <- K_test[ , I] %*% dh 
    # y_predicted_2 <- K_test %*% ch
    ## pure value testing
    y_predicted_1 <- R %*% dh
    y_predicted_2 <- K_train %*% ch
    # find the testing MSE
    #result[t,"mse"] <- mean((y_predicted_1-y_test)^2)
    result[t,"mse"] <- mean((y_predicted_1-y_pure)^2)
    result[t,"rep"] <- rep
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"type"] <- "nys_ada_dh"
    t <- t + 1
    
    #result[t,"mse"] <- mean((y_predicted_2-y_test)^2)
    result[t,"mse"] <- mean((y_predicted_2-y_pure)^2)
    result[t,"rep"] <- rep
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"type"] <- "nys_ada_ch"
    
    m_selected[k, rep] <-  length(I)# save the selected m
    t <- t + 1
    
    # # nystrom with random selection 
    # 
    # dt <- system.time({I_ran <- sample(n, length(I))
    # # # train the model based on I
    # R <- new_gauss_kern(X, X[I_ran, ], b=d)
    # D <- R[I_ran, ]
    # lambda_opt <- select_lambda(R, D, y)
    # #lambda_opt <- 1
    # dh <- train_krr(y, R, D, 1)$dh
    # ch <- R %*% solve(crossprod(R) + diag(length(I)), D %*% dh)
    # })["elapsed"] 
    # 
    # ## get the test error
    # ## out-sample test
    # # y_predicted_1 <- K_test[ , I_ran] %*% dh 
    # # y_predicted_2 <- K_test %*% ch
    # ## pure value testing
    # y_predicted_1 <- R %*% dh
    # y_predicted_2 <- K_train %*% ch
    # # find the testing MSE
    # # result[t,"mse"] <- mean((y_predicted_1-y_test)^2)
    # # result[t,"mse"] <- mean((y_predicted_1-y_pure)^2)
    # # result[t,"rep"] <- rep
    # # result[t,"dt"] <- dt
    # # result[t,"n"] <- n
    # # result[t,"type"] <- "nys_ran_dh"
    # # t <- t + 1
    # 
    # #result[t,"mse"] <- mean((y_predicted_2-y_test)^2)
    # result[t,"mse"] <- mean((y_predicted_2-y_pure)^2)
    # result[t,"rep"] <- rep
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"type"] <- "nys_ran_ch"
    # t <- t + 1
    
    # ## linear regression
    # # put X and y in a data frame
    # data_df <- data.frame(X,y)
    # # put X and y in a formula
    # data_formula <- formula(paste(tail(names(data_df), 1), 
    #                               paste(names(data_df)[-ncol(data_df)], collapse = " + "), 
    #                               sep = " ~ "))
    # # put the scaled version of X and y in the LM and record the time
    # dt <- system.time( m_lm <- lm(data_formula, data = data_df) )["elapsed"]
    # # dt <- system.time( m_lm <- lm(X~y, data = data_df) )["elapsed"]
    # result[t,"type"] <- "lm"
    # ## out-sample test
    # #y_predicted <- predict.lm(m_lm, newdata = data.frame(X_test_scaled))
    # ## pure value test
    # y_predicted <- predict.lm(m_lm)
    # # find the testing MSE
    # # result[t,"mse"] <- mean((y_predicted-y_test)^2)
    # result[t,"mse"] <- mean((y_predicted-y_pure)^2)
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"rep"] <- rep
    # t <- t + 1
    
    # ## original KRLS
    # dt <- system.time( m_krls <- krls(X = X, y=y, b = d))["elapsed"]
    # ## out-sample test
    # #y_predicted = predict.krls2(object=m_krls, newdata = X_test)$fit
    # ## pure value test
    # y_predicted = m_krls$fitted
    # # result[t,"mse"] <- mean((y_predicted-y_test)^2)
    # result[t,"mse"] <- mean((y_predicted-y_pure)^2)
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"type"] <- "KRLS"
    # result[t,"rep"] <- rep
    # t <- t + 1
    # 
    # # random forest
    # dt <- system.time( m_rf <- randomForest(data_formula, data = data_df)  )["elapsed"]
    # result[t,"type"] <- "rf"
    # ## out-sample test
    # # y_predicted <- predict(m_rf, newdata = data.frame(X_test_scaled))
    # ## pure value test
    # y_predicted <- predict(m_rf)
    # #result[t,"mse"] <- mean((y_predicted-y_test)^2)
    # result[t,"mse"] <- mean((y_predicted-y_pure)^2)
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"rep"] <- rep
    # t <- t + 1
  }
}

## MSE plot
result2 <- result %>% group_by(type,n) %>% summarize(mean_mse = mean(mse), mean_dt = mean(dt))
#ggplot(subset(result2, type != "lm" & type != "rf"), aes(x=n,y=mean_mse,color=type)) + 
ggplot(result2, aes(x=n,y=mean_mse,color=type)) +    
  geom_line(size=1.1) +
  xlab("n (log scale)") + ylab("MSE (log scale)") + scale_y_log10() + scale_x_log10()+
theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.8),
        legend.text = element_text(size=14),
        text = element_text(size=20)) 
