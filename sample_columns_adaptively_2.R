# search for m adaptively version 2
# out sample prediction using Rd directly
rm(list = ls())
library(dplyr)
library(ggplot2)
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

ns <- round(10^seq(3.1, 4.3, 0.3))
#ns <- c(200, 500, 1000, 1500, 2000)
#ns <- 1000
t <- 1
Nrep <- 5 #number of simulated training data set
lambdas = 10^(seq(-6, 2, length.out = 10))
result <- data.frame(type=factor("lm", levels=c("lm","rf","nys", "truc", "col", "ske","KRLS")), n=1, mse=0, dt = 0, rep=1)
m_selected <- matrix(nrow = length(ns), ncol = Nrep) # save the selected m for each n
row.names(m_selected) <- paste0("n=", ns)
colnames(m_selected) <- 1:Nrep

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
    
    
    m_0 <- 50 # number of observations for each batch
    n_batch <- 5 # number of batches
    index_all <- sample(n, m_0*n_batch) # get all the index of batches
    alpha <- 0.25 # percentile for worst MSE
    
    result_batches <- vector(mode="list", length=2)
    names(result_batches) <- c("MSE", "I")
    result_batches[[1]] <- c()
    result_batches[[2]] <- list()
    
    dt <- system.time({for(i in 1:n_batch){
      # index for batch i
      sub_index_i <- ((i-1)*m_0+1) : (i*m_0)
      index_i <- index_all[sub_index_i]
      # KRLS on batch
      KRLS_result <- krls(X = X[index_i, ], y=y[index_i]) # do full KRLS
      krls.predict <- predict.krls2(object=KRLS_result, newdata = X[-index_i, ]) # prediction on the rest of data
      yh <- krls.predict$fit # predicted value
      rs <- (yh - y[-index_i])^2 # residual squared
      # select observations that give (1-alpha) percentile MSE
      threshold <- max(rs)*alpha
      index_left <- (1:n)[-index_i]
      I <- c(index_i, index_left[rs > threshold])
      # result <- list(MSE = mean(rs),
      #                I = I)
      result_batches[[1]][i] <- mean(rs)
      result_batches[[2]][[i]] <- I
      #return(result)
    }
    
    # select the batch that gives the smallest MSE
    I <- result_batches[[2]][[which.min(result_batches[[1]])]]
     
    # # train the model based on I
    R <- new_gauss_kern(X, X[I, ], b=d)
    D <- R[I, ]
    lambda_opt <- select_lambda(R, D, y)
    dh_train <- train_krr(y, R, D, lambda_opt)$dh
    })["elapsed"] 
    
    ## get the test error
    # get back to the original scale for the predicted y
    y_predicted <- K_test[ , I] %*% dh_train *
      attributes(y)$`scaled:scale` + attributes(y)$`scaled:center`
    # find the testing MSE
    result[t,"mse"] <- mean((y_predicted-y_test)^2)
    result[t,"rep"] <- rep
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"type"] <- "nys"
    m_selected[k, rep] <-  length(I)# save the selected m
    t <- t + 1
  }
}

## MSE plot
result2 <- result %>% group_by(type,n) %>% summarize(mean_mse = mean(mse), mean_dt = mean(dt))
#ggplot(subset(result2, type != "lm" & type != "rf"), aes(x=n,y=mean_mse,color=type)) + 
ggplot(result2, aes(x=n,y=mean_mse,color=type)) +    
  geom_line(size=1.1) +
  xlab("n(log)") + ylab("MSE(log)") #+ scale_y_log10() + scale_x_log10()+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.8),
        legend.text = element_text(size=14),
        text = element_text(size=20)) 
