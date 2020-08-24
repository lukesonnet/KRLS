# simple version of cross validation
rm(list = ls())
library(randomForest)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(123)
source("select_m.R")

## load full data and clean it
dat <- readRDS('C:/Users/Linfan/Downloads/idhs_final_transformed_full.rds')
# consider include country as variable
# drop Bangladesh and India

# Recodings
dat <- dat %>% mutate(infmort = as.numeric(mortality.under12m)-1,
                      urban = 2-(as.numeric(urban)),
                      kidsexgirl = as.numeric(kidsex)-1,
                      prev_death_full = as.numeric(prev_death_full)-1,
                      watergood = as.numeric(drinkwtr_new)-1,
                      pregtermin = as.numeric(pregtermin)-1,
                      hhousehead_fem = as.numeric(hheadsex)-1,
                      hwwthtpct = hwwthtpct/1e4,
                      hwhtapct = hwhtapct/1e4
)

# One-hot coding of wealth quintiles
dat <- dat %>% mutate(wealth_q1 = ifelse(wealthq2=='firstq', 1, 0),
                      wealth_q2 = ifelse(wealthq2=='secondq', 1, 0),
                      wealth_q3 = ifelse(wealthq2=='thirdq', 1, 0),
                      wealth_q4 = ifelse(wealthq2=='fourthq', 1, 0),
                      wealth_q5 = ifelse(wealthq2=='fifthq', 1, 0),
                      toilet_flush = ifelse(toilettype_new=="flush",1,0),
                      toilet_pit = ifelse(toilettype_new=="pit",1,0),
                      toilet_unimproved = ifelse(toilettype_new=="unimproved",1,0)
)
yvar = "hwwthtpct"
xvars = c("country", "urban","kidsexgirl","edyrtotal",
          "maternal_age_month", "prev_death_full","wealth_q1","wealth_q2","wealth_q3","wealth_q4",
          "watergood","pregtermin","hhousehead_fem","toilet_flush","toilet_pit")
dat2 = na.omit(dat[, c(xvars, yvar)])
# drop Bangladesh and India
dat2 <- subset(dat2, subset = country != "Bangladesh"  & country != "India")
# create dummy var for country
dummies <- dummies::dummy(dat2$country, sep = "_")
dat2 <- cbind(dummies, dat2[, -1]) # combine the dummy variable to the data
names(dat2) <- make.names(names(dat2)) # make legitimate names
xvars <- names(dat2)[-ncol(dat2)] # update the name of X variable

## Test data
n_test = 5*1e3
idx_test <- sample(nrow(dat2), n_test)
X_test <- dat2[idx_test, -ncol(dat2)]
y_test <- dat2[idx_test, ncol(dat2)]
d <- ncol(X_test)

## set up the saving vector
ns <-round(10^seq(3.1, log10(120000), length.out = 6))
#ns <- c(200, 500, 1000, 1500, 2000)
idx_train_max <- sample((1:nrow(dat2))[-idx_test], max(ns))
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
    idx_train <- idx_train_max[1:n]
    ## training data
    X_init <- dat2[idx_train, -ncol(dat2)]
    y_init <- dat2[idx_train, ncol(dat2)]
    X <- scale(X_init) # do scaling at start and end only
    y <- scale(y_init)
    # scale the testing data according to the training
    X_test_scaled <- scale(X_test, 
                           center = attributes(X)$`scaled:center`,
                           scale = attributes(X)$`scaled:scale`)
    
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
    y_predicted <- new_gauss_kern(X_test_scaled, X[I, ], b=d)%*% dh_train * 
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
    y_predicted <- new_gauss_kern(X_test_scaled, X[I, ], b=d)%*% dh_train * 
      attributes(y)$`scaled:scale` + attributes(y)$`scaled:center`
    # find the testing MSE
    result[t,"mse"] <- mean((y_predicted-y_test)^2)
    result[t,"dt"] <- dt
    result[t,"n"] <- n
    result[t,"rep"] <- rep
    t <- t + 1
    
    
    # ## sketching
    # dt <- system.time({
    #   K <- kern_gauss(X, d)
    #   S <- matrix(rnorm(m*n), nrow = n, ncol = m)/sqrt(m)
    #   R <- K %*% S
    #   D <- t(S) %*% K %*% S
    #   lambda_opt <- select_lambda(R, D, y)
    #   dh_train <- train_krr(y, R, D, lambda_opt)$dh
    # })["elapsed"] + time_m
    # result[t,"type"] <- "ske"
    # # get back to the original scale for the predicted y
    # y_predicted <- new_gauss_kern(X_test_scaled, X[I, ], b=d) %*% S %*% dh_train *
    # attributes(y)$`scaled:scale`+ attributes(y)$`scaled:center`
    # # find the testing MSE
    # result[t,"mse"] <- mean((y_predicted-y_test)^2)
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"rep"] <- rep
    # t <- t + 1
    # 
    # ## truncation
    # dt <- system.time( m_krls <- krls(X = X_init, y = y_init, b = d, epsilon = 0.001))["elapsed"]
    # y_predicted = predict.krls2(object=m_krls, newdata = X_test)$fit
    # result[t,"mse"] <- mean((y_predicted-y_test)^2)
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"type"] <- "truc"
    # result[t,"rep"] <- rep
    # t <- t + 1
    # 
    # ## original KRLS
    # dt <- system.time( m_krls <- krls(X = X_init, y=y_init, b = d))["elapsed"]
    # y_predicted = predict.krls2(object=m_krls, newdata = X_test)$fit
    # result[t,"mse"] <- mean((y_predicted-y_test)^2)
    # result[t,"dt"] <- dt
    # result[t,"n"] <- n
    # result[t,"type"] <- "KRLS"
    # result[t,"rep"] <- rep
    # t <- t + 1
    
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
  xlab("n") + ylab("Time (log scale)") + scale_y_log10() + scale_x_log10()+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.4),
        legend.text = element_text(size=14),
        text = element_text(size=20)) 

