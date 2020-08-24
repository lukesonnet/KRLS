rm(list = ls())
library(randomForest)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(123)


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
## assign labels to the training set
n_train <- nrow(dat2) - n_test
train_groups <- sample(1:5, size =n_train, replace = T)
idx_train_all <- (1:nrow(dat2))[-idx_test]

result <- data.frame()
t <- 1
for (i in 1:5){
  ## training data
  idx_train <- idx_train_all[train_groups <= i]
  n <- length(idx_train)
  X_init <- dat2[idx_train, -ncol(dat2)]
  y <- dat2[idx_train, ncol(dat2)]
  X <- scale(X_init) 
  # scale the testing data according to the training
  X_test_scaled <- scale(X_test, 
                         center = attributes(X)$`scaled:center`,
                         scale = attributes(X)$`scaled:scale`)
  Ktest <- createKXY(X_test_scaled, X)
  
  # # KRLS
  # fit_krls <- krls2(y, X)
  # result[t,"mse"] <- mean((fit_krls$yh-y_pure)^2)
  # result[t,"rep"] <- rep
  # result[t,"n"] <- n
  # result[t,"type"] <- "KRLS"
  # t <- t + 1
  
  # Nystrom with adaptive sampling
  fit_nys_ada <- krls_nys(y, X)
  I <- fit_nys_ada$I
  m_selected[k, rep] <- length(I)# save the selected m
  yh <- Ktest %*% fit_nys_ada$ch
  result[t,"mse"] <- mean((yh-y_test)^2)
  result[t,"rep"] <- rep
  result[t,"n"] <- n
  result[t,"type"] <- "nys_ada"
  t <- t + 1
  
  # Nystrom with random sampling
  fit_nys_ran <- krls_nys(y, X, I = sample(n, length(I)))
  yh <- Ktest %*% fit_nys_ran$ch
  result[t,"mse"] <- mean((yh-y_test)^2)
  result[t,"rep"] <- rep
  result[t,"n"] <- n
  result[t,"type"] <- "nys_ran"
  t <- t + 1
  
  ## linear regression
  # put X and y in a data frame
  data_df <- data.frame(X,y)
  # put X and y in a formula
  data_formula <- formula(paste(tail(names(data_df), 1), 
                                paste(names(data_df)[-ncol(data_df)], collapse = " + "), 
                                sep = " ~ "))
  m_lm <- lm(data_formula, data = data_df)
  result[t,"type"] <- "lm"
  yh <- X_test_scaled %*%  m_lm$coefficients
  result[t,"mse"] <- mean((yh-y_pure)^2)
  result[t,"n"] <- n
  result[t,"rep"] <- rep
  t <- t + 1
  
  # random forest
  m_rf <- randomForest(data_formula, data = data_df)
  result[t,"type"] <- "rf"
  ## out-sample test
  yh <- predict(m_rf, newdata = data.frame(X_test_scaled))
  result[t,"mse"] <- mean((yh-y_test)^2)
  #result[t,"mse"] <- mean((y_predicted-y_pure)^2)
  result[t,"n"] <- n
  result[t,"rep"] <- rep
  t <- t + 1
}

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