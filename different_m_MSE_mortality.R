# For a fixed largish n, see the MSE with different m in real data

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
## assign labels to the training set
n_train <- nrow(dat2) - n_test
train_groups <- sample(1:5, size =n_train, replace = T)
idx_train_all <- (1:nrow(dat2))[-idx_test]

result <- data.frame(type = factor("nys", levels = c("nys","col")))
ms <- c(seq(10,200,10),seq(250,1000,50))
t <- 1
for (i in 1:5){
  ## training data
  idx_train <- idx_train_all[train_groups == i]
  ## training data
  X_init <- dat2[idx_train, -ncol(dat2)]
  y_init <- dat2[idx_train, ncol(dat2)]
  X <- scale(X_init) # do scaling at start and end only
  y <- scale(y_init)
  # scale the testing data according to the training
  X_test_scaled <- scale(X_test, 
                         center = attributes(X)$`scaled:center`,
                         scale = attributes(X)$`scaled:scale`)
  K_test <- new_gauss_kern(X_test_scaled, X, b=d)
  
  ## nystrom
  I_max <- sample(length(idx_train), max(ms))
  # construct the largest kernal 
  K_max <- kern_gauss(X[I_max, ], b = d)
  K_cv <- new_gauss_kern(X[-I_max, ], X[I_max, ], b = d)
  y_cv <- y[-I_max]
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
