# Truncation, linear regression, Nystrom together

# library(RcppArmadillo)
# # source("regbench.R")
# Rcpp::sourceCpp('cpp/krls_solve.cpp')
# Rcpp::sourceCpp('cpp/krls_kernels.cpp')
library(randomForest)
library(dplyr)
source("kernel_mod.R")
set.seed(1)

lineplot <- function(x,y) {
  idx <- order(x)
  plot(x[idx],y[idx],type="l")
}

ns <- round(seq(10,2000,length.out = 10))
#ns <- c(500, 1000, 1500, 2000)
num_n <- length(ns)



# Need to fix the function 
Ntru <- 50
p <- 5
d <- p
Xtru <- scale( matrix(runif(Ntru*p), nrow = Ntru, ncol = p) )
#Xtru <- matrix(runif(Ntru*p), nrow = Ntru, ncol = p)
coeff <- rnorm(Ntru)
ftru <- function(X) {
  K_X_Xtru <- new_gauss_kern(Xtru, X, d)
  t(K_X_Xtru) %*% coeff
}
sgm <- sqrt(2)*sd(ftru(Xtru))


compute_mse <- function(y,yhat) mean((y-yhat)^2)

t <- 1
result <- data.frame(type=factor("lm", levels=c("lm","rf","nys_ada","nys_rnd","KRLS")), n=1, mse=0, dt = 0, rep=1)
Nrep <- 5
for (rep in 1:Nrep) {
  cat('.')
  for (k in 1:num_n) {
      ##### Naive simulation #####
      n <- ns[k] # total num of observations
      
      X <- scale( matrix(runif(n*p), nrow = n, ncol = p) )
      #X <- matrix(runif(n*p), nrow = n, ncol = p)
      # X <- scale( matrix( sort(runif(n*p)), nrow = n, ncol = p) )
      y_pure  <- ftru(X)
      y <- y_pure + rnorm(n, sd = sqrt(2))
      #X <- scale(X)
      
      # did not select optimal lambda?
      lambda <- 1
      dt <- system.time( { 
        I <- krls_select_m(y, X, lambda, m = 25, ncv=5,  press_pct = 0.5);
        krls_model <- krls2(y,X,lambda,I) } )["elapsed"]
      m_ada <- length(I)
      result[t,"type"] <- "nys_ada"
      result[t,"mse"] <- compute_mse(predict(krls_model), y_pure)
      result[t,"dt"] <- dt
      result[t,"n"] <- n
      result[t,"m"] <- m_ada
      result[t,"rep"] <- rep
      t <- t+1
      
      dt <- system.time( { 
        J <- sample(n, m_ada);
        krls_model2 <- krls2(y,X,lambda,J) } )["elapsed"]
      result[t,"type"] <- "nys_rnd"
      result[t,"mse"] <- compute_mse(predict(krls_model2), y_pure)
      result[t,"dt"] <- dt
      result[t,"n"] <- n
      result[t,"m"] <- m_ada
      result[t,"rep"] <- rep
      t <- t+1
      
      data_df <- data.frame(X, y)
      # #linear regression
      
      dt <- system.time( lm_model <- lm(y~X, data = data_df) )["elapsed"]
      result[t,"type"] <- "lm"
      result[t,"mse"] <- compute_mse(predict.lm(lm_model), y_pure)
      result[t,"dt"] <- dt
      result[t,"n"] <- n
      result[t,"rep"] <- rep
      # lm_model <- lm(X~y, lm_df, subset = train_idx)
      # test_MSE_lm[k] <- mean((predict.lm(lm_model, lm_df[-train_idx, ])-y_test)^2)
      t <- t+1
      
      # random forest
      dt <- system.time( m_rf <- randomForest(y~. , data = data_df)  )["elapsed"]
      result[t,"type"] <- "rf"
      result[t,"mse"] <- compute_mse(predict(m_rf), y_pure) 
      result[t,"dt"] <- dt
      result[t,"n"] <- n
      result[t,"rep"] <- rep
      t <- t+1
      
      # m_rf <- randomForest(y~. , data = data_df, subset = train_idx)
      # m_rf <- randomForest(y~. , data = lm_df, subset = train_idx)
      # pred_rf <- predict(m_rf, lm_df[-train_idx, ])
      # test_MSE_rf[k] <- mean((pred_rf-y_test)^2)
      ## original KRLS
      dt <- system.time( m_krls <- krls(X = X, y=y, b = d))["elapsed"]
      ## out-sample test
      #y_predicted = predict.krls2(object=m_krls, newdata = X_test)$fit
      ## pure value test
      y_predicted = m_krls$fitted
      # result[t,"mse"] <- mean((y_predicted-y_test)^2)
      result[t,"mse"] <- mean((y_predicted-y_pure)^2)
      result[t,"dt"] <- dt
      result[t,"n"] <- n
      result[t,"type"] <- "KRLS"
      result[t,"rep"] <- rep
      t <- t + 1
  }
}
# df_MSE <- reshape::melt(data.frame(n = ns,
#                                    nystrom = test_MSE_nys, 
#                                    truncation = test_MSE_tru,
#                                    lm = MSE_lm,
#                                    rf = MSE_rf),
#                         id = "n")

library(ggplot2)
# ggplot(result, aes(x=as.factor(n),y=mse,fill=type)) + geom_boxplot() +
#    theme_bw() + xlab("n") + ylab("MSE")
# 
# ggplot(result, aes(x=as.factor(n),y=mse,color=type)) + geom_line() +
#    theme_bw() + xlab("n") + ylab("MSE")
  

result2 <- result %>% group_by(type,n) %>% summarize(mean_mse = mean(mse), mean_dt = mean(dt), mean_m = mean(m))
ggplot(result2, aes(x=n,y=mean_mse,color=type)) + geom_line(size=1.1) +
  xlab("n (log scale)") + ylab("MSE (log scale)") + scale_y_log10() +  scale_x_log10() +
  #scale_x_continuous(limits=c(ns[1], ns[num_n]), expand = c(0.02, 0)) +
  theme_bw() +
  theme(#axis.text.x = element_blank(), 
    #axis.ticks.x = element_blank(),
    legend.background=element_blank(), 
    legend.title=element_blank(), 
    legend.position = c(0.85, 0.3),
    legend.text = element_text(size=14),
    text = element_text(size=16)) #+
  # guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch")) 
#ggsave('mse.png',width = 6,height = 5)

ggplot(result2, aes(x=n,y=mean_dt,color=type)) + geom_line(size=1.1) +
  xlab("n (log scale)") + ylab("Time (log scale)") + scale_y_log10() + scale_x_log10() +
  #scale_x_continuous(limits=c(ns[1], ns[num_n]), expand = c(0.02, 0)) +
  theme_bw() +
  theme(#axis.text.x = element_blank(), 
    #axis.ticks.x = element_blank(),
    legend.background=element_blank(), 
    legend.title=element_blank(), 
    legend.position = c(0.85, 0.3),
    legend.text = element_text(size=14),
    text = element_text(size=16))  
#ggsave('time.png',width = 6,height = 5)

fit <- lm(log10(mean_mse) ~ log10(n), data= result2 %>% filter(type == "rf"))
summary(fit)

ggplot(filter(result2, type == "nys_ada"), aes(x=n,y=mean_m)) + geom_line(size=1.1) +
  xlab("n (log scale)") + ylab("m (log scale)") + scale_y_log10() + scale_x_log10() +
  #scale_x_continuous(limits=c(ns[1], ns[num_n]), expand = c(0.02, 0)) +
  theme_bw() +
  theme(#axis.text.x = element_blank(), 
    #axis.ticks.x = element_blank(),
    legend.background=element_blank(), 
    legend.title=element_blank(), 
    legend.position = c(0.85, 0.3),
    legend.text = element_text(size=14),
    text = element_text(size=16))  
#ggsave('m.png',width = 6,height = 5)

# ggplot(df_MSE, aes(x=n, y=value, colour=variable)) +
#     geom_line(size = 1.1) +
#     ylab("Test MSE") +
#     ggtitle("Test MSE using Nystrom, truncation, lm, random forest (Simple Simulation)")+
#     theme(legend.position="bottom", legend.title=element_blank(),
#           plot.title = element_text(hjust = 0.5))
