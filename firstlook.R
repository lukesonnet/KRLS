library(dplyr)
library(tidyr)

# dat <- readRDS("~/Dropbox/ScalingKRLS/application/earlymortality/idhs_final_transformed_full.rds")
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

# HWWTHTPCT is (Weight for height percentile (CDC standards)). Consider hwhtapct as well.
# The main variables we've thought about as predictors and transformed are:
# kidbord, pregtermin, toilettype_new, drinkwtr_new, floor_new,cookfuel_new,wealthq urban,kidsex,
# bednetnum,agefrstmar,hheadsex,hheadagehh,edyrtotal,religion,maternal_age_month,
# kidbirthmo,prev_death,prev_death_full

# Choose outcome var:
yvar = "hwwthtpct"
#yvar = "hwhtapct"

# Predictor vars -- this is a work in progress; we can add more.
xvars = c("country", "urban","kidsexgirl","edyrtotal",
          "maternal_age_month", "prev_death_full","wealth_q1","wealth_q2","wealth_q3","wealth_q4",
          "watergood","pregtermin","hhousehead_fem","toilet_flush","toilet_pit")
# xvars = c("urban","kidsexgirl","edyrtotal",
#           "maternal_age_month", "prev_death_full","wealth_q1","wealth_q2","wealth_q3","wealth_q4",
#           "watergood","pregtermin","hhousehead_fem","toilet_flush","toilet_pit") 
#omitting wealth q5 and toilet_unimproved for intercepts. 

dat2 = na.omit(dat[, c(xvars, yvar)])
# drop Bangladesh and India
# dat2 <- subset(dat2, subset = country != "Bangladesh"  & country != "India")
# # create dummy var for country
# dummies <- dummies::dummy(dat2$country, sep = "_") 
# dat2 <- cbind(dummies, dat2[, -1]) # combine the dummy variable to the data 
# names(dat2) <- make.names(names(dat2)) # make legitimate names
# xvars <- names(dat2)[-ncol(dat2)] # update the name of X variable

ntrain = 5*1e3
ntest = 15*1e3
set.seed(12345)
train.ind = sample(1:nrow(dat2), ntrain)
dat2.train = dat2[train.ind,]
dat2.test = dat2[-train.ind,]

test.ind = sample(1:nrow(dat2.test), ntest)
dat2.test.small = dat2.test[test.ind,]

form = formula(paste(yvar, paste(xvars, collapse = " + "), 
              sep = " ~ "))
               
lm.small.out = lm(form, data = dat2.train)
summary(lm.small.out)
summary(lm.small.out)$r.squared

lm.predict = predict.lm(lm.small.out, newdata = dat2.test.small)
cor(lm.predict, dat2.test.small[,1])^2

### KRLS
krls.small.out = KRLS2::krls(X = dat2.train[,1:(ncol(dat2)-1)], y=dat2.train[,ncol(dat2)], epsilon = .001)
krls.small.out$lastkeeper
yhat = krls.small.out$fitted
cor(yhat, dat2.train[,1])^2
mean((yhat-dat2.train[,1])^2)

krls.predict = KRLS2:::predict.krls2(object=krls.small.out, newdata = dat2.test.small[,-ncol(dat2)])

cor(krls.predict$fit, dat2.test.small[,1])^2
mean((krls.predict$fit - dat2.test.small[,1])^2)
