
setwd("C:/Dropbox/Research/krlogit/")

#----------
# Installing/Documenting the package
#----------

# You need the devtools package
#install.packages("devtools")

install <- TRUE

if (install) {
  ## This takes the roxygen code in the .cpp and .r files and build the .Rd,
  ## NAMESPACE, and DESCRIPTION Files (as well as important c++ files that
  ## install also does)
  devtools::document("KRLS2")
  ## Installs the package as if you had pulled it from CRAN or github
  devtools::install("KRLS2")
}

#----------
# Showing it works
#----------

#----------
# KRLS
#----------

## Careful if you load both KRLS and KRLS2, they will override the main function
## for now whenever I want krls from the old package I use KRLS::krls() and never
## explicitly load the old package
library(KRLS2)
n=50
betas=as.matrix(c(1,0, -2, 3))
p=length(betas)

X <- matrix(rnorm(n*p, mean = 2, sd = 3), ncol = p)
y <- rbinom(n, 1, prob = 1/(1+exp(-X%*%betas + X[, 1]/sqrt(abs(X[, 2])) + sin(X[, 3]))))

# Truncate false by default, outcome treated as continuous by default
krls.out <- krls(X=X, y=y)
krls.out.old <- KRLS::krls(X=X, y=y)

# Compare the coeffs rotated back to columns of full K to old coeffs
cbind(krls.out$coeffs, krls.out.old$coeffs)
# Not completely identical, shuttling back and forth to C++ makes sense
table(krls.out$coeffs == krls.out.old$coeffs)
sum(abs(krls.out$coeffs - krls.out.old$coeffs)) # Not so worried

# Cool
c(krls.out$var.avgderiv)
c(krls.out.old$var.avgderivatives)


# What aboud speed?
#install.packages("profvis")
n=500
betas=as.matrix(c(1,0, -2, 3))
p=length(betas)

X <- matrix(rnorm(n*p, mean = 2, sd = 3), ncol = p)
y <- rbinom(n, 1, prob = 1/(1+exp(-X%*%betas + X[, 1]/sqrt(abs(X[, 2])) + sin(X[, 3]))))
profvis::profvis({kout <- krls(X=X, y=y)})
profvis::profvis({kout <- KRLS::krls(X=X, y=y)})
## Generally about half of the time or less (without truncation!)
## Still have to add a couple features, but looking at about half the time needed

#----------
# KRLogit
#----------
# Wont work without truncate right now
ntrain <- 400
Xtrain <- X[1:ntrain, ]
ytrain <- y[1:ntrain]
Xtest <- X[(ntrain+1):nrow(X), ]
ytest <- y[(ntrain+1):length(y)]

krlogit.out=krls(X=Xtrain,y=ytrain, epsilon = .01, lambdarange = seq(0.001, 0.2, by = 0.02), method = "logistic", truncate = T)
# use old krls below bc haven't updated predict fn yet
krls.out=KRLS::krls(X=Xtrain, y=ytrain, derivative = F, vcov = F)
predkrlog <- predict(krlogit.out, Xtest)$fit
predkrls <- KRLS::predict.krls(krls.out, Xtest)$fit
sqrt(sum((ytest - predkrlog)^2))
sqrt(sum((ytest - predkrls)^2))
# did some sims, we are outperforming OOS pretty consistently in this toy example
