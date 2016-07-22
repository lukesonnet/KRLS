## This file contains functions for predicting outcome variables
## Functions:
##            predict.krls2
##            logistic


## a predict function for class 'krlogit'
#' @export
predict.krls2 <- function(object, newdata, se.fit = FALSE, ...) {
  if (class(object) != "krls2") {
    warning("Object not of class 'krls2'")
    UseMethod("predict")
    return(invisible(NULL))
  }
  
  newdata <- as.matrix(newdata)
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted krls object")
  }
  
  newdataK <- newKernel(X = object$X, newData = newdata, whichkernel = object$kernel)
  
  #if(object$truncate){
  #  newdataK <- newdataK%*%object$Utrunc
  #}
  
  if(object$loss == "logistic") {
    yfitted <- logistic(K = newdataK, coeff = object$coeffs, beta0 = object$beta0hat)
    vcov.fit <- se.fit <- NULL
  } else if (object$loss == "leastsquares") {
    
    yfitted <- newdataK %*% object$coeffs
    
    # ses for fitted   
    if(se.fit){
      # transform to variance of c's on standarized scale
      vcov.c.raw <-  object$vcov.c * as.vector((1/var(object$y)))
      vcov.fitted <- tcrossprod(newdataK%*%vcov.c.raw,newdataK)          
      vcov.fit <- (apply(object$y,2,sd)^2)*vcov.fitted
      se.fit <- matrix(sqrt(diag(vcov.fit)),ncol=1)
    } else {
      vcov.fit <- se.fit <- NULL
    }
    
    # bring back to original scale
    yfitted <- (yfitted * apply(object$y,2,sd))+mean(object$y)
    
  }
  return(list(fit = yfitted,
              se.fit = se.fit, vcov.fit = vcov.fit,# newdata = newdata, 
              newdataK = newdataK))
}

## The logistic function that takes values for coeff, b0, and a K or Ktilde
#' @export
logistic <- function(K, coeff, beta0) {
  
  yhat <- 1 / (1 + exp(-(beta0+K%*%coeff)))
  
  return(yhat)
}
