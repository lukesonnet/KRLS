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
  
  if(is.null(object$vcov.c) & se.fit)
    stop("se.fit requires the vcov.c. Use summary object with predict if you want stand error for the fits")
    # todo: we should probably call inference here if they want it?
  
  newdata <- as.matrix(newdata)
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted krls object")
  }
  
  newdataK <- newKernel(X = object$X, newData = newdata, whichkernel = object$kernel, b = object$b)

  #if(object$truncate){
  #  newdataK <- newdataK%*%object$U
  #}
  
  if(object$loss == "logistic") {
    yfitted <- logistic(K = newdataK, coeff = object$coeffs, beta0 = object$beta0hat)
    
    if(se.fit){
      # use truncation so that we can use vcov.db0 because vcov.c does not have uncertainty of b0. Could change this and return vcov.c padded with UDinv %*% vcov.db0[1:nrow(UDinv), ncol(vcov.db0)] instead
      newU <- newdataK %*% mult_diag(object$U, 1/object$D)
      # todo: could move to cpp if slow
      partiallogit <- partial_logit(newU, object$dhat, object$beta0hat)
      # each row is dy/dd with dy/db in the last column
      deriv.logit <- cbind(t(mult_diag(t(newU), partiallogit)),
                           partiallogit)
      vcov.fit <- tcrossprod(deriv.logit %*% object$vcov.db0, deriv.logit)
      se.fit <- matrix(sqrt(diag(vcov.fit)),ncol=1)
    } else {
      vcov.fit <- se.fit <- deriv.logit <- NULL
    }
    
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
  fit <- list(fit = yfitted,
              se.fit = se.fit, vcov.fit = vcov.fit,
              newdataK = newdataK)
  
  if(object$loss == "logistic") {
    fit$deriv.logit <- deriv.logit
  }
  
  return(fit)
}

## The logistic function that takes values for coeff, b0, and a K or Ktilde
#' @export
logistic <- function(K, coeff, beta0) {
  
  yhat <- 1 / (1 + exp(-(beta0+K%*%coeff)))
  
  return(yhat)
}
