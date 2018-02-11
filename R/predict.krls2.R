## This file contains functions for predicting outcome variables
## Functions:
##            predict.krls2
##            logistic


## a predict function for class 'krls2'
#' @export
predict.krls2 <- function(object, newdata, se.fit = FALSE, ...) {
  
  if(is.null(object$vcov.c) & se.fit)
    stop("se.fit requires the vcov.c. Use summary object with predict if you want stand error for the fits")
    # todo: we should probably call inference here if they want it?
  
  newdata <- as.matrix(newdata)
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted krls object")
  }
  
  newdataK <- newKernel(X = object$X, newData = newdata, whichkernel = object$kernel, b = object$b)

  if (object$loss == "logistic") {
    fit <- predict_logistic(
      newdataK = newdataK,
      dhat = object$dhat,
      coeffs = object$coeffs,
      beta0hat = object$beta0hat,
      U = object$U,
      D = object$D,
      vcov.d = object$vcov.d,
      se.fit = se.fit
    )
  } else if (object$loss == "leastsquares") {
    fit <- predict_leastsquares(
      newdataK = newdataK,
      y = object$y,
      coeffs = object$coeffs,
      vcov.c = object$vcov.c,
      se.fit = se.fit
    )
  }
  
  return(fit)
}

predict_logistic <- function(newdataK, dhat, coeffs, beta0hat, U, D, vcov.d, se.fit) {
  
  yfitted <- logistic(K = newdataK, coeff = coeffs, beta0 = beta0hat)
  
  if(se.fit){
  
    newU <- newdataK %*% mult_diag(U, 1 / D)
    # todo: could move to cpp if slow
  
    partiallogit <- partial_logit(newU, dhat, beta0hat)
    # each row is dy/dd with dy/db in the last column
    deriv.logit <- cbind(t(mult_diag(t(newU), partiallogit)),
                         partiallogit)
    
    vcov.fit <- tcrossprod(deriv.logit %*% vcov.d, deriv.logit)
    se.fit <- matrix(sqrt(diag(vcov.fit)),ncol=1)
  } else {
    vcov.fit <- se.fit <- deriv.logit <- NULL
  }

  fit <- list(
    fit = yfitted,
    se.fit = se.fit,
    vcov.fit = vcov.fit,
    newdataK = newdataK,
    deriv.logit = deriv.logit
  )
}
  
predict_leastsquares <- function(newdataK, coeffs, vcov.c, y, se.fit) {
  
  yfitted <- newdataK %*% coeffs
  
  # ses for fitted
  if (se.fit) {
    # transform to variance of c's on standarized scale
    vcov.c.raw <-  vcov.c * as.vector((1 / var(y)))
    vcov.fitted <-
      tcrossprod(newdataK %*% vcov.c.raw, newdataK)
    vcov.fit <- (apply(y, 2, sd) ^ 2) * vcov.fitted
    se.fit <- matrix(sqrt(diag(vcov.fit)), ncol = 1)
  } else {
    vcov.fit <- se.fit <- NULL
  }
  
  # bring back to original scale
  yfitted <- (yfitted * apply(y, 2, sd)) + mean(y)

  fit <- list(
    fit = yfitted,
    se.fit = se.fit,
    vcov.fit = vcov.fit,
    newdataK = newdataK
  )

}
