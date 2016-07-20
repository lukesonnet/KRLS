# This file contains functions that return hessians, scores, and more
fdskrls <- 
  function(object,...){ 
    d <- ncol(object$X)
    n <- nrow(object$X)
    lengthunique    <- function(x){length(unique(x))}
    
    fdderivatives           =object$derivatives
    fdavgderivatives        =object$avgderivatives
    fdvar.var.avgderivatives=object$var.avgderivatives    
    
    # vector with positions of binary variables
    binaryindicator <-which(apply(object$X,2,lengthunique)==2)    
    if(length(binaryindicator)==0){
      # no binary vars in X; return derivs as is  
    } else {
      print("found binary")
      # compute marginal differences from min to max 
      est  <- se <- matrix(NA,nrow=1,ncol=length(binaryindicator))
      diffsstore <- matrix(NA,nrow=n,ncol=length(binaryindicator))
      for(i in 1:length(binaryindicator)){
        X1 <- X0 <- object$X
        # test data with D=Max
        X1[,binaryindicator[i]] <- max(X1[,binaryindicator[i]])
        # test data with D=Min
        X0[,binaryindicator[i]] <- min(X0[,binaryindicator[i]])
        Xall      <- rbind(X1,X0)
        # contrast vector
        h         <- matrix(rep(c(1/n,-(1/n)),each=n),ncol=1)
        # fitted values
        pout      <- predict(object,newdata=Xall,se.fit=TRUE)
        # store FD estimates
        est[1,i] <- t(h)%*%pout$fit        
        # SE (multiply by sqrt2 to correct for using data twice )
        se[1,i] <- as.vector(sqrt(t(h)%*%pout$vcov.fit%*%h))*sqrt(2)
        # all
        diffs <- pout$fit[1:n]-pout$fit[(n+1):(2*n)]          
        diffsstore[,i] <- diffs 
      }        
      # sub in first differences
      object$derivatives[,binaryindicator] <- diffsstore
      object$avgderivatives[,binaryindicator] <- est
      object$var.avgderivatives[binaryindicator] <- se^2
      object$binaryindicator[,binaryindicator] <- TRUE
    }  
    
    return(invisible(object))
    
  }


## The functions below are not used because their cpp alternatives are preferred
## for speed

## The Kernel Regularized Logit target function to be minimized
## Parameters:
##   'par' - the parameters to be optimized, contains C and beta0, the unpenalized constant
##           because demeaning the data here does not make sense.
##   'K' - the Kernel matrix
##   'y' - the outcome variable
##   'lambda' - the regularizing parameter
## Values:
##   'r' - The penalized log-likelihood given 'y' and 'K'
#' @export
krlogit.fn <- function(par, K, y, lambda = 0.5) {
  
  coef <- par[1:ncol(K)]
  beta0 <- par[ncol(K)+1]
  
  Kc <- crossprod(K, coef)
  
  r <- sum(y * log(1 + exp( -(beta0 +  Kc))) + (1 - y) * log(1 + exp(beta0 + Kc))) +
    lambda*crossprod(Kc, coef)
  return(r)
}

## This has the wrong norm
#' @export
krlogit.fn.trunc <- function(par, K, Utrunc, y, lambda = 0.5) {
  
  coef <- par[1:ncol(K)]
  beta0 <- par[ncol(K)+1]
  
  Kc <- crossprod(t(K), coef)
  
  r <- sum(y * log(1 + exp( -(beta0 +  Kc))) + (1 - y) * log(1 + exp(beta0 + Kc))) +
    lambda*tcrossprod(coef, Utrunc)%*%Kc
  return(r)
}

## The analytic gradient. Pretty sure this is correct, was previously confirmed
## with the numDeriv package
#' @export
krlogit.gr <- function(par, K, y, lambda){
  
  coef <- par[1:ncol(K)]
  beta0 <- par[ncol(K)+1]
  
  Kc <- crossprod(K, coef)
  
  dc <- -crossprod(K, y - 1 / (1 + exp(-(beta0 + Kc)))) +
    2*lambda * Kc
  db0 <- -sum(y - 1 / (1 + exp(-(beta0 + Kc))))
  
  return(c(dc, db0))
}

## this has the wrong norm
#' @export
krlogit.gr.trunc <- function(par, K, Utrunc, y, lambda){
  
  coef <- par[1:ncol(K)]
  beta0 <- par[ncol(K)+1]
  
  Kc <- crossprod(t(K), coef)
  
  dc <- -crossprod(K, y - 1 / (1 + exp(-(beta0 + Kc)))) +
    2*lambda*crossprod(Utrunc, Kc)
  db0 <- -sum(y - 1 / (1 + exp(-(beta0 + Kc))))
  
  return(c(dc, db0))
}

#' @export
krlogit.hess.trunc <- function(par, Utrunc, D, y, lambda) {
  
  coef <- par[1:ncol(Utrunc)]
  beta0 <- par[ncol(Utrunc)+1]
  
  Ud <- Utrunc %*% coef
  
  meat <- exp(-Ud - beta0) / (1 + exp(-Ud - beta0))^2
  
  dcdc <- mult_diag(t(Utrunc), meat) %*% Utrunc + diag(2 * lambda / D)
  dcdb <- crossprod(Utrunc, meat)
  dbdb <- sum(meat)
  
  print(dcdc[1:3, 1:3])
  
  ret = matrix(nrow = length(par), ncol = length(par))
  ret[1:length(coef), 1:length(coef)] <- dcdc
  ret[length(par), 1:length(coef)] <- dcdb
  ret[1:length(coef), length(par)] <- dcdb
  ret[length(par), length(par)] <- dbdb
  
  print(ret[1:3, 1:3])
  
  return(ret)
}