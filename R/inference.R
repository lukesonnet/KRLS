## This file contains functions that generate the variance covariance matrices
## derivatives, hessians and more
## Functions:
##            inference.krls2
##            fdskrls
## Also in CPP:
##            See krlogit_fn, krlogit_fn_trunc, krlogit_hess_trunc etc...
##            Krlogit.fn, krlogit.hess.trunc, etc. are deprecated in
##            favor of these alternatives and can be found in deprecated/deprecated.R

#' @export
inference.krls2 <- function(obj,
                            sandwich = ifelse(obj$loss == "leastsquares", F, T),
                            robust = FALSE,
                            clusters = NULL,
                            returnmoreinf = FALSE,
                            vcov = TRUE,
                            derivative = TRUE,
                            cpp = TRUE) {
  
  if(derivative & !vcov) {
    stop("Derivatives require vcov = T")
  }
  
  if(!sandwich & !is.null(clusters)) {
    stop("Clusters require using sandwich == T")
  }
  
  if(!obj$truncate & obj$loss == "logistic") {
    stop("Derivatives only available for binary models with truncation")
  }
  
  if(!is.null(clusters) & !is.list(clusters)) {
    if(length(clusters) != nrow(obj$X)) stop("Clusters must be a vector the same length as X")
    clusters <- lapply(unique(clusters), function(clust) which(clusters == clust))
  } else if(is.null(clusters) & robust) {
    clusters <- as.list(1:nrow(obj$X))
  }
  
  ## Scale data
  y.init <- obj$y
  X.init <- obj$X
  X.init.sd <- apply(X.init, 2, sd)
  X <- scale(X.init, center = TRUE, scale = X.init.sd)    
  
  if (obj$loss == "leastsquares") {
    y.init.sd <- apply(y.init,2,sd)
    y.init.mean <- mean(y.init)
    y <- scale(y.init,center=y.init.mean,scale=y.init.sd)
    
    yfitted <- (obj$fitted - y.init.mean) / y.init.sd
    
  } else {
    yfitted <- obj$fitted
    y <- obj$y
  }
  
  n <- length(y.init)
  d <- ncol(X.init)
  
  ###----------------------------------------
  ## Getting vcov matrices
  ###----------------------------------------
  
  vcov.c <- NULL
  score <- NULL
  invhessian <- NULL
  vcov.db0 <- NULL
    
  if (vcov) {
    
    if (obj$loss == "leastsquares") {
      if(!sandwich) {
        sigmasq <- as.vector((1/n) * crossprod(y-yfitted))
        
        vcov.c <- tcrossprod(mult_diag(obj$U,sigmasq*(obj$D+obj$lambda)^-2),obj$U)
      } else {
        ## get reconstructed K
        if(!obj$truncate){
          stop("Sandwich estimators only available with truncated KRLS.")
          #invhessian <- krls_hess_inv(obj$K, obj$lambda)
        }
        else {
          UDinv <- mult_diag(obj$U, 1/obj$D)
          invhessian <- krls_hess_trunc_inv(obj$U, obj$D, obj$lambda)
        }
        
        ## todo: should we implement straight up hessian inversion? no... then we should force sandwich when logistic
        
        if(is.null(clusters)) {
          score <- matrix(nrow = n, ncol = length(obj$dhat))
          for(i in 1:n) {
            score[i, ] = krls_gr_trunc(obj$U[i, , drop = F], obj$D, y[i], yfitted[i], obj$dhat, obj$lambda/n)
          }
        } else {
          score <- matrix(nrow = length(clusters), ncol = length(obj$dhat))
          for(j in 1:length(clusters)){
            score[j, ] = krls_gr_trunc(obj$U[clusters[[j]], , drop = F], obj$D, y[clusters[[j]]], yfitted[clusters[[j]]], obj$dhat, length(clusters[[j]]) * obj$lambda/n)
          }
          
        }
        
        vcov.d <- invhessian %*% crossprod(score) %*% invhessian
        vcov.c <- tcrossprod(UDinv%*%vcov.d, UDinv)
        
      }
      
    } else { # if loss == 'logistic'
      ## Vcov of d
      if (vcov & obj$truncate) {
        
        UDinv <- mult_diag(obj$U, 1/obj$D)
        
        invhessian <- krlogit_hess_trunc_inv(c(obj$dhat, obj$beta0hat), obj$U, obj$D, y, obj$w, obj$lambda)
 
        if(sandwich) {
          if(is.null(clusters)) {
            score <- matrix(nrow = n, ncol = length(c(obj$dhat, obj$beta0hat)))
            for(i in 1:n) {
              score[i, ] = t(-1 * krlogit_gr_trunc(c(obj$dhat, obj$beta0hat), obj$U[i, , drop=F], obj$D, y[i, drop = F], obj$w[i], obj$lambda/n))
            }
          } else {
            score <- matrix(nrow = length(clusters), ncol = length(c(obj$dhat, obj$beta0hat)))
            for(j in 1:length(clusters)){
              score[j, ] <-  t(-1 * krlogit_gr_trunc(c(obj$dhat, obj$beta0hat), obj$U[clusters[[j]], , drop=F], obj$D, y[clusters[[j]], drop = F], obj$w[clusters[[j]]], length(clusters[[j]]) * obj$lambda/n))
            }
            
          }
          
          vcov.db0 <- invhessian %*% crossprod(score) %*% invhessian
          
        } else {
          vcov.db0 <- invhessian
        }
        vcov.c <- tcrossprod(UDinv%*%vcov.db0[1:length(obj$dhat), 1:length(obj$dhat)], UDinv)
        ## Just work with truncation later
      } else {vcov.db0 = NULL}
      
    }
    
  }
  
  ###----------------------------------------
  ## Getting pwmfx
  ###----------------------------------------
  
  if (derivative==TRUE) {
    
    obj$binaryindicator=matrix(FALSE,1,d)
    colnames(obj$binaryindicator) <- colnames(X)
    
    derivatives <- matrix(NA, ncol=d, nrow=n,
                          dimnames=list(NULL, colnames(X)))
    var.avgderivatives <- matrix(NA,1,d)
    
    if(obj$loss == "leastsquares") {
      tau <- rep(2, n)
    } else if (obj$loss == "logistic") {
      tau <- 2 * yfitted*(1-yfitted)
    }
    
    #construct coefhat=c for no truncation and coefhat = U*c
    #to take the truncation into the coefficients. 
    
    if(cpp) {
      
      if(!obj$truncate){
        derivout <- pwmfx(obj$K, X, obj$coeffs, vcov.c, tau, obj$b)
      } else {
        derivout <- pwmfx(tcrossprod(mult_diag(obj$U, obj$D), obj$U), X, obj$coeffs, vcov.c, tau, obj$b)
      }
      derivatives <- derivout[1:n, ]
      var.avgderivatives <- derivout[n+1, ]
      
    } else {
      rows <- cbind(rep(1:n, each = n), 1:n)
      distances <- X[rows[,1],] - X[ rows[,2],] 
      
      for(k in 1:d){
        print(paste("Computing derivatives, variable", k))
        if(d==1){
          distk <-  matrix(distances,n,n,byrow=TRUE)
        } else {
          distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
        }
        L <-  distk*obj$K  #Kfull rather than K here as truncation handled through coefs
        derivatives[,k] <- tau*(-1/obj$b)*(L%*%obj$coeffs)
        if(obj$truncate==FALSE) {
          var.avgderivatives = NULL
        } else {
          var.avgderivatives[1,k] <- 1/(obj$b * n)^2 * sum(crossprod(tau^2, crossprod(L,vcov.c%*%L)))
        }
      }
    }
    
    ## Rescale quantities of interest
    if (obj$loss == "leastsquares") {
      avgderivatives <- colMeans(as.matrix(derivatives))
      derivatives <- scale(y.init.sd*derivatives,center=FALSE,scale=X.init.sd)
      attr(derivatives,"scaled:scale")<- NULL
      avgderivatives <- scale(y.init.sd*matrix(avgderivatives, nrow = 1),center=FALSE,scale=X.init.sd)
      attr(avgderivatives,"scaled:scale")<- NULL
      var.avgderivatives <- (y.init.sd/X.init.sd)^2*var.avgderivatives
      attr(var.avgderivatives,"scaled:scale")<- NULL
      
      vcov.c <- vcov.c * (y.init.sd^2)
      
    } else {
      derivatives <- scale(derivatives, center=F, scale = X.init.sd)
      avgderivatives <- t(colMeans(as.matrix(derivatives)))
      var.avgderivatives <- (1/X.init.sd)^2*var.avgderivatives
    }
    
    
  } else {
    derivatives=NULL
    avgderivatives=NULL
    var.avgderivatives=NULL
  }
  
  
  
  ###----------------------------------------
  ## Returning results
  ###----------------------------------------
  
  colnames(derivatives) <- colnames(obj$X)
  colnames(avgderivatives) <- colnames(obj$X)
  names(var.avgderivatives) <- colnames(obj$X)
  
  
  z <- c(obj,
         list(vcov.c = vcov.c,
              vcov.db0 = vcov.db0,
            derivatives = derivatives,
            avgderivatives = avgderivatives,
            var.avgderivatives = var.avgderivatives)
  )
  
  class(z) <- "krls2"
  
  z <- fdskrls(z)
  
  if(returnmoreinf) {
    z$score = score
    z$invhessian = invhessian
  }
  
  return(z)
}

# Chad's pwmfxR function which is mostly cribbed from above
pwmfxR <- function() {
  rows <- cbind(rep(1:n, each = n), 1:n)
  distances <- X[rows[,1],] - X[ rows[,2],] 
  
  for(k in 1:d){
    print(paste("Computing derivatives, variable", k))
    if(d==1){
      distk <-  matrix(distances,n,n,byrow=TRUE)
    } else {
      distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
    }
    L <-  distk*Kdat$K  #Kfull rather than K here as truncation handled through coefs
    derivatives[,k] <- tau*(-1/b)*(L%*%coefhat)
    if(truncate==FALSE) {
      var.avgderivatives = NULL
    } else {
      var.avgderivatives[1,k] <- 1/(b * n)^2 * sum(crossprod(tau^2, crossprod(L,vcov.c%*%L)))
    }
  }
}

# First differences
#' @export
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
        if(object$loss == "leastsquares") {
          # SE (multiply by sqrt2 to correct for using data twice )
          se[1,i] <- as.vector(sqrt(t(h)%*%pout$vcov.fit%*%h))*sqrt(2)
        } else {
          #deriv.avgfd.logit <- colMeans(pout$deriv.logit[1:n, ] - pout$deriv.logit[(n+1):(2*n), ])
          deriv.avgfd.logit <- crossprod(h, pout$deriv.logit)
          vcov.avgfd <- tcrossprod(deriv.avgfd.logit %*% object$vcov.db0, deriv.avgfd.logit)
          se[1,i] <- as.vector(sqrt(vcov.avgfd)) *sqrt(2)
        }
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