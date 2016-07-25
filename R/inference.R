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
  hessian <- NULL
  vcov.cb0 <- NULL
  vcov.cb0.trunc <- NULL
  
    
  if (vcov) {
    
    if (obj$loss == "leastsquares") {
      if(!sandwich) {
        sigmasq <- as.vector((1/n) * crossprod(y-yfitted))
        
        
        if(!obj$truncate) {
          vcov.c <- tcrossprod(mult_diag(obj$eigobj$vectors,sigmasq*(obj$eigobj$values+obj$lambda)^-2),obj$eigobj$vectors)
        } else {
          vcov.c <- tcrossprod(mult_diag(obj$Utrunc,sigmasq*(obj$D+obj$lambda)^-2),obj$Utrunc)
        }
      } else {
        ## get reconstructed K
        hessian <- krls_hess_sample(obj$K, obj$lambda)
        
        Ainv <- solve(hessian)
        
        score <- matrix(nrow = n, ncol = length(obj$coeffs))
        B <- matrix(0, nrow = length(obj$coeffs), ncol = length(obj$coeffs))
        
        ## todo: doing this in cpp or with matrices could be faster
        for(i in 1:n) {
          score[i, ] <- krls_gr(obj$K[, i, drop = F], y[i], yfitted[i], yfitted, obj$lambda/n)
          B <- B + tcrossprod(score[i, ])
        }
        
        if(!is.null(clusters)) {
          B <- matrix(0, nrow = length(obj$coeffs), ncol = length(obj$coeffs))
          for(j in 1:length(clusters)){
            B <- B + tcrossprod(apply(score[clusters[[j]], ], 2, sum))
          }
        }
        vcov.c <- Ainv %*% B %*% Ainv
      }
      
    } else {
      ## Vcov of d
      if (vcov & obj$truncate) {
        
        UDinv <- mult_diag(obj$Utrunc, 1/obj$D)
        
        hessian <- krlogit_hess_trunc(c(obj$chat, obj$beta0hat), obj$Utrunc, obj$D, y, obj$lambda)
        vcov.cb0 <- solve(hessian)
        
        score <- matrix(nrow = n, ncol = length(c(obj$chat, obj$beta0hat)))
        B <- matrix(0, nrow = length(c(obj$chat,obj$beta0hat)), ncol = length(c(obj$chat,obj$beta0hat)))
        
        for(i in 1:n) {
          score[i, ] = t(-1 * krlogit_gr_trunc(c(obj$chat, obj$beta0hat), obj$Utrunc[i, , drop=F], obj$D, y[i, drop = F], obj$lambda/n))
          B <- B + tcrossprod(score[i, ])
        }
        if(sandwich) {
          if(!is.null(clusters)) {
            B <- matrix(0, nrow = length(c(obj$chat,obj$beta0hat)), ncol = length(c(obj$chat,obj$beta0hat)))
            for(j in 1:length(clusters)){
              B <- B + tcrossprod(apply(score[clusters[[j]], ], 2, sum))
            }
            
          }
          vcov.cb0 <- vcov.cb0 %*% B %*% vcov.cb0
          
        }
        vcov.cb0.trunc <- vcov.cb0
        vcov.c <- tcrossprod(UDinv%*%vcov.cb0[1:length(obj$chat), 1:length(obj$chat)], UDinv)
        ## Just work with truncation later
      } else {vcov.cb0 = NULL}
      
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
    
    #construct coefhat=c for no truncation and coefhat = Utrunc*c
    #to take the truncation into the coefficients. 
    
    if(cpp) {
      
      if(!obj$truncate){
        derivout <- pwmfx(obj$K, X, obj$coeffs, vcov.c, tau, obj$sigma)
      } else {
        derivout <- pwmfx(tcrossprod(mult_diag(obj$Utrunc, obj$D), obj$Utrunc), X, obj$coeffs, vcov.c, tau, obj$sigma)
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
        derivatives[,k] <- tau*(-1/obj$sigma)*(L%*%obj$coeffs)
        if(obj$truncate==FALSE) {
          var.avgderivatives = NULL
        } else {
          var.avgderivatives[1,k] <- 1/(obj$sigma * n)^2 * sum(crossprod(tau^2, crossprod(L,vcov.c%*%L)))
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
              vcov.cb0.trunc = vcov.cb0.trunc,
            derivatives = derivatives,
            avgderivatives = avgderivatives,
            var.avgderivatives = var.avgderivatives)
  )
  
  class(z) <- "krls2"
  
  z <- fdskrls(z)
  
  if(returnmoreinf) {
    z$score = score
    z$hessian = hessian
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
    derivatives[,k] <- tau*(-1/sigma)*(L%*%coefhat)
    if(truncate==FALSE) {
      var.avgderivatives = NULL
    } else {
      var.avgderivatives[1,k] <- 1/(sigma * n)^2 * sum(crossprod(tau^2, crossprod(L,vcov.c%*%L)))
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
        if(object$loss == "leastsquares") {
          # SE (multiply by sqrt2 to correct for using data twice )
          se[1,i] <- as.vector(sqrt(t(h)%*%pout$vcov.fit%*%h))*sqrt(2)
        } else {
          #deriv.avgfd.logit <- colMeans(pout$deriv.logit[1:n, ] - pout$deriv.logit[(n+1):(2*n), ])
          deriv.avgfd.logit <- crossprod(h, pout$deriv.logit)
          vcov.avgfd <- tcrossprod(deriv.avgfd.logit %*% object$vcov.cb0.trunc, deriv.avgfd.logit)
          se[i,1] <- as.vector(sqrt(vcov.avgfd)) *sqrt(2)
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