# This file contains functions that generate the variance covariance matrix
#' @export
inference.krls2 <- function(obj,
                            sandwich = ifelse(obj$loss == "leastsquares", F, T),
                            clusters = NULL,
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
  
  if(!is.list(clusters)) {
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
  
  if (vcov) {
    
    if (obj$loss == "leastsquares") {
      if(!sandwich) {
        sigmasq <- as.vector((1/n) * crossprod(y-yfitted))
        
        vcov.c <- tcrossprod(mult_diag(obj$eigobj$vectors,sigmasq*(obj$eigobj$values+obj$lambda)^-2),obj$eigobj$vectors)
      } else {
        
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
        vcov.c <- tcrossprod(UDinv%*%vcov.cb0[1:length(obj$chat), 1:length(obj$chat)], UDinv)
        ## todo: return var b0
      } else {vcov.cb0 = NULL}
      
    }
    
  }
  
  ###----------------------------------------
  ## Getting pwmfx
  ###----------------------------------------
  
  if (derivative==TRUE) {
    
    derivatives <- matrix(NA, ncol=d, nrow=n,
                          dimnames=list(NULL, colnames(X)))
    var.avgderivatives <- matrix(NA,1,d)
    
    if(obj$loss == "leastsquares") {
      tau <- rep(2, n)
    } else if (obj$loss == "logistic") {
      tau <- yfitted*(1-yfitted)
    }
    
    #construct coefhat=c for no truncation and coefhat = Utrunc*c
    #to take the truncation into the coefficients. 
    
    if(cpp) {
      
      derivout <- pwmfx(obj$K, X, obj$coeffs, vcov.c, tau, obj$sigma)
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
      avgderivatives <- colMeans(derivatives)
      derivatives <- scale(y.init.sd*derivatives,center=FALSE,scale=X.init.sd)
      attr(derivatives,"scaled:scale")<- NULL
      avgderivatives <- scale(y.init.sd*matrix(avgderivatives, nrow = 1),center=FALSE,scale=X.init.sd)
      attr(avgderivatives,"scaled:scale")<- NULL
      var.avgderivatives <- (y.init.sd/X.init.sd)^2*var.avgderivatives
      attr(var.avgderivatives,"scaled:scale")<- NULL
    } else {
      derivatives <- scale(derivatives, center=F, scale = X.init.sd)
      avgderivatives <- t(colMeans(derivatives))
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
  
  z <- list(vcov.cb0 = vcov.cb0,
            vcov.c = vcov.c,
            score = score,
            hessian = hessian,
            derivatives = derivatives,
            avgderivatives = avgderivatives,
            var.avgderivatives = var.avgderivatives
  )
  
  return(z)
}