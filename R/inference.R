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
                            derivative = TRUE,
                            vcov = TRUE,
                            robust = FALSE,
                            clusters = NULL,
                            return_components = FALSE,
                            cpp = TRUE) {
  
  if (obj$kernel != 'gaussian') {
    stop(
      "Currently only supportmarginal effects and standard errors for ",
      "gaussian kernels"
    )
  }
  
  clustered <- !is.null(clusters)
  if (clustered) {
    robust <- TRUE
  }
  
  if (!vcov) {
    warning("Standard errors only available if `vcov = TRUE`")
  }
  
  if (!obj$truncate) {
    if (obj$loss == 'logistic' && (vcov || robust)) {
      warning(
        "Must use truncation to get standard errors with logistic loss; ",
        "returning average marginal effects without standard errors."
      )
      vcov <- FALSE
    } else if (obj$loss == "leastsquares" && robust) {
      warning(
        "Must use truncation to get robust standard errors with leastsquares ",
        "loss; returning average marginal effects without standard errors."
      )
      vcov <- FALSE
      robust <- FALSE
      cluster <- NULL
    }
  }
  
  weight <- length(unique(obj$w)) > 1

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
  
  if (clustered && !is.list(clusters)) {
    if (length(clusters) != length(obj$y)) {
      stop("`clusters` must be a vector the same length as `y`")
    }
    
    clusters <- lapply(
      unique(clusters), 
      function(clust) which(clusters == clust)
    )
  }
  
  ###----------------------------------------
  ## Getting vcov matrices
  ###----------------------------------------
  
  vcov.c <- NULL
  vcov.d <- NULL
  score <- NULL
  invhessian <- NULL

  if (robust && !clustered) {
    clusters <- seq_len(n)
  }
  
  if (vcov) {
    
    if (obj$loss == "leastsquares") {
      
      if (weight || robust) {
        invhessian <- krls_hess_trunc_inv(obj$U, obj$D, obj$w, obj$lambda)
      } 
      
      if (robust) {
        
        # actually t(score)
        score <- sapply(clusters, function(clust_ids) {
          krls_gr_trunc(
            obj$U[clust_ids, , drop = FALSE],
            obj$D,
            y[clust_ids],
            obj$w[clust_ids],
            yfitted[clust_ids],
            obj$dhat,
            obj$lambda / n
          )
        })
        
        vcov.d <- invhessian %*% tcrossprod(score) %*% invhessian
        
      } else {
        sigmasq <- as.vector((1 / n) * crossprod(y - yfitted))
        
        if (weight) {
          vcov.d <- 
            invhessian %*% 
            crossprod(mult_diag(obj$U, sigmasq * obj$w^2), obj$U) %*% 
            invhessian
        } else {
          vcov.c <- 
            tcrossprod(
              mult_diag(obj$U, sigmasq * (obj$D + obj$lambda)^-2),
              obj$U
            )
        }
      }
      
      if (is.null(vcov.c)) {
        UDinv <- mult_diag(obj$U, 1/obj$D)
        
        vcov.c <- tcrossprod(UDinv %*% vcov.d, UDinv)
      }
      
    } else if (obj$truncate) { # if loss == 'logistic'
      
      invhessian <- krlogit_hess_trunc_inv(c(obj$dhat, obj$beta0hat), obj$U, obj$D, y, obj$w, obj$lambda)

      if (robust) {
        
        # actually t(score)
        score <- sapply(clusters, function(clust_ids) {
          score_lambda <- length(clust_ids) * obj$lambda / n
          krlogit_gr_trunc(
            c(obj$dhat, obj$beta0hat),
            obj$U[clust_ids, , drop = FALSE],
            obj$D,
            y[clust_ids, drop = FALSE],
            obj$w[clust_ids],
            score_lambda
          ) * -1
        })
        
        
        vcov.d <- invhessian %*% tcrossprod(score) %*% invhessian
        
      } else {
        vcov.d <- invhessian
      }
      
      UDinv <- mult_diag(obj$U, 1/obj$D)
      # vcov.d has b0 in the last column and row for logistic loss
      # Need to drop for vcov.c
      n_d <- length(obj$dhat)
      vcov.c <- tcrossprod(
        UDinv %*% vcov.d[-(n_d + 1), -(n_d + 1), drop = FALSE],
        UDinv
      )
    }
  }

  ###----------------------------------------
  ## Getting pwmfx
  ###----------------------------------------
  
  if (derivative) {
    
    derivatives <- matrix(
      NA, 
      ncol=d,
      nrow=n,
      dimnames=list(NULL, colnames(X))
    )
    
    var.avgderivatives <- rep(NA, d)
    
    binaryindicator <- 
      apply(obj$X, 2, function(x) length(unique(x))) == 2
    
    if (any(!binaryindicator)) {
      
      if(obj$loss == "leastsquares") {
        tau <- rep(2, n)
      } else if (obj$loss == "logistic") {
        tau <- 2 * yfitted*(1-yfitted)
      }
      tau2 <- tau^2
      
      #construct coefhat=c for no truncation and coefhat = U*c
      #to take the truncation into the coefficients. 
      #Kfull rather than K here as truncation handled through coefs
      derivout <- as.matrix(apply(
        X[, !binaryindicator, drop = FALSE],
        2, 
        function(x) pwmfx(obj$K, x, obj$coeffs, vcov.c, tau, tau2, obj$b)
      ))

      derivatives[, !binaryindicator] <- derivout[1:n, ]
      var.avgderivatives[!binaryindicator] <- derivout[n+1, ]
      
      ## Rescale quantities of interest
      if (obj$loss == "leastsquares") {
        avgderivatives <- colMeans(as.matrix(derivatives))
        derivatives <- scale(y.init.sd*derivatives,center=FALSE,scale=X.init.sd)
        attr(derivatives,"scaled:scale")<- NULL
        avgderivatives <- as.vector(scale(y.init.sd*matrix(avgderivatives, nrow = 1),center=FALSE,scale=X.init.sd))
        attr(avgderivatives,"scaled:scale")<- NULL
        var.avgderivatives <- (y.init.sd/X.init.sd)^2*var.avgderivatives
        attr(var.avgderivatives,"scaled:scale")<- NULL
        

        
      } else {
        derivatives <- scale(derivatives, center=FALSE, scale = X.init.sd)
        avgderivatives <- as.vector(colMeans(as.matrix(derivatives)))
        var.avgderivatives <- (1/X.init.sd)^2*var.avgderivatives
      }
    }
    
    if (obj$loss == "leastsquares" && !is.null(vcov.c)) {
      vcov.c <- vcov.c * (y.init.sd^2)
    }
    
    if (any(binaryindicator)) {
      # Contrast vector
      h <- rep(c(1/n, -(1/n)), each=n)
      
      fdout <- as.matrix(sapply(
        which(binaryindicator),
        function(p) {
          firstdiffs(
            object = obj,
            n = n,
            p = p,
            h = h,
            vcov.c = vcov.c,
            vcov.d = vcov.d
          )
        }
      ))

      derivatives[, binaryindicator] <- fdout[1:n, ]
      avgderivatives[binaryindicator] <- fdout[n+2, ]
      var.avgderivatives[binaryindicator] <- fdout[n+1, ]
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
  names(avgderivatives) <- colnames(obj$X)
  names(var.avgderivatives) <- colnames(obj$X)
  
  
  z <- append(
    obj,
    list(
      vcov.c = vcov.c,
      vcov.d = vcov.d,
      derivatives = derivatives,
      avgderivatives = avgderivatives,
      var.avgderivatives = var.avgderivatives
    )
  )
  
  class(z) <- "krls2"

  if(return_components) {
    z$score <- score
    z$invhessian <- invhessian
  }
  
  return(z)
}


# First differences
firstdiffs <- function(object, n, p, h, vcov.c, vcov.d){ 
  
  n <- nrow(object$X)
  
  # compute marginal differences from min to max 
  Xall <- rbind(object$X, object$X)
  Xall[1:n, p] <- max(object$X[, p])
  Xall[(n+1):(2*n), p] <- min(object$X[, p])

  getvar <- !is.null(vcov.c)
  
  newdataK <- newKernel(X = object$X, newData = Xall, whichkernel = object$kernel, b = object$b)
  if (object$loss == "leastsquares") {
    pout <- predict_leastsquares(
      newdataK = newdataK,
      coeffs = object$coeffs,
      vcov.c = vcov.c,
      y = object$y,
      se.fit = getvar
    )
    
    if (getvar) {
      # multiply by 2 to correct for using data twice
      var <- as.vector(t(h) %*% pout$vcov.fit %*% h) * 2
    } else {
      var <- NA
    }
    
  } else {
    pout <- predict_logistic(
      newdataK = newdataK,
      dhat = object$dhat,
      coeffs = object$coeffs,
      beta0hat = object$beta0hat,
      U = object$U,
      D = object$D,
      vcov.d = vcov.d,
      se.fit = getvar
    )
    
    if (getvar) {
      deriv.avgfd.logit <- crossprod(h, pout$deriv.logit)
      vcov.avgfd <-
        tcrossprod(deriv.avgfd.logit %*% vcov.d, deriv.avgfd.logit)
      # multiply by 2 to correct for using data twice
      var <- as.vector(vcov.avgfd) * 2
    } else {
      var <- NA
    }
  }
  
  # store FD estimates
  est <- t(h) %*% pout$fit
  
  # all
  diffs <- pout$fit[1:n] - pout$fit[(n + 1):(2 * n)]
  
  fd <- c(diffs, var, est)
  
  return(fd)
}