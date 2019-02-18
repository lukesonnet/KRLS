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
  hessian <- NULL
  sigmasq_epsilon <- NULL

  if (robust && !clustered) {
    clusters <- seq_len(n)
  }
  
  if (vcov) {
    
    if (obj$loss == "leastsquares") {
      
      if (weight || robust) {
        hessian <- krls_hess_trunc(obj$U, obj$D, obj$w, obj$lambda)
        invhessian <- solve(hessian)
      } 
      
      if (robust) {
        
        # actually t(score)
        score <- do.call(
          cbind,
          lapply(clusters, function(clust_ids) {
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
        )
        
        vcov.d <- invhessian %*% tcrossprod(score) %*% invhessian
        
      } else {
        sigmasq_epsilon <- as.vector((1 / n) * crossprod(y - yfitted))

        if (obj$truncate) {
          # With truncation, we cannot "flip the middle",
          # we instead use the Woodbury matrix identity=
          
          if (weight) {
            stop("TODO")
          } else {
            Ulambdainv <- KRLS2:::mult_diag(obj$U, rep(1 / obj$lambda, nrow(obj$U)))
            
            UDUt_lambda_inv <- 
              diag(1 / obj$lambda, nrow(obj$U)) - 
              tcrossprod(
                Ulambdainv %*% solve(diag(1 / obj$D) + crossprod(Ulambdainv, obj$U)), 
                Ulambdainv
              )
  
            vcov.c <- sigmasq_epsilon * (UDUt_lambda_inv %*% UDUt_lambda_inv)
          }
        } else {
          if (weight) {
            stop("TODO")
          } else {
            # Without truncation or weights, we can "flip the middle"
            vcov.c <- 
              tcrossprod(
                mult_diag(obj$U, sigmasq_epsilon * (obj$D + obj$lambda)^-2),
                obj$U
              )
          }
        }
      }
      
      if (is.null(vcov.c)) {
        UDinv <- mult_diag(obj$U, 1/obj$D)
        
        vcov.c <- tcrossprod(UDinv %*% vcov.d, UDinv)
      }
      
    } else if (obj$truncate) { # if loss == 'logistic'
      
      hessian <- krlogit_hess_trunc(c(obj$dhat, obj$beta0hat), obj$U, obj$D, y, obj$w, obj$lambda)
      invhessian <- solve(hessian)

      if (robust) {
        
        # actually t(score)
        score <- do.call(
          cbind,
          lapply(clusters, function(clust_ids) {
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
        )
        
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
    
    derivatives <- var.derivatives <- matrix(
      NA, 
      ncol=d,
      nrow=n,
      dimnames=list(NULL, colnames(X))
    )
    
    var.avgderivatives <- rep(NA, d)
    
    binaryindicator <- 
      apply(obj$X, 2, function(x) length(unique(x))) == 2
    
    # Get actual pwmfx
    if (any(!binaryindicator)) {
      
      if (obj$loss == "leastsquares") {
        tau <- rep(1, n)
      } else if (obj$loss == "logistic") {
        tau <- yfitted * (1 - yfitted)
      }
      
      #construct coefhat=c for no truncation and coefhat = U*c
      #to take the truncation into the coefficients. 
      #Kfull rather than K here as truncation handled through coefs
      for (i in seq_len(d)) {
        if (!binaryindicator[i]) {
          if (is.numeric(vcov.c)) {
            deriv_list <- pwmfx(
              obj$K, 
              X[, i, drop = FALSE], 
              obj$coeffs, 
              vcov.c,
              tau, 
              obj$b
            )
            var.derivatives[, i] <- diag(deriv_list$var_deriv)
            var.avgderivatives[i] <- deriv_list$var_avg_deriv
          } else {
            deriv_list <- pwmfx(
              obj$K, 
              X[, i, drop = FALSE], 
              obj$coeffs,
              matrix(0),
              tau, 
              obj$b
            )
          }
          
          derivatives[, i] <- deriv_list$deriv
        }
      }
      
      ## Rescale quantities of interest
      if (obj$loss == "leastsquares") {
        
        avgderivatives <- colMeans(as.matrix(derivatives))
        derivatives <- scale(y.init.sd*derivatives,center=FALSE,scale=X.init.sd)
        attr(derivatives,"scaled:scale")<- NULL
        var.derivatives <- scale(y.init.sd^2 * var.derivatives, center=FALSE, scale=X.init.sd^2)
        attr(avgderivatives,"scaled:scale")<- NULL
        avgderivatives <- as.vector(
          scale(y.init.sd * matrix(avgderivatives, nrow = 1), center=FALSE, scale=X.init.sd)
        )
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
    
    # Get FDs
    if (any(binaryindicator)) {
      # Contrast 
      if (obj$loss == "leastsquares") {
        h <- as.matrix(rep(1/n, each=n))
      } else {
        h <- as.matrix(rep(c(1/n, -(1/n)), each=n))
      }
      
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
  colnames(var.derivatives) <- colnames(obj$X)
  names(avgderivatives) <- colnames(obj$X)
  names(var.avgderivatives) <- colnames(obj$X)
  
  
  z <- append(
    obj,
    list(
      vcov.c = vcov.c,
      vcov.d = vcov.d,
      derivatives = derivatives,
      avgderivatives = avgderivatives,
      var.derivatives = var.derivatives,
      var.avgderivatives = var.avgderivatives
    )
  )
  
  class(z) <- "krls2"

  if(return_components) {
    z$sigmasq_epsilon <- sigmasq_epsilon
    z$score <- score
    z$hessian <- hessian
    z$invhessian <- invhessian
  }
  
  return(z)
}

# First differences
firstdiffs <- function(object, n, p, h, vcov.c, vcov.d){ 
  
  n <- nrow(object$X)
  
  # compute marginal differences from min to max 
  X1 <- X0 <- object$X
  # Given we know it's binary
  X1[, p] <- max(object$X[, p])
  X0[, p] <- min(object$X[, p])
  
  getvar <- !is.null(vcov.c)
  
  if (object$loss == "leastsquares") {
    M <- 
      newKernel(X = object$X, newData = X1, whichkernel = object$kernel, b = object$b) -
      newKernel(X = object$X, newData = X0, whichkernel = object$kernel, b = object$b)
    
    diffs <- M %*% object$coeffs * sd(object$y)
    est <- crossprod(h, diffs)
    
    if (getvar) {
      # multiply by 2 to correct for using data twice
      
      var <- crossprod(h, M %*% vcov.c %*% crossprod(M, h))
                       
    } else {
      var <- NA
    }
    
  } else {
    newdataK <- newKernel(X = object$X, newData = rbind(X0, X1), whichkernel = object$kernel, b = object$b)
    
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
    
    est <- crossprod(h, pout$fit)
    diffs <- pout$fit[1:n] - pout$fit[(n + 1):(2 * n)]
    if (getvar) {
      deriv.avgfd.logit <- crossprod(h, pout$deriv.logit)
      vcov.avgfd <-
        tcrossprod(deriv.avgfd.logit %*% vcov.d, deriv.avgfd.logit)
      var <- as.vector(vcov.avgfd)
    } else {
      var <- NA
    }
  }
  
  fd <- c(diffs, var, est)
  
  return(fd)
}
