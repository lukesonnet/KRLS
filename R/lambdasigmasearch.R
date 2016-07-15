## This file contains functions that handle lambda and sigma searches

## This function searches for both lambda and sigma
#' @export
lambdasigmasearch <- function(y,
                              X.init,
                              hyperctrl,
                              control) {
  
  if(is.null(hyperctrl$lambdarange)) {
    if (!is.null(hyperctrl$sigmarange)) stop("Grid search for sigma only works with fixed lambda or lambda grid search.")
    
    fit.hyper <- optim(par=log(c(hyperctrl$lambdastart, 2.2*control$d)), lambdasigma.fn,
                       X=X.init, y=y, folds=length(hyperctrl$chunks),
                       chunks = hyperctrl$chunks, truncate=control$truncate,
                       epsilon=control$epsilon, lastkeeper=control$lastkeeper,
                       vcov = FALSE,
                       control=list(trace = T, abstol = 1e-4), method="BFGS")
    
    lambda <- exp(fit.hyper$par[1])
    sigma <- exp(fit.hyper$par[2])
  } else {
    if (hyperctrl$optimsigma) stop("Optimizing sigma does not work with a grid search for lambda.")
    hypergrid <- as.matrix(expand.grid(hyperctrl$lambdarange, hyperctrl$sigmarange))
    hyperMSE <- NULL
    for(i in 1:nrow(hypergrid)){
      hyperMSE[i] <- lambdasigma.fn(par = log(hypergrid[i, ]), X=X.init,
                                    y=y, folds=length(hyperctrl$chunks),
                                    chunks = hyperctrl$chunks, truncate=control$truncate,
                                    epsilon=control$epsilon, lastkeeper=control$lastkeeper)
    }
    hyper <- hypergrid[which.min(hyperMSE), ]
    lambda <- hyper[1]
    sigma <- hyper[2]
  }
  
  return(list(lambda=lambda,
              sigma = sigma))
  
}

## This function governs searching for lambda with a fixed sigma
#' @export
lambdasearch <- function(y,
                         X.init,
                         Kdat,
                         hyperctrl,
                         control,
                         sigma) {
  
  if (control$loss == "leastsquares") {
    if(control$truncate) {
      #todo: no way this lambdasearch is right. The Looe can't be right, the bounds can't be right... but it works
      #stop("Must specify lambda for truncated least squares for now.")
      lambda <- lambdasearch(y=y, D=kdat$eigvals, Utrunc=Kdat$Utrunc,
                             Eigenobject=Kdat$eigobj, truncate=control$truncate)#,eigtrunc=eigtrunc,noisy=noisy,L=L,U=U)
    } else {
      lambda <- lambdasearch(y=y, K=kdat$K, Eigenobject=kdat$eigobj,
                             truncate=control$truncate)#,eigtrunc=eigtrunc,noisy=noisy,L=L,U=U)
    }
  } else {
    
    if(is.null(hyperctrl$lambdarange)) {
      
      fit.lambda <- optim(par=log(hyperctrl$lambdastart), lambdasigma.fn,
                          X=X.init, y=y, folds=length(hyperctrl$chunks),
                          chunks = hyperctrl$chunks, truncate=control$truncate,
                          epsilon=control$epsilon, lastkeeper=control$lastkeeper,
                          sigma = sigma,
                          vcov = FALSE,
                          control=list(trace = T, abstol = 1e-4), method="BFGS")
      
      lambda <- exp(fit.lambda$par)
      
    } else {
      
      lambdaMSE <- NULL
      for(i in 1:length(lambdarange)){
        lambdaMSE[i] <- lambdasigma.fn(par = log(hyperctrl$lambdarange[i]), X=X.init,
                                       y=y, folds=hyperfolds, folds=length(hyperctrl$chunks),
                                       chunks = hyperctrl$chunks, truncate=control$truncate,
                                       epsilon=control$epsilon, lastkeeper=control$lastkeeper,
                                       sigma = sigma)
      }
      lambda <- hyperctrl$lambdarange[which.min(lambdaMSE)]
    }
  }
  
  return(lambda)
  
}

## This function searches for just sigma
#' @export
sigmasearch <- function(y,
                              X.init,
                              hyperctrl,
                              control,
                              lambda) {
  if(is.null(hyperctrl$sigmarange)) {
    fit.sigma <- optim(par=log(2.2*control$d), lambdasigma.fn,
                       X=X.init, y=y, folds=length(hyperctrl$chunks),
                       chunks = hyperctrl$chunks, truncate=control$truncate,
                       epsilon=control$epsilon, lastkeeper=control$lastkeeper,
                       lambda = lambda,
                       vcov = FALSE,
                       control=list(trace = T, abstol = 1e-4), method="BFGS")
    sigma <- exp(fit.sigma$par)
  } else {
    sigmaMSE <- NULL
    for(i in 1:length(hyperctrl$sigmarange)){
      if(hyperctrl$sigmarange[i] < 0) stop("All sigmas must be positive")
      sigmaMSE[i] <- lambdasigma.fn(par = log(hyperctrl$sigmarange[i]), X=X.init,
                                    y=y, folds=length(hyperctrl$chunks),
                                    chunks = hyperctrl$chunks, truncate=control$truncate,
                                    epsilon=control$epsilon, lastkeeper=control$lastkeeper,
                                    #sigma = sigma,
                                    lambda = lambda)
      sigma <- hyperctrl$sigmarange[which.min(sigmaMSE)]
    }
  }

  return(sigma)
  
}