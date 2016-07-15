## This file contains functions that handle lambda and sigma searches

## This function governs searching for lambda with a fixed sigma
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
      lambda <- lambdarange[which.min(lambdaMSE)]
    }
  }
  
  return(lambda)
  
}