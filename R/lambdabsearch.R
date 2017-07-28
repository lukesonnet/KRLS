## This file contains functions that handle lambda and b searches
## Functions:
##            lambdabsearch
##            bsearch
##            lambdasearch
##            lambdaline (lambdasearch in the old KRLS)
##            looloss

## This function searches for both lambda and b
#' @export
lambdabsearch <- function(y,
                          X,
                          control) {

  if(is.null(control$lambdarange)) {
    if (!is.null(control$brange)) stop("Grid search for b only works with fixed lambda or lambda grid search.")

    fit.hyper <- optim(par=log(c(control$lambdastart, 2.2*control$d)), lambdab.fn,
                       X=X, y=y, folds=length(control$chunks),
                       chunks = control$chunks, ctrl = control,
                       vcov = FALSE,
                       control=list(trace = T, abstol = 1e-4), method="BFGS")

    lambda <- exp(fit.hyper$par[1])
    b <- exp(fit.hyper$par[2])
  } else {
    if (control$optimb) stop("Optimizing b does not work with a grid search for lambda.")
    hypergrid <- as.matrix(expand.grid(control$lambdarange, control$brange))
    hyperMSE <- NULL
    for(i in 1:nrow(hypergrid)){
      hyperMSE[i] <- lambdab.fn(par = log(hypergrid[i, ]), X=X,
                                    y=y, folds=length(control$chunks),
                                    chunks = control$chunks, ctrl = control)
    }
    hyper <- hypergrid[which.min(hyperMSE), ]
    lambda <- hyper[1]
    b <- hyper[2]
  }

  return(list(lambda=lambda,
              b = b))

}

## This function governs searching for lambda with a fixed b
#' @export
lambdasearch <- function(y,
                         X,
                         Kdat,
                         control,
                         b) {

  if (control$loss == "leastsquares") {
    if(control$weight){
      print('here')
      lambda <- lambdaline(y=y, D=Kdat$D, U=Kdat$U, w=control$w, tol=control$tol, noisy = !control$quiet)
    } else {
      lambda <- lambdaline(y=y, D=Kdat$D, U=Kdat$U, tol=control$tol, noisy = !control$quiet)
    }
  } else {

    if(is.null(control$lambdarange)) {
      
      start_val <- log(control$lambdastart)
        
      #fit.lambda <- optim(par=start_val, lambdab.fn,
      #                   Kdat = Kdat, y=y, folds=length(control$chunks),
      #                    chunks = control$chunks, ctrl = control,
      #                    b = b,
      #                    vcov = FALSE,
      #                    control=list(trace = T, abstol = 1e-4, REPORT = 1), method="Nelder-Mead")

      fit.lambda <- optimize(lambdab.fn, interval=log(c(10^-12, 4)),
                          Kdat = Kdat, y=y, folds=length(control$chunks),
                          chunks = control$chunks, ctrl = control,
                          b = b,
                          vcov = FALSE)
        
      fit.lambda$par=fit.lambda$minimum
  
        
      if(!control$quiet) {
        print(fit.lambda)
      }
        
      lambda <- exp(fit.lambda$par)
      
    } else {

      lambdaMSE <- NULL
      for(i in 1:length(control$lambdarange)){
        lambdaMSE[i] <- lambdab.fn(par = log(control$lambdarange[i]), Kdat = Kdat,
                                       y=y, folds=length(control$chunks),
                                       chunks = control$chunks, ctrl = control,
                                       b = b)
      }
      lambda <- control$lambdarange[which.min(lambdaMSE)]
    }
  }

  return(lambda)

}

## This function searches for just b
#' @export
bsearch <- function(y,
                              X,
                              control,
                              lambda) {
  if(is.null(control$brange)) {
    fit.b <- optim(par=log(2.2*control$d), lambdab.fn,
                       X=X, y=y, folds=length(control$chunks),
                       chunks = control$chunks, ctrl = control,
                       lambda = lambda,
                       vcov = FALSE,
                       control=list(trace = T, abstol = 1e-4), method="BFGS")
    b <- exp(fit.b$par)
  } else {
    bMSE <- NULL
    for(i in 1:length(control$brange)){
      if(control$brange[i] < 0) stop("All bs must be positive")
      bMSE[i] <- lambdab.fn(par = log(control$brange[i]), X=X,
                                    y=y, folds=length(control$chunks),
                                    chunks = control$chunks, ctrl = control,
                                    #b = b,
                                    lambda = lambda)
      b <- control$brange[which.min(bMSE)]
    }
  }

  return(b)

}

# todo: Unless the folds are much smaller, lastkeeper will probably be the same
# or only one or two smaller, should we pass lastkeeper and just use it?
## This computes the RMSPE using 'folds' cross-validation, fitting a new model
## for each subset of the data
#' @export
lambdab.fn <- function(par = NULL,
                           y = NULL,
                           X = NULL,
                           Kdat = NULL,
                           folds = NULL,
                           chunks = NULL,
                           b = NULL,
                           lambda = NULL,
                           vcov = NULL,
                           ctrl = NULL) {

  if(is.null(c(lambda, b))) {
    lambda <- exp(par[1])
    b <- exp(par[2])
    Kdat <- generateK(X=X,
                      b=b,
                      control=ctrl)
  } else {
    if (is.null(lambda)) {
      lambda <- exp(par[1])
    }
    if (is.null(b)) {
      b <- exp(par)
      Kdat <- generateK(X=X,
                        b=b,
                        control=ctrl)
    }
  }
  
  pars <- rep(0, ncol(Kdat$U) + 1)

  loss <- 0
  for(j in 1:folds){
    fold <- chunks[[j]]
    UFold <- Kdat$U[-fold, ]

    suppressWarnings({suppressMessages({invisible(capture.output({
    out <- getDhat(par = pars,
                   y=y[-fold],
                   U=UFold,
                   D=Kdat$D,
                   w=ctrl$w[-fold],
                   lambda=lambda,
                   ctrl = ctrl)
    }))})})
    if(ctrl$loss == "leastsquares") {
      yhat <- Kdat$U[fold, ] %*% out$dhat
      loss <- loss + sum((y[fold] - yhat)^2)

    } else if (ctrl$loss == "logistic") {
      yhat <- logistic(Kdat$U[fold, ], out$dhat, out$beta0hat)
      yhat[yhat==1]=1-1e-12
      yhat[yhat==0]=0+1e-12
      loss <- loss - sum(y[fold]*log(yhat)+(1-y[fold])*log(1-yhat))
    }

  }
  
  if(ctrl$loss == "leastsquares") {
    loss <- sqrt(loss/length(y))
  } else if(ctrl$loss == 'logistic') {
    loss <- loss/length(y)
  }

  print(loss)
  return(loss)
}

## Lambda search for KRLS
#' @export
lambdaline <-
  function(Lbound=NULL,
           Ubound=NULL,
           y=NULL,
           U=NULL,
           D=NULL,
           w=NULL,
           tol=NULL,
           noisy=FALSE){

    n <- length(y)
    if(is.null(tol)){
      tol <- 10^-3 * n
    } else {
      stopifnot(is.vector(tol),
                length(tol)==1,
                is.numeric(tol),
                tol>0)
    }

    # get upper bound starting value
    if(is.null(Ubound)){
      Ubound <- n

      while(sum(D / (D + Ubound)) < 1){
        Ubound <- Ubound-1
      }
    } else {
      stopifnot(is.vector(Ubound),
                length(Ubound)==1,
                is.numeric(Ubound),
                Ubound>0)
    }

    # get lower bound starting value
    if(is.null(Lbound)){
      q <- which.min(abs(D - (max(D)/1000)))

      #Lbound <- 0
      Lbound = .Machine$double.eps  #CJH: to avoid Inf in next statement

      #while(sum(D / (D + Lbound)) > q){
      #  Lbound <- Lbound+.05
      #}
    }  else {
      stopifnot(is.vector(Lbound),
                length(Lbound)==1,
                is.numeric(Lbound),
                Lbound>=0)
    }
    # create new search values
    X1 <- Lbound + (.381966)*(Ubound-Lbound)
    X2 <- Ubound - (.381966)*(Ubound-Lbound)

    # starting LOO losses
    S1 <- looloss(lambda=X1,y=scale(y),U=U, D=D,w=w)
    S2 <- looloss(lambda=X2,y=scale(y),U=U, D=D,w=w)

    if(noisy){cat("Lbound:",Lbound,"X1:",X1,"X2:",X2,"Ubound:",Ubound,"S1:",S1,"S2:",S2,"\n") }

    while(abs(S1-S2)>tol){ # terminate if difference between S1 and S2 less than tolerance

      # update steps and use caching
      if(S1 < S2){
        Ubound  <- X2
        X2 <- X1
        X1 <- Lbound + (.381966)*(Ubound-Lbound)
        S2 <- S1
        S1 <- looloss(lambda=X1,y=scale(y),U=U, D=D,w=w)

      } else { #S2 < S1
        Lbound  <- X1
        X1 <- X2
        X2 <- Ubound - (.381966)*(Ubound-Lbound)
        S1 <- S2
        S2 <- looloss(lambda=X2,y=scale(y),U=U, D=D,w=w)
      }

      if(noisy){cat("Lbound:",Lbound,"X1:",X1,"X2:",X2,"Ubound:",Ubound,"S1:",S1,"S2:",S2,"\n") }
    }
    out <- ifelse(S1<S2,X1,X2)
    if(noisy){cat("Lambda:",out,"\n")}
    return(invisible(out))
  }


## looloss for krls
#' @export
looloss <-
  function(y=NULL,K=NULL,D=NULL,U=NULL,w=NULL,lambda=NULL){
    print(D[1]);print(y[1])
      if(is.null(w)) {
        Le <- solve_for_d_ls(y=y,D=D,U=U,lambda=lambda)$Le
      } else {
        Le <- solve_for_d_ls_w(y=y,D=D,U=U,w=w,lambda=lambda)$Le
      }
  }
