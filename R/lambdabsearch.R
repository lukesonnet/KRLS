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
                              hyperctrl,
                              control) {
  
  if(is.null(hyperctrl$lambdarange)) {
    if (!is.null(hyperctrl$brange)) stop("Grid search for b only works with fixed lambda or lambda grid search.")
    
    fit.hyper <- optim(par=log(c(hyperctrl$lambdastart, 2.2*control$d)), lambdab.fn,
                       X=X, y=y, folds=length(hyperctrl$chunks),
                       chunks = hyperctrl$chunks, ctrl = control,
                       vcov = FALSE,
                       control=list(trace = T, abstol = 1e-4), method="BFGS")
    
    lambda <- exp(fit.hyper$par[1])
    b <- exp(fit.hyper$par[2])
  } else {
    if (hyperctrl$optimb) stop("Optimizing b does not work with a grid search for lambda.")
    hypergrid <- as.matrix(expand.grid(hyperctrl$lambdarange, hyperctrl$brange))
    hyperMSE <- NULL
    for(i in 1:nrow(hypergrid)){
      hyperMSE[i] <- lambdab.fn(par = log(hypergrid[i, ]), X=X,
                                    y=y, folds=length(hyperctrl$chunks),
                                    chunks = hyperctrl$chunks, ctrl = control)
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
                         hyperctrl,
                         control,
                         b) {
  
  if (control$loss == "leastsquares") {
    if(control$truncate) {
      #todo: no way this lambdasearch is right. The Looe can't be right, the bounds can't be right... but it works
      #stop("Must specify lambda for truncated least squares for now.")
      lambda <- lambdaline(y=y, D=Kdat$eigvals, Utrunc=Kdat$Utrunc,
                             Eigenobject=Kdat$eigobj, truncate=control$truncate, noisy = !control$quiet)#,eigtrunc=eigtrunc,noisy=noisy,L=L,U=U)
    } else {
      lambda <- lambdaline(y=y, K=Kdat$K, Eigenobject=Kdat$eigobj,
                             truncate=control$truncate, noisy = !control$quiet)#,eigtrunc=eigtrunc,noisy=noisy,L=L,U=U)
    }
  } else {
    
    if(is.null(hyperctrl$lambdarange)) {
      
      fit.lambda <- optim(par=log(hyperctrl$lambdastart), lambdab.fn,
                          Kdat = Kdat, y=y, folds=length(hyperctrl$chunks),
                          chunks = hyperctrl$chunks, ctrl = control,
                          b = b,
                          vcov = FALSE,
                          control=list(trace = T, abstol = 1e-4), method="BFGS")
      
      lambda <- exp(fit.lambda$par)
      
    } else {
      
      lambdaMSE <- NULL
      for(i in 1:length(hyperctrl$lambdarange)){
        lambdaMSE[i] <- lambdab.fn(par = log(hyperctrl$lambdarange[i]), Kdat = Kdat,
                                       y=y, folds=length(hyperctrl$chunks),
                                       chunks = hyperctrl$chunks, ctrl = control,
                                       b = b)
      }
      lambda <- hyperctrl$lambdarange[which.min(lambdaMSE)]
    }
  }
  
  return(lambda)
  
}

## This function searches for just b
#' @export
bsearch <- function(y,
                              X,
                              hyperctrl,
                              control,
                              lambda) {
  if(is.null(hyperctrl$brange)) {
    fit.b <- optim(par=log(2.2*control$d), lambdab.fn,
                       X=X, y=y, folds=length(hyperctrl$chunks),
                       chunks = hyperctrl$chunks, ctrl = control, 
                       lambda = lambda,
                       vcov = FALSE,
                       control=list(trace = T, abstol = 1e-4), method="BFGS")
    b <- exp(fit.b$par)
  } else {
    bMSE <- NULL
    for(i in 1:length(hyperctrl$brange)){
      if(hyperctrl$brange[i] < 0) stop("All bs must be positive")
      bMSE[i] <- lambdab.fn(par = log(hyperctrl$brange[i]), X=X,
                                    y=y, folds=length(hyperctrl$chunks),
                                    chunks = hyperctrl$chunks, ctrl = control,
                                    #b = b,
                                    lambda = lambda)
      b <- hyperctrl$brange[which.min(bMSE)]
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
      lambda <- exp(par)
    }
    if (is.null(b)) {
      b <- exp(par)
      Kdat <- generateK(X=X,
                        b=b,
                        control=ctrl)
    }
  }
  
  pars <- rep(0, ifelse(ctrl$truncate, ncol(Kdat$Utrunc), ncol(Kdat$K)) + 1)
  
  mse <- 0
  for(j in 1:folds){
    fold <- chunks[[j]]
    if(ctrl$truncate) {
      #cjh added lastkeeper to this, to reuse it rather than re-search
      KFold <- NULL
      UtruncFold <- Kdat$Utrunc[-fold, ]
      D <- Kdat$eigvals
    } else {
      KFold <- Kdat$K[-fold, -fold]
      UtruncFold <- NULL
      D <- NULL
    }
    
    if(ctrl$loss == "leastsquares") {
      if(ctrl$truncate) {
        parshat <- solve_for_c_ls_trunc(y = y[-fold], Utrunc = UtruncFold, D = D, lambda = lambda)
        yhat <- Kdat$Utrunc[fold, ] %*% parshat$chat
      } else {
        print(dim(KFold))
        print(dim(Kdat$K))
        print(length(y[-fold]))
        parshat <- solve_for_c_ls(y = y[-fold], K = KFold, lambda = lambda)
        yhat <- Kdat$K[fold, -fold] %*% parshat$chat
      }

    } else if (ctrl$loss == "logistic") {
      parshat <- solveForC(par = pars, y=y[-fold], K=KFold, Utrunc=UtruncFold, D=D, lambda=lambda)
      
      if(ctrl$truncate) {
        yhat <- logistic(Kdat$Utrunc[fold, ], parshat$chat, parshat$beta0hat)
      } else {
        yhat <- logistic(Kdat$K[fold, -fold], parshat$chat, parshat$beta0hat)
      }
    }
    
    #K <- newKernel(X[-fold, ], newData = X[fold, ])
    ## Is this transformation right?

    if(ctrl$loss == "logistic") {

    } else {

    }
    
    mse <- mse + sum((y[fold] - yhat)^2)
  }
  rmspe <- sqrt(mse/length(y))
  
  return(rmspe)
}

## Lambda search for KRLS
#' @export
lambdaline <-
  function(L=NULL,
           U=NULL,
           y=NULL,
           K=NULL,
           D=NULL,
           Utrunc=NULL,
           Eigenobject=NULL,
           tol=NULL,
           noisy=FALSE,
           eigtrunc=NULL,
           truncate=NULL){
    
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
    if(is.null(U)){
      U <- n
      
      while(sum(Eigenobject$values / (Eigenobject$values + U)) < 1){
        U <- U-1    
      }
    } else {
      stopifnot(is.vector(U),
                length(U)==1,
                is.numeric(U),
                U>0)
    }
    
    # get lower bound starting value
    if(is.null(L)){
      q <- which.min(abs(Eigenobject$values - (max(Eigenobject$values)/1000)))    
      
      #L <- 0
      L = .Machine$double.eps  #CJH: to avoid Inf in next statement
      
      while(sum(Eigenobject$values / (Eigenobject$values + L)) > q){
        L <- L+.05    
      }
    }  else {
      stopifnot(is.vector(L),
                length(L)==1,
                is.numeric(L),
                L>=0)
    }
    # create new search values    
    X1 <- L + (.381966)*(U-L)
    X2 <- U - (.381966)*(U-L)
    
    # starting LOO losses
    S1 <- looloss(lambda=X1,y=y,K=K, D=D, Utrunc=Utrunc,eigtrunc=eigtrunc,truncate=truncate)
    S2 <- looloss(lambda=X2,y=y,K=K, D=D, Utrunc=Utrunc,eigtrunc=eigtrunc,truncate=truncate)
    
    if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }
    
    while(abs(S1-S2)>tol){ # terminate if difference between S1 and S2 less than tolerance
      
      # update steps and use caching
      if(S1 < S2){
        U  <- X2
        X2 <- X1
        X1 <- L + (.381966)*(U-L)
        S2 <- S1
        S1 <- looloss(lambda=X1,y=y,K=K, D=D, Utrunc=Utrunc,eigtrunc=eigtrunc,truncate=truncate)
        
      } else { #S2 < S1
        L  <- X1
        X1 <- X2
        X2 <- U - (.381966)*(U-L)
        S1 <- S2
        S2 <- looloss(lambda=X2,y=y,K=K, D=D, Utrunc=Utrunc,eigtrunc=eigtrunc,truncate=truncate)
      }
      
      if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }
    }
    out <- ifelse(S1<S2,X1,X2)
    if(noisy){cat("Lambda:",out,"\n")}  
    return(invisible(out))
  }


## looloss for krls
#' @export
looloss <-
  function(y=NULL,K=NULL,D=NULL,Utrunc=NULL,lambda=NULL,eigtrunc=NULL,truncate=NULL){
    if (truncate) {
      Le <- solve_for_c_ls_trunc(y=y,D=D,Utrunc=Utrunc,lambda=lambda)$Le
      return(Le)
    } else {
      Le <- solve_for_c_ls(y=y,K=K,lambda=lambda)$Le
      return(Le)
    }
  }