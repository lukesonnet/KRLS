#----------
# Roxygen Commands
#----------

#' @import methods
#' @import rARPACK
NULL
#' @useDynLib KRLS2
#' @importFrom Rcpp sourceCpp
NULL

#############
# Functions #
#############

## This is the main function. Fits the model by preparing the data, searching
## for lambda through CV, minimizing the target function, and returning an
## object that contains all of the information about the fitting that is
## necessary for future analysis. I borrow liberally from the KRLS package

#' @export
krls <- function(# Data arguments
                    X = NULL,
                    y = NULL,
                    # Family
                    loss = "leastsquares",
                    # Kernel arguments
                    whichkernel = "gaussian", 
                    sigma = NULL,
                    optimsigma = FALSE,
                    # Lambda arguments
                    lambda = NULL,
                    hyperfolds = 5,
                    lambdastart = 0.5,
                    L = NULL,
                    U = NULL,
                    # Truncation arguments
                    truncate = FALSE,
                    lastkeeper = NULL,
                    epsilon = 0.01,
                    # Optimization arguments
                    con = list(maxit=500),
                    printout = TRUE,
                    returnopt = TRUE,
                    # Standard error arguments
                    vcov = TRUE,
                    derivative = TRUE,
                    cpp = TRUE,
                    clusters = NULL) {
  
  ## Prepare the data
  X <- as.matrix(X)
  y <- as.matrix(y)
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (is.numeric(X)==FALSE) {
    stop("X must be numeric")
  }      
  if (is.numeric(y)==FALSE) {
    stop("y must be numeric")
  }
  if (sd(y)==0) {
    stop("y is a constant")        
  }  
  if (sum(is.na(X))>0){
    stop("X contains missing data")
  }
  if (sum(is.na(y))>0){
    stop("y contains missing data")
  }
  if (n!=nrow(y)){
    stop("nrow(X) not equal to number of elements in y")
  }
  
  #if (truncate & loss == "leastsquares")
  #  stop("Truncation is not allowed with least squares for now")
  
  stopifnot(
    is.logical(derivative),
    is.logical(vcov)#,
    #is.logical(binary)
  )
  
  if(derivative==TRUE){
    if(vcov==FALSE){
      stop("derivative==TRUE requires vcov=TRUE")
    }
    if(truncate==FALSE & loss == "logistic") {
      stop("derivative==TRUE requires truncate==TRUE if loss == 'logistic'")
    }
  }
  
  # column names
  if (is.null(colnames(X))) {
    colnames(X) <- paste("x",1:d,sep="")
  }
  
  ## Scale data
  X.init <- X
  X.init.sd <- apply(X.init, 2, sd)
  if (sum(X.init.sd == 0)){
    stop("at least one column in X is a constant, please remove the constant(s)")
  }
  X <- scale(X, center = TRUE, scale = X.init.sd)    
  
  if (loss == "leastsquares") {
    ## Scale y
    y.init <- y
    y.init.sd <- apply(y.init,2,sd)
    y.init.mean <- mean(y.init)
    y <- scale(y,center=y.init.mean,scale=y.init.sd)
    
  } else if (loss == "logistic") {
    if (!all(y %in% c(0,1))) {
      stop("y is not binary data")
    }
  }
  
  lambdarange <- NULL
  sigmarange <- NULL
  ## Default sigma to the number of features
  if(is.null(sigma)){
    if (!optimsigma) {sigma <- d}
  }  else{
    if(length(sigma) > 1) {
      sigmarange <- sigma
      sigma <- NULL
    } else {
      stopifnot(is.vector(sigma),
                is.numeric(sigma),
                sigma>0)
    }
  }
  
  
  ## todo: unpack truncdat and add checks for if truncate is true, do not have
  ## two seperate processes but instead checks at pertinent steps of process and
  ## try to use similar variable names to economize on code
  
  ## todo: build distance matrix E first so that we can more quickly compute
  ## K below
  
  if(!is.null(lambda)) {
    # Allow range in lambda argument
    if (length(lambda) > 1) {
      lambdarange <- lambda
      lambda <- NULL
    } else {
      # check user specified lambda
      stopifnot(is.vector(lambda),
                is.numeric(lambda),
                lambda>0) 
    }
  }
  
  ## set chunks for any ranges
  if(!is.null(c(lambdarange, sigmarange))) {
    chunks <- chunk(sample(n), hyperfolds)
  }
  
  ## lastkeeper is recomputed for now, might want to use a guess instead
  if(!is.null(lambda)) {
    
    if(is.null(sigma)) {
      if(is.null(sigmarange)) {
        fit.sigma <- optim(par=lambdastart, lambdasigma.fn,
                           X=X.init, y=y, folds=hyperfolds, chunks = chunks,
                           truncate=truncate, epsilon=epsilon, lastkeeper=lastkeeper,
                           #sigma = sigma,
                           lambda = lambda,
                           vcov = FALSE,
                           control=list(trace = T, abstol = 1e-5), method="BFGS")
        sigma <- exp(fit.sigma$par)
        message(sprintf("Sigma selected: %s", sigma))
      } else {
        for(i in 1:length(sigmarange)){
          if(sigmarange[i] < 0) stop("All sigmas must be positive")
          sigmaMSE[i] <- lambdasigma.fn(par = log(sigmarange[i]), X=X.init,
                                    y=y, folds=hyperfolds, chunks=chunks,
                                    truncate=truncate, epsilon=epsilon,
                                    lastkeeper = lastkeeper,
                                    #sigma = sigma,
                                    lambda = lambda)
        }
        sigma <- sigmarange[which.min(sigmaMSE)]
        message(sprintf("Sigma selected: %s", sigma))
      }

      
    }
  } else {
    
    if (loss == "leastsquares") {
      if(is.null(sigma)) stop("Cannot simultaneously search for both lambda and sigma with leastsquares loss")
      lambda <- lambdasearch(#L=L,U=U,
        y=y,K=K,Eigenobject=eigobj,truncate=truncate)#,eigtrunc=eigtrunc,noisy=noisy)
    } else {
      
      if(is.null(lambdarange)) {
        if (!is.null(sigmarange)) stop("Grid search for sigma only works with fixed lambda or lambda grid search.")
        if (is.null(sigma)) {
          fit.hyper <- optim(par=c(lambdastart, 2.2*d), lambdasigma.fn,
                             X=X.init, y=y, folds=hyperfolds, chunks = chunks,
                             truncate=truncate, epsilon=epsilon, lastkeeper=lastkeeper,
                             #sigma = sigma,
                             #lambda = lambda,
                             vcov = FALSE,
                             control=list(trace = T, abstol = 1e-5), method="BFGS")
          
          lambda <- exp(fit.hyper$par[1])
          sigma <- exp(fit.hyper$par[2])
          
          message(sprintf("Lambda selected: %s", lambda))
          message(sprintf("Sigma selected: %s", sigma))
        } else {
          fit.lambda <- optim(par=log(lambdastart), lambdasigma.fn,
                             X=X.init, y=y, folds=hyperfolds, chunks = chunks,
                             truncate=truncate, epsilon=epsilon, lastkeeper=lastkeeper,
                             sigma = sigma,
                             #lambda = lambda,
                             vcov = FALSE,
                             control=list(trace = T, abstol = 1e-5), method="BFGS")
          
          lambda <- exp(fit.lambda$par)

          message(sprintf("Lambda selected: %s", lambda))
        }

      } else {
        if (is.null(sigma)) stop("Grid search for lambda does not work with optimizing sigma. Set optimsigma = F or use a grid search for sigma.")
        
        lambdaMSE <- NULL
        for(i in 1:length(lambdarange)){
          lambdaMSE[i] <- lambdasigma.fn(par = log(lambdarange[i]), X=X.init,
                                         y=y, folds=hyperfolds, chunks=chunks,
                                         truncate=truncate, epsilon=epsilon,
                                         lastkeeper = lastkeeper,
                                         sigma = sigma)
                                         #lambda = lambda)
        }
        lambda <- lambdarange[which.min(lambdaMSE)]
        message(sprintf("Lambda selected: %s", lambda))
      }
    }
  }
 
  ## Compute kernel matrix
  K <- NULL
  if(whichkernel=="gaussian"){ K <- kern_gauss(X, sigma)}
  if(whichkernel=="linear"){ K <- tcrossprod(X) }
  if(whichkernel=="poly2"){ K <- (tcrossprod(X)+1)^2 }
  if(whichkernel=="poly3"){ K <- (tcrossprod(X)+1)^3 }
  if(whichkernel=="poly4"){ K <- (tcrossprod(X)+1)^4 }
  if(is.null(K)){ stop("No valid Kernel specified") }
  
  Utrunc <- NULL
  if(truncate) {
    full <- FALSE
    if (loss == "leastsquares") {
      full <- TRUE
    }
    truncDat <- Ktrunc(K = K, sigma=sigma, lastkeeper=lastkeeper, epsilon=epsilon,
                       full = full)
    Utrunc <- truncDat$Utrunc
    eigvals <- truncDat$eigvals
    if (full) {
      eigobj <- truncDat$eigobj
    }
  } else {
    if (loss == "leastsquares") {
      eigobj <- eigen(K)
    }
    Utrunc <- NULL
    eigvals <- NULL
  }
  
  #rm(Kfull) ## remove Kfull because it is probably quite a large object
  #CJH: we'll need to keep Kfull for derivative
  
  lastkeeper <- NULL
  lastkeeper <- if(truncate) ncol(Utrunc)
  
  
  ## Solve
  vcovmatc=NULL
  
  if (loss == "leastsquares") {
    if (truncate){
      out <- solve_for_c_ls_trunc(y=y, K=K,Utrunc =Utrunc, lambda=lambda)
      print(Utrunc%*%out$coeffs)
    } else {
      out <- solve_for_c_ls(y=y, K=K, lambda=lambda)
      chat <- out$coeffs
    }
    
    ## getting c_p from d
    coefhat=NULL
    if(truncate==FALSE){coefhat=chat} else {
      UDinv = mult_diag(Utrunc, 1/eigvals)
      coefhat=UDinv %*% chat
    }
    
    yfitted <- K%*%coefhat
    
    if (vcov) {
      sigmasq <- as.vector((1/n) * crossprod(y-yfitted))
      
      vcovmatc <- tcrossprod(mult_diag(eigobj$vectors,sigmasq*(eigobj$values+lambda)^-2),eigobj$vectors)   
      
    }
    beta0hat <- NULL
  } else {
    out <- solveForC(y=y, K=K, Utrunc=Utrunc, D=eigvals, lambda=lambda, con=con, vcov=vcov, printout = printout, returnopt = returnopt)
    chat <- out$chat
    beta0hat <- out$beta0hat
    hessian <- out$hessian
    opt <- if(returnopt) out$opt else NULL
    
    ## getting c_p from d
    coefhat=NULL
    if(truncate==FALSE){coefhat=chat} else {
      UDinv = mult_diag(Utrunc, 1/eigvals)
      coefhat=UDinv %*% chat
    }
    
    yfitted <- logistic(K=K, coeff=coefhat, beta0 = beta0hat)
    
    ## Vcov of d
    if (vcov & truncate) {
      vcov.cb0 = solve(hessian)
      vcovmatc = tcrossprod(UDinv%*%vcov.cb0[1:length(chat), 1:length(chat)], UDinv)
      ## todo: return var b0
      nclust <- length(clusters)
      score <- matrix(NA, nclust, length(c(chat, beta0hat)))
      for (j in 1:nclust) {
        score[j, ] = krlogit_gr_trunc2(par=c(chat,beta0hat), Utrunc[clusters[[j]], , drop = F], eigvals, y[clusters[[j]], drop = F], lambda, n/length(clusters[[j]]))
      }
      
      print(nclust)
      #score = krlogit_gr_trunc(par=c(chat,beta0hat), Utrunc, eigvals, y, lambda)
      meat <- (nclust/(nclust-1)) * crossprod(score)
      vcov.rob.cb0 = vcov.cb0 %*% meat %*% vcov.cb0
      vcovmatc.rob = tcrossprod(UDinv%*%vcov.rob.cb0[1:length(chat), 1:length(chat)], UDinv)
      
      #vcovrobmatc = 
    } else {vcov.cb0 = NULL}
    
  }
  
  ###----------------------------------------
  ## Getting pwmfx
  ###----------------------------------------
  if (derivative==TRUE){
  derivmat <- matrix(NA, ncol=d, nrow=n,
                     dimnames=list(NULL, colnames(X)))
  varavgderivmat <- matrix(NA,1,d)
  
  if(loss == "leastsquares") {
    p1p0 <- rep(2, n)
  } else if (loss == "logistic") {
    p1p0 <- yfitted*(1-yfitted)
  }
  
  #construct coefhat=c for no truncation and coefhat = Utrunc*c
  #to take the truncation into the coefficients. 

  if(cpp) {

    derivout <- pwmfx(K, X, coefhat, vcovmatc, p1p0, sigma)
    derivmat <- derivout[1:n, ]
    varavgderivmat <- derivout[n+1, ]

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
      L <-  distk*Kfull  #Kfull rather than K here as truncation handled through coefs
      derivmat[,k] <- p1p0*(-1/sigma)*(L%*%coefhat)
      if(truncate==FALSE) {
        varavgderivmat = NULL
      } else {
        varavgderivmat[1,k] <- 1/(sigma * n)^2 * sum(crossprod(p1p0^2, crossprod(L,vcovmatc%*%L)))
      }
    }
    
  }
    #donebelow
    #avgderiv <- matrix(colMeans(derivmat),nrow=1)
    #colnames(avgderiv) <- colnames(X)
    
    
    #print(colMeans(X.init))
    #print(X.init.sd)
    #print(derivmat)
    #effectmat <- t((1/X.init.sd)*t(derivmat))
    # or, equivalently, and probably more efficiently
    ## todo: smart checker if Y was rescaled?
    if (loss == "leastsquares") {
      avgderivmat <- colMeans(derivmat)
      derivmat <- scale(y.init.sd*derivmat,center=FALSE,scale=X.init.sd)
      attr(derivmat,"scaled:scale")<- NULL
      avgderiv <- scale(y.init.sd*matrix(avgderivmat, nrow = 1),center=FALSE,scale=X.init.sd)
      attr(avgderiv,"scaled:scale")<- NULL
      varavgderivmat <- (y.init.sd/X.init.sd)^2*varavgderivmat
      attr(varavgderivmat,"scaled:scale")<- NULL
      yfitted     <- yfitted*y.init.sd+y.init.mean
    } else {
      derivmat <- scale(derivmat, center=F, scale = X.init.sd)
      avgderiv <- colMeans(derivmat)
      varavgderivmat <- (1/X.init.sd)^2*varavgderivmat
    }
  

  } else {
    derivmat=NULL
    avgderiv=NULL
    varavgderivmat=NULL
  }
  
  if (truncate) {
    scorei <- matrix(NA, n, length(c(chat, beta0hat)))
    for(i in 1:n) {
      scorei[i, ] = krlogit_gr_trunc2(par=c(chat,beta0hat), Utrunc[i, , drop = F], eigvals, y[i, drop = F], lambda, n)
    }
    score2 = krlogit_gr_trunc(par=c(chat,beta0hat), Utrunc, eigvals, y, lambda)
  } 
  #else {score = krlogit.gr(par=c(chat,beta0hat), K=K, y=y,lambda=lambda)}
  
  # return
  z <- list(K=K,
            Utrunc=Utrunc,
            lastkeeper = lastkeeper,
            truncate = truncate,
            coeffs=coefhat,
            chat=chat,
            beta0hat = beta0hat,
            fitted=yfitted,
            X=X.init,
            y=y,
            sigma=sigma,
            lambda=lambda,
            kernel=whichkernel,
            derivmat = derivmat,
            avgderiv = avgderiv,
            var.avgderiv = varavgderivmat,
            #lastkeeper = ncol(Ktilde)  
            score=score,
            scorei=scorei,
            score2=score2,
            vcovmatc.rob = vcovmatc.rob,
            vcovmatc = vcovmatc,
            vcov.cb0 = vcov.cb0,
            vcov.rob.cb0 = vcov.rob.cb0
            #vcov.cb0 = vcov.cb0,
            #opt = opt
  )
    
  class(z) <- "krlogit"  
  
  return(z)
}

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

#krlogit.hess <- function(par, K, y, lambda){
#  cb <- crossprod(K, (exp(-K%*%chat - beta0)/(1 + exp(-K%*%chat - beta0))^2))
#  bb <- sum(exp(-K%*%chat - beta0)/(1 + exp(-K%*%chat - beta0))^2)
#  cc <- can't get this to match what optim returns, see krlogit_testingmath
  

# todo: Unless the folds are much smaller, lastkeeper will probably be the same
# or only one or two smaller, should we pass lastkeeper and just use it?
## This computes the RMSPE using 'folds' cross-validation, fitting a new model
## for each subset of the data
#' @export
lambdasigma.fn <- function(par = NULL,
                      X = NULL,
                      y = NULL,
                      folds = NULL,
                      chunks = NULL,
                      sigma = NULL,
                      lambda = NULL,
                      truncate = NULL,
                      vcov = NULL,
                      epsilon = NULL,
                      lastkeeper = NULL) {
  
  if(is.null(c(lambda, sigma))) {
    lambda <- exp(par[1])
    sigma <- exp(par[2])
  } else {
    if (is.null(lambda)) {
      lambda <- exp(par)
    }
    if (is.null(sigma)) {
      sigma <- exp(par)
    }
  }
  
  mse <- 0
  for(j in 1:folds){
    fold <- chunks[[j]]
    if(truncate) {
      #cjh added lastkeeper to this, to reuse it rather than re-search
      truncDat <- Ktrunc(X = X[-fold,], sigma=sigma, epsilon=epsilon,lastkeeper=lastkeeper) 
      K <- NULL
      Utrunc <- truncDat$Utrunc
      eigvals <- truncDat$eigvals
    } else {
      K <- kern_gauss(scale(X[-fold,]), sigma)
      Utrunc <- NULL
      eigvals <- NULL
    }
    pars <- rep(0, ifelse(truncate, ncol(Utrunc), ncol(K)) + 1)
    parshat <- solveForC(par = pars, y=y[-fold], K=K, Utrunc=Utrunc, D=eigvals, lambda=lambda, vcov = FALSE)

    K <- newKernel(X[-fold, ], newData = X[fold, ])
    ## Is this transformation right?
    if(truncate) {
      coefhat <- mult_diag(Utrunc, 1/eigvals) %*% parshat$chat
    } else {
      coefhat <- parshat$chat
    }
    yhat <- logistic(K=K, coefhat, beta0=parshat$beta0hat)
    mse <- mse + sum((y[fold] - yhat)^2)
  }
  rmspe <- sqrt(mse/length(y))
    
  return(rmspe)
}

## Function that returns truncated versions of the data if given
## todo: throws warning whenever n < 500 because it uses eigen, should we suppress?
#' @export
Ktrunc <- function(X=NULL, K=NULL, sigma=NULL, epsilon=NULL, lastkeeper=NULL, full=FALSE){
  if(is.null(K)){
    X.sd <- apply(X, 2, sd)
    X.mu <- colMeans(X)
    K <- kern_gauss(scale(X), ifelse(is.null(sigma), ncol(X), sigma))
  }
  
  ## Todo: maybe later allow for full return of eigs_sym but not full
  ## eigen decomposition. Might be useful for lambda search for logistic
  if(full){
    eigobj <- eigen(K)
    
    lastkeeper <- ifelse(is.null(lastkeeper),
                        min(which(cumsum(eigobj$values)/nrow(K) > (1-epsilon))),
                        lastkeeper)
    
    Ktilde <- mult_diag(eigobj$vectors[, 1:lastkeeper], eigobj$values)
    
    print(paste("lastkeeper=",lastkeeper))
    
    return(
      list(Ktilde=Ktilde,
           Utrunc=eigobj$vectors[, 1:lastkeeper],
           eigobj=eigobj)
    )
    
  }
  if(is.null(lastkeeper)){
    #denoms <- c(10, 6, 3, 1) # we could also let people set this to speed up the algorithm
    #CJH: even at 10, with N=5000 it is too big. Let's try this: 
    if (nrow(K) <= 500) numvectorss=nrow(K) else numvectorss=c(50, 100, 200, 300,500,1000)
  
    enoughvar=FALSE
    j=1
    while (enoughvar==FALSE){
      numvectors=numvectorss[j]
      print(paste("trying",numvectors,"vectors"))
      eigobj <- NULL
      eigobj <- try({eigs_sym(K, numvectors, which="LM")}, silent = T)
    
      #for now, letting it throw an error in certain failure cases.
      totalvar=sum(eigobj$values)/nrow(K)
    
      #notice I'm assuming epsilon is the remaining bit, e.g. .01 not .99.
      #changed default to match.
      if (totalvar>=(1-epsilon)){
        lastkeeper = min(which(cumsum(eigobj$values)/nrow(K) > (1-epsilon)))
        Ktilde <- mult_diag(eigobj$vectors[, 1:lastkeeper], eigobj$values)
        
        print(paste("lastkeeper=",lastkeeper))
        enoughvar=TRUE
        
      }
      j=j+1
      # Will eventually add warning here for var < 0.99 at numvector= 300
      # Allow to proceed w/out truncation and SEs or to try full eigen
    } #end while gotenough==FALSE      
  } else { #end if is.null(lastkeeper)
    eigobj <- eigs_sym(K, lastkeeper, which="LM")
    Ktilde <- mult_diag(eigobj$vectors, eigobj$values)
  }
  
  return(
         list(Utrunc=eigobj$vectors[, 1:lastkeeper],
              eigvals=eigobj$values[1:lastkeeper])
        )
}


## Computes the Gaussian kernel matrix from the data matrix X
## Parameters:
##   'X' - the data matrix
##   'sigma' - the kernel bandwitch, recommended by HH2013 to be ncol(X)
## Values:
##   'K' - the Gaussian kernel matrix
#' @export
gaussKernel <- function(X=NULL, sigma=NULL) {
  return( exp(-1*as.matrix(dist(X)^2)/(2*sigma)) )
}


## Optimize C in 'krlogit.fn' given the data and lambda
## Parameters:
##   'par' - starting parameters
##   'y' - the outcome variable
##   'K' - the kernel Matrix
##   'lambda' - the regularizing parameter
## Values:
##   a list containing three objects:
##     'chat' - the estimated values for c, a vector of length n
##     'fullopt' - the full object that optim returns, for diagnostic purposes
#' @export
solveForC <- function(par= NULL,
                      y = NULL,
                      K = NULL,
                      Utrunc = NULL,
                      D = NULL,
                      lambda = NULL,
                      con = list(),
                      vcov = NULL,
                      printout = FALSE,
                      returnopt = TRUE) {
  
  ncoeffs <- ifelse(is.null(D), ncol(K), ncol(Utrunc))
  par <- rep(0, ncoeffs + 1)
  
  if(is.null(D)){
    opt <- optim(par, krlogit.fn, gr=krlogit.gr, K=K, y=y, lambda=lambda,
                 method="BFGS", control=con, hessian = vcov)
  } else {
    opt <- optim(par, krlogit_fn_trunc, gr=krlogit_gr_trunc, Utrunc=Utrunc, D = D, y=y, lambda=lambda,
                 method="BFGS", control=con, hessian = vcov)
  }
  
  chat <- opt$par[1:ncoeffs]
  beta0hat <- opt$par[ncoeffs+1]
  
  if (!returnopt) opt <- NULL
  if (printout) {
    print(paste("Calls to function:", opt$counts[1], ", calls to gradient:", opt$counts[2]))
    print(paste("Converged? ", as.logical(1 - opt$convergence)))
  }
  
  return(list(chat = chat,
              beta0hat = beta0hat,
              hessian = opt$hessian,
              opt = opt))
}


## a predict function for class 'krlogit'
#' @export
predict.krlogit <- function(object, newdata, ...) {
  if (class(object) != "krlogit") {
    warning("Object not of class 'krlogit'")
    UseMethod("predict")
    return(invisible(NULL))
  }

  newdata <- as.matrix(newdata)
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted krls object")
  }
  
  newdataK <- newKernel(X = object$X, newData = newdata, whichkernel = object$kernel)
  
  #if(object$truncate){
  #  newdataK <- newdataK%*%object$Utrunc
  #}
  
  yfitted <- logistic(K = newdataK, coeff = object$coeffs, beta0 = object$beta0hat)
  return(list(fit = yfitted,
              #se.fit = se.fit, vcov.fit = vcov.fit, newdata = newdata, 
              newdataK = newdataK))
    #if (se.fit) {
    #  vcov.c.raw <- object$vcov.c * as.vector((1/var(object$y)))
    #  vcov.fitted <- tcrossprod(newdataK %*% vcov.c.raw, newdataK)
    #  vcov.fit <- (apply(object$y, 2, sd)^2) * vcov.fitted
    #  se.fit <- matrix(sqrt(diag(vcov.fit)), ncol = 1)
    #}
    #else {
    #  vcov.fit <- se.fit <- NULL
    #}
    #yfitted <- (yfitted * apply(object$y, 2, sd)) + mean(object$y)
  
}

## The logistic function that takes values for coeff, b0, and a K or Ktilde
#' @export
logistic <- function(K, coeff, beta0) {
  
  yhat <- 1 / (1 + exp(-(beta0+K%*%coeff)))
  
  return(yhat)
}

## Produce a new Kernel from new data and data that was used to fit model for
## prediction
## Parameters:
##   'X' - the original data matrix, UNSCALED
##   'newData' - the new data matrix (with same features), unscaled
##   'whichkernel' - which kernel was used on the original data
## Values:
##   'newK' - The new Kernel to be used for prediction
#' @export
newKernel <- function(X, newData, whichkernel = "gaussian") {
  
  # Get values of oldData to scale newData
  Xmeans <- colMeans(X)
  Xsd <- apply(X, 2, sd)
  X <- scale(X, center = Xmeans, scale = Xsd)

  # scale new data by means and sd of old data
  newData <- scale(newData, center = Xmeans, scale = Xsd)      
  
  # predict based on new kernel matrix
  nn <- nrow(newData)
  
  ## Compute kernel matrix
  K <- NULL
  if(whichkernel=="gaussian"){ 
    newK <- new_gauss_kern(newData, X, ncol(X))
  }
  
  if(whichkernel=="linear"){ K <- tcrossprod(rbind(newData, X)) }
  if(whichkernel=="poly2"){ K <- (tcrossprod(rbind(newData, X))+1)^2 }
  if(whichkernel=="poly3"){ K <- (tcrossprod(rbind(newData, X))+1)^3 }
  if(whichkernel=="poly4"){ K <- (tcrossprod(rbind(newData, X))+1)^4 }
  if(is.null(K) & is.null(newK)){ stop("No valid Kernel specified") }
  
  if(whichkernel != "gaussian"){
    newK <- matrix(K[1:nn, (nn+1):(nn+nrow(X))],
                   nrow=nrow(newData),
                   byrow=FALSE)
  }
  
  return(newK)
}

## Function that splits a vector in to n chunks
#' @export
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

#Function to multiply a square matrix, X, with a diagonal matrix, diag(d)
#' @export
multdiag <- function(X,d){	
  R=matrix(NA,nrow=dim(X)[1],ncol=dim(X)[2])		
  for (i in 1:dim(X)[2]){
    R[,i]=X[,i]*d[i]	
  }
  return(R)
}

## Lambda search for KRLS
#' @export
lambdasearch <-
  function(L=NULL,
           U=NULL,
           y=NULL,
           K=NULL,
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
    S1 <- looloss(lambda=X1,y=y,K=K, eigtrunc=eigtrunc,truncate=truncate)
    S2 <- looloss(lambda=X2,y=y,K=K, eigtrunc=eigtrunc,truncate=truncate)
    
    if(noisy){cat("L:",L,"X1:",X1,"X2:",X2,"U:",U,"S1:",S1,"S2:",S2,"\n") }
    
    while(abs(S1-S2)>tol){ # terminate if difference between S1 and S2 less than tolerance
      
      # update steps and use caching
      if(S1 < S2){
        U  <- X2
        X2 <- X1
        X1 <- L + (.381966)*(U-L)
        S2 <- S1
        S1 <- looloss(lambda=X1,y=y,K=K, eigtrunc=eigtrunc,truncate=truncate)
        
      } else { #S2 < S1
        L  <- X1
        X1 <- X2
        X2 <- U - (.381966)*(U-L)
        S1 <- S2
        S2 <- looloss(lambda=X2,y=y,K=K,eigtrunc=eigtrunc,truncate=truncate)
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
  function(y=NULL,K=NULL,Utrunc=NULL,lambda=NULL,eigtrunc=NULL,truncate=NULL){
    if (truncate) {
      return(solve_for_c_ls_trunc(y=y,K=K,Utrunc=Utrunc,lambda=lambda)$Le)
    } else {
      return(solve_for_c_ls(y=y,K=K,lambda=lambda)$Le)
    }
  }

