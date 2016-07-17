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
                    sandwich = ifelse(loss == "leastsquares", FALSE, TRUE),
                    clusters = NULL) {
  
  ###----------------------------------------
  ## Input validation and housekeeping
  ###----------------------------------------
  
  ## Prepare the data
  X <- as.matrix(X)
  y <- as.matrix(y)
  y.init <- y
  
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
  
  stopifnot(
    is.logical(derivative),
    is.logical(vcov)
    # todo: add more checks here
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

    y.init.sd <- apply(y.init,2,sd)
    y.init.mean <- mean(y.init)
    y <- scale(y,center=y.init.mean,scale=y.init.sd)
    
  } else if (loss == "logistic") {
    if (!all(y %in% c(0,1))) {
      stop("y is not binary data")
    }
  }
  
  ## Initialize hyper-parameter variables
  lambdarange <- NULL
  sigmarange <- NULL
  
  ## Default sigma to the number of features
  if(is.null(sigma)){
    if (!optimsigma) {sigma <- d}
  } else{
    if(length(sigma) > 1) {
      sigmarange <- sigma
      sigma <- NULL
    } else {
      stopifnot(is.vector(sigma),
                is.numeric(sigma),
                sigma>0)
    }
    if(optimsigma) warning("optimsigma ignored when sigma argument is used.")
    
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
  if(is.null(lambda) | is.null(sigma)) {
    chunks <- chunk(sample(n), hyperfolds)

    hyperctrl <- list(chunks = chunks,
                      lambdastart = lambdastart,
                      lambdarange = lambdarange,
                      sigmarange = sigmarange,
                      optimsigma = optimsigma,
                      L = L,
                      U = U)
  } else {
    chunks <- NULL
  }
  
  ## Carry vars in list
  control = list(d=d,
                 loss=loss,
                 whichkernel=whichkernel,
                 truncate=truncate,
                 lastkeeper=lastkeeper,
                 epsilon=epsilon,
                 quiet=FALSE)
  
  ###----------------------------------------
  ## Setting hyperparameters sigma and lambda
  ###----------------------------------------
  
  ## If sigma ix fixed, compute all of the kernels
  if(!is.null(sigma)) {
    
    ## Compute kernel matrix and truncated objects
    Kdat <- generateK(X=X,
                     sigma=sigma,
                     control=control)
    
    if(is.null(lambda)) {
      lambda <- lambdasearch(y=y,
                             X.init=X.init,
                             Kdat=Kdat,
                             hyperctrl=hyperctrl,
                             control=control,
                             sigma=sigma)
    }
  } else { # if(is.null(sigma)
    if(loss == "leastsquares") {
      message("Warning: Cannot simultaneously search for both lambda and sigma with leastsquares loss.")
      message(sprintf("Setting sigma to %s", d))
      sigma <- d
    }
    
    if(is.null(lambda)) {

      hyperOut <- lambdasigmasearch(y=y,
                                    X.init=X.init,
                                    hyperctrl=hyperctrl,
                                    control=control)
      lambda <- hyperOut$lambda
      sigma <- hyperOut$sigma
      
    } else { # lambda is set
      sigma <- sigmasearch(y=y,
                           X.init=X.init,
                           hyperctrl=hyperctrl,
                           control=control,
                           lambda=lambda)
    }

    Kdat <- generateK(X=X,
                      sigma=sigma,
                      control=control)
    
  }

  message(sprintf("Lambda selected: %s", lambda))
  message(sprintf("Sigma selected: %s", sigma))
  
  ###----------------------------------------
  ## Estimating choice coefficients (solving)
  ###----------------------------------------
  
  vcovmatc=NULL
  score <- NULL
  hessian <- NULL
  vcov.cb0 <- NULL
  if (loss == "leastsquares") {
    if (truncate){
      out <- solve_for_c_ls_trunc(y=y, D = Kdat$eigvals, Utrunc =Kdat$Utrunc, lambda=lambda)
      chat <- out$coeffs
    } else {
      out <- solve_for_c_ls(y=y, K=Kdat$K, lambda=lambda)
      chat <- out$coeffs
    }
    
    ## getting c_p from d
    coefhat=NULL
    if(truncate==FALSE){coefhat=chat} else {
      UDinv = mult_diag(Kdat$Utrunc, 1/Kdat$eigvals)
      coefhat=UDinv %*% chat
    }
    
    yfitted <- Kdat$K%*%coefhat
    
    if (vcov) {
      sigmasq <- as.vector((1/n) * crossprod(y-yfitted))


      vcovmatc <- tcrossprod(mult_diag(Kdat$eigobj$vectors,sigmasq*(Kdat$eigobj$values+lambda)^-2),Kdat$eigobj$vectors)   
    }
    beta0hat <- NULL
  } else {
    out <- solveForC(y=y, K=Kdat$K, Utrunc=Kdat$Utrunc, D=Kdat$eigvals, lambda=lambda, con=con, vcov=vcov, printout = printout, returnopt = returnopt)
    chat <- out$chat
    beta0hat <- out$beta0hat
    hessian <- out$hessian
    opt <- if(returnopt) out$opt else NULL
    
    ## getting c_p from d
    coefhat=NULL
    if(truncate==FALSE){coefhat=chat} else {
      UDinv = mult_diag(Kdat$Utrunc, 1/Kdat$eigvals)
      coefhat=UDinv %*% chat
    }
    
    yfitted <- logistic(K=Kdat$K, coeff=coefhat, beta0 = beta0hat)
    
    ## Vcov of d
    if (vcov & truncate) {
      #optimHess = hessian
      #fixedHess = krlogit_hess_trunc(c(chat, beta0hat), Kdat$Utrunc, Kdat$eigvals, y, lambda)
      #fixedHessR <- krlogit.hess.trunc(c(chat, beta0hat), Kdat$Utrunc, Kdat$eigvals, y, lambda)
      vcov.cb0 = solve(hessian)
      
      score <- matrix(nrow = n, ncol = length(c(chat,beta0hat)))
      B <- matrix(0, nrow = length(c(chat,beta0hat)), ncol = length(c(chat,beta0hat)))
      
      for(i in 1:n) {
        score[i, ] = t(-1 * krlogit_gr_trunc(par=c(chat,beta0hat), Kdat$Utrunc[i, , drop=F], Kdat$eigvals, y[i, drop = F], lambda/n))
        B <- B + tcrossprod(score[i, ])
      }
      if(sandwich) {
        if(!is.null(clusters)) {
          B <- matrix(0, nrow = length(c(chat,beta0hat)), ncol = length(c(chat,beta0hat)))
          for(j in 1:length(clusters)){
            B <- B + tcrossprod(apply(score[clusters[[j]], ], 2, sum))
          }
          
        }
        vcov.cb0 <- vcov.cb0 %*% B %*% vcov.cb0
        
      }
      vcovmatc = tcrossprod(UDinv%*%vcov.cb0[1:length(chat), 1:length(chat)], UDinv)
      ## todo: return var b0
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

    derivout <- pwmfx(Kdat$K, X, coefhat, vcovmatc, p1p0, sigma)
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
      L <-  distk*Kdat$K  #Kfull rather than K here as truncation handled through coefs
      derivmat[,k] <- p1p0*(-1/sigma)*(L%*%coefhat)
      if(truncate==FALSE) {
        varavgderivmat = NULL
      } else {
        varavgderivmat[1,k] <- 1/(sigma * n)^2 * sum(crossprod(p1p0^2, crossprod(L,vcovmatc%*%L)))
      }
    }
    
  }

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
  
  #score <- matrix(nrow = n, ncol = length(c(chat,beta0hat)))
  #B <- matrix(0, nrow = length(c(chat,beta0hat)), ncol = length(c(chat,beta0hat)))
  #if (truncate & loss == 'logistic') {
  #  for(i in 1:n) {
  #    score[i, ] = t(-1 * krlogit_gr_trunc(par=c(chat,beta0hat), Kdat$Utrunc[i, , drop=F], Kdat$eigvals, y[i, drop = F], lambda/n))
  #    B <- B + tcrossprod(score[i, ])
  #  }
  #} else {
  #  score = -2 * K %*% (y - yfitted) + 2 * lambda * yfitted
  #}
  #if {score = krlogit.gr(par=c(chat,beta0hat), K=K, y=y,lambda=lambda)}
  
  
  #vcov.cb0.sand <- vcov.cb0 %*% B %*% vcov.cb0
  
  ###----------------------------------------
  ## Returning results
  ###----------------------------------------
  
  z <- list(K=Kdat$K,
            Utrunc=Kdat$Utrunc,
            lastkeeper = lastkeeper,
            truncate = truncate,
            coeffs=coefhat,
            chat=chat,
            beta0hat = beta0hat,
            fitted=yfitted,
            X=X.init,
            y=y.init,
            sigma=sigma,
            lambda=lambda,
            kernel=whichkernel,
            derivmat = derivmat,
            avgderiv = avgderiv,
            var.avgderiv = varavgderivmat,
            loss=loss,
            score=score,
            hessian=hessian,
            #fixedHessR=fixedHessR
            #lastkeeper = ncol(Ktilde)  
            #score = score
            vcov.cb0 = vcov.cb0
            #vcov.cb0.sand = vcov.cb0.sand
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

## This has the wrong norm
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

## this has the wrong norm
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

#' @export
krlogit.hess.trunc <- function(par, Utrunc, D, y, lambda) {
  
  coef <- par[1:ncol(Utrunc)]
  beta0 <- par[ncol(Utrunc)+1]
  
  Ud <- Utrunc %*% coef
  
  meat <- exp(-Ud - beta0) / (1 + exp(-Ud - beta0))^2
  
  dcdc <- mult_diag(t(Utrunc), meat) %*% Utrunc + diag(2 * lambda / D)
  dcdb <- crossprod(Utrunc, meat)
  dbdb <- sum(meat)

  print(dcdc[1:3, 1:3])
  
  ret = matrix(nrow = length(par), ncol = length(par))
  ret[1:length(coef), 1:length(coef)] <- dcdc
  ret[length(par), 1:length(coef)] <- dcdb
  ret[1:length(coef), length(par)] <- dcdb
  ret[length(par), length(par)] <- dbdb
  
  print(ret[1:3, 1:3])
  
  return(ret)
}

#krlogit.hess <- function(par, K, y, lambda){
#  cb <- crossprod(K, (exp(-K%*%chat - beta0)/(1 + exp(-K%*%chat - beta0))^2))
#  bb <- sum(exp(-K%*%chat - beta0)/(1 + exp(-K%*%chat - beta0))^2)
#  cc <- can't get this to match what optim returns, see krlogit_testingmath
  



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
  
  if(object$loss == "logistic") {
    yfitted <- logistic(K = newdataK, coeff = object$coeffs, beta0 = object$beta0hat)
  } else if (object$loss == "leastsquares") {
    yfitted <- ((newdataK %*% object$coeffs)* apply(object$y,2,sd))+mean(object$y)
  }
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
      return(solve_for_c_ls_trunc(y=y,D=D,Utrunc=Utrunc,lambda=lambda)$Le)
    } else {
      return(solve_for_c_ls(y=y,K=K,lambda=lambda)$Le)
    }
  }

