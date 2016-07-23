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
                    X,
                    y,
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
                    quiet = TRUE){
                    #cpp = TRUE,
                    #sandwich = ifelse(loss == "leastsquares", FALSE, TRUE),
                    #clusters = NULL) {
  
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
  
  # column names in case there are none
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
  
  ## Carry vars in list
  control = list(d=d,
                 loss=loss,
                 whichkernel=whichkernel,
                 truncate=truncate,
                 lastkeeper=lastkeeper,
                 epsilon=epsilon,
                 quiet=quiet)
  
  ###----------------------------------------
  ## Preparing to search for hyperparameters
  ###----------------------------------------
  
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
  
  ###----------------------------------------
  ## searching for hyperparameters
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
    yfitted <- yfitted*y.init.sd+y.init.mean
    
    beta0hat <- NULL
    
  } else {
    out <- solveForC(y=y, K=Kdat$K, Utrunc=Kdat$Utrunc, D=Kdat$eigvals, lambda=lambda, con=con, printout = printout, returnopt = returnopt)
    chat <- out$chat
    beta0hat <- out$beta0hat
    opt <- if(returnopt) out$opt else NULL
    
    ## getting c_p from d
    coefhat=NULL
    if(truncate==FALSE){coefhat=chat} else {
      UDinv = mult_diag(Kdat$Utrunc, 1/Kdat$eigvals)
      coefhat=UDinv %*% chat
    }
    
    yfitted <- logistic(K=Kdat$K, coeff=coefhat, beta0 = beta0hat)
    
  }
  
  z <- list(K=Kdat$K,
            Utrunc=Kdat$Utrunc,
            D=Kdat$eigvals,
            eigobj=Kdat$eigobj,
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
            loss=loss
  )
  
  class(z) <- "krls2"  
  
  return(z)
}


