#----------
# Roxygen Commands
#----------

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
                    w = NULL, # weights
                    # Family
                    loss = "leastsquares",
                    # Kernel arguments
                    whichkernel = "gaussian", 
                    b = NULL,
                    optimb = FALSE,
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
  
  ## todo: it might be faster to have different versions with and without weights
  ## because it saves a bunch of pointless multiplications. However this will
  ## increase the amount of code.
  if (is.null(w)) {
    w <- rep(1, n)
  } else if (length(w) != n) {
    stop("w is not the same length as y")
    w <- n * (w / sum(w))
  }
  ## Carry vars in list
  control = list(w=w,
                 d=d,
                 loss=loss,
                 whichkernel=whichkernel,
                 truncate=truncate,
                 lastkeeper=lastkeeper,
                 epsilon=epsilon,
                 quiet=quiet,
                 returnopt=returnopt)
  
  ###----------------------------------------
  ## Preparing to search for hyperparameters
  ###----------------------------------------
  
  ## Initialize hyper-parameter variables
  lambdarange <- NULL
  brange <- NULL
  
  ## Default b to be 2*the number of features
  if(is.null(b)){
    if (!optimb) {
      b <- 2*d
    } else if (loss == "leastsquares") {
      message("Warning: Cannot simultaneously search for both lambda and b with leastsquares loss.")
      message(sprintf("Setting b to %s", 2*d))
      b <- 2*d
    }
  } else{
    if(length(b) > 1) {
      brange <- b
      b <- NULL
    } else {
      stopifnot(is.vector(b),
                is.numeric(b),
                b>0)
    }
    if(optimb) warning("optimb ignored when b argument is used.")
    
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
  if(is.null(lambda) | is.null(b)) {
    chunks <- chunk(sample(n), hyperfolds)
    
    hyperctrl <- list(chunks = chunks,
                      lambdastart = lambdastart,
                      lambdarange = lambdarange,
                      brange = brange,
                      optimb = optimb,
                      Lbound = L,
                      Ubound = U)
  } else {
    chunks <- NULL
  }
  
  ###----------------------------------------
  ## searching for hyperparameters
  ###----------------------------------------

  ## If b ix fixed, compute all of the kernels
  if(!is.null(b)) {
    
    ## Compute kernel matrix and truncated objects
    Kdat <- generateK(X=X,
                      b=b,
                      control=control)
    
    if(is.null(lambda)) {
      lambda <- lambdasearch(y=y,
                             X=X,
                             Kdat=Kdat,
                             hyperctrl=hyperctrl,
                             control=control,
                             b=b)
    }
  } else { # if(is.null(b)
    
    if(is.null(lambda)) {
      
      hyperOut <- lambdabsearch(y=y,
                                    X=X,
                                    hyperctrl=hyperctrl,
                                    control=control)
      lambda <- hyperOut$lambda
      b <- hyperOut$b
      
    } else { # lambda is set
      b <- bsearch(y=y,
                           X=X,
                           hyperctrl=hyperctrl,
                           control=control,
                           lambda=lambda)
    }
    
    Kdat <- generateK(X=X,
                      b=b,
                      control=control)
    
  }
  
  message(sprintf("Lambda selected: %s", lambda))
  message(sprintf("b selected: %s", b))
  
  ###----------------------------------------
  ## Estimating choice coefficients (solving)
  ###----------------------------------------

  out <- getDhat(y = y,
                 U = Kdat$U,
                 D = Kdat$D,
                 w = control$w,
                 lambda = lambda,
                 con = con,
                 ctrl = control)

  UDinv <- mult_diag(Kdat$U, 1/Kdat$D)

  coefhat <- UDinv %*% out$dhat
  
  # todo: replace with internal predict dispatcher that takes getDhat object
  if (loss == "leastsquares") {

    yfitted <- Kdat$K%*%coefhat
    yfitted <- yfitted*y.init.sd+y.init.mean
  
  } else {
    
    opt <- if(returnopt) out$opt else NULL
    
    yfitted <- logistic(K=Kdat$K, coeff=coefhat, beta0 = out$beta0hat)
    
  }
  
  z <- list(K=Kdat$K,
            U=Kdat$U,
            D=Kdat$D,
            w=control$w,
            lastkeeper = Kdat$lastkeeper,
            truncate = truncate,
            coeffs=coefhat,
            dhat=out$dhat,
            beta0hat = out$beta0hat,
            fitted=yfitted,
            X=X.init,
            y=y.init,
            b=b,
            lambda=lambda,
            kernel=whichkernel,
            loss=loss
  )
  
  class(z) <- "krls2"  
  
  return(z)
}


# check ncol U reliance