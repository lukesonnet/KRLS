#----------
# Roxygen Commands
#----------
#' Main function of KRLS
#' @description This is the primary fitting function.
#' By default it uses squared loss
#' (loss="leastsquares") as one would for a
#' continuous outcome, but now also implements
#' logistic regression with the loss="logistic" option.
#' It also allows faster computation and larger
#' training sets than prior versions by an option
#' to approximate the kernel matrix with
#' lower dimensional approximation using the
#' truncate argument.
#'
#' The workflow for using KRLS mimics that
#' of lm and similar functions:
#' a krls object of class KRLS2 is fitted in one step,
#' then can later be examined using summary().
#' The krls object contains all the information
#' that may be needed at summary time,
#' including information required to estimate
#' pointwise partial derivatives, their average
#' for each covariate, standard errors, etc.
#' via the summary() function. See summary.krls2().
#'
#' @import RSpectra
#' @useDynLib KRLS2
#' @importFrom Rcpp sourceCpp


#############
# Functions #
#############

## This is the main function. Fits the model by preparing the data, searching
## for lambda through CV, minimizing the target function, and returning an
## object that contains all of the information about the fitting that is
## necessary for future analysis. I borrow liberally from the KRLS package

#' @description Function implements Kernel-Based Regularized Least Squares (KRLS), a machine learning method described in Hainmueller and Hazlett (2014) that allows users to solve regression and classification problems without manual specification search and strong functional form assumptions. KRLS finds the best fitting function by minimizing a Tikhonov regularization problem with a squared loss, using Gaussian Kernels as radial basis functions. KRLS reduces misspecification bias since it learns the functional form from the data. Yet, it nevertheless allows for interpretability and inference in ways similar to ordinary regression models. In particular, KRLS provides closed-form estimates for the predicted values, variances, and the pointwise partial derivatives that characterize the marginal effects of each independent variable at each data point in the covariate space. The distribution of pointwise marginal effects can be used to examine effect heterogeneity and or interactions.
#'
#'  @param X \emph{N} by \emph{D} data numeric matrix that contains the values of \emph{D} predictor variables for \eqn{i=1,\ldots,N} observations. The matrix may not contain missing values or constants. Note that no intercept is required for the least squares or logistic loss. In the case of least squares, the function operates on demeaned data and subtracting the mean of \emph{y} is equivalent to including an (unpenalized) intercept into the model. In the case of logistic loss, we automatically estimate an unpenalized intercept in the linear component of the model.
#'  @param y \var{N} by \var{1} data numeric matrix or vector that contains the values of the response variable for all observations. This vector may not contain missing values, and in the case of logistic loss should be a vector of \var{0}s and \var{1}s.
#'  @param w \var{N} by \var{1} data numeric matrix or vector that contains the weights that should applied to each observation. These need not sum to one.
#'  @param loss String vector that specifies the loss function. For KRLS, use \code{leastsquares} and for KRLogit, use \code{logistic}.
#'  @param whichkernel String vector that specifies which kernel should be used. Must be one of \code{gaussian}, \code{linear}, \code{poly1}, \code{poly2}, \code{poly3}, or \code{poly4} (see details). Default is \code{gaussian}.
#'  @param b A positive scalar (formerly \code{sigma}) that specifies the bandwidth of the Gaussian kernel (see \code{\link{gausskernel}} for details). By default, the bandwidth is set equal to \var{2D} (twice the number of dimensions) which typically yields a reasonable scaling of the distances between observations in the standardized data that is used for the fitting.
#'  @param optimb A boolean that if \code{TRUE} will numerically estimate \code{b} using cross validation error instead of setting \code{b} to the default.
#'  @param lambda A positive scalar that specifies the \eqn{\lambda}{lambda} parameter for the regularizer (see details). It governs the tradeoff between model fit and complexity. By default, this parameter is chosen by minimizing the sum of the squared leave-one-out errors for KRLS and by minimizing the sum of squared cross-validation errors for KRLogit, with the number of folds set by \code{hyperfolds}. When using logistic loss, \code{lambda} can also be a numeric vector of positive scalars in which case a linesearch over these values will be used to choose \code{lambda}.
#'  @param hyperfolds A positive scalar that sets the number of folds used in selecting \code{lambda} via cross-validation error.
#'  @param lambdastart A positive scalar that specifices the starting value for a numerical optimization of \code{lambda}. Only is used when \code{lambda} is \code{NULL} and with logistic loss.
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
                    quiet = TRUE,
                    sigma = NULL){
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
  if (is.null(w)){ #| all(w==1)) {
    w <- rep(1, n)
    weight = F
  } else if (length(w) != n) {
    stop("w is not the same length as y")
  } else {
    if(loss=="leastsquares" & !truncate)
      stop("For now, weighted KRLS only works with truncation")
    w <- n * (w / sum(w))
    weight = T
  }
  ## Carry vars in list
  control = list(w=w,
                 weight=weight,
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

  ## Legacy support for KRLS code with fixed sigma
  if(!is.null(sigma) & is.null(b)) {
    b <- sigma
  }

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
