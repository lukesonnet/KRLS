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
#' @useDynLib KRLS2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics arrows par plot
#' @importFrom stats as.formula optim optimize predict pt quantile sd var

#############
# Functions #
#############

## This is the main function. Fits the model by preparing the data, searching
## for lambda through CV, minimizing the target function, and returning an
## object that contains all of the information about the fitting that is
## necessary for future analysis. I borrow liberally from the KRLS package

#' @param X \emph{N} by \emph{D} data numeric matrix that contains the values of \emph{D} predictor variables for \eqn{i=1,\ldots,N} observations. The matrix may not contain missing values or constants. Note that no intercept is required for the least squares or logistic loss. In the case of least squares, the function operates on demeaned data and subtracting the mean of \emph{y} is equivalent to including an (unpenalized) intercept into the model. In the case of logistic loss, we automatically estimate an unpenalized intercept in the linear component of the model.
#' @param y \var{N} by \var{1} data numeric matrix or vector that contains the values of the response variable for all observations. This vector may not contain missing values, and in the case of logistic loss should be a vector of \var{0}s and \var{1}s.
#' @param w \var{N} by \var{1} data numeric matrix or vector that contains the weights that should applied to each observation. These need not sum to one.
#' @param loss String vector that specifies the loss function. For KRLS, use \code{leastsquares} and for KRLogit, use \code{logistic}.
#' @param whichkernel String vector that specifies which kernel should be used. Must be one of \code{gaussian}, \code{linear}, \code{poly1}, \code{poly2}, \code{poly3}, or \code{poly4} (see details). Default is \code{gaussian}.
#' @param b A positive scalar (formerly \code{sigma}) that specifies the bandwidth of the Gaussian kernel (see \code{\link{gausskernel}} for details). By default, the bandwidth is set equal to \var{2D} (twice the number of dimensions) which typically yields a reasonable scaling of the distances between observations in the standardized data that is used for the fitting. You can also pass a numeric vector to do a grid search over possible b values.
#' @param bstart A positive scalar that is the starting value for a numerical estimation of the \code{b} parameter using cross validation error, overriding the default. If \code{b} is specified as an argument, \code{bstart} is ignored.
#' @param binterval A numeric vector of length two that specifies the minimum and maxmum \code{b} values to look over with \code{optimize} when optimizing only the \code{b} hyperparameter. Both values must be strictly positive. Only used with logistic loss and when \code{b} is NULL. Defaults to \code{c(10^-8, 500*p)}. This is for use with numerical optimization, if you want to do a grid search, instead pass a numerical vector to the \code{lambda} argument.
#' @param lambda A positive scalar that specifies the \eqn{\lambda}{lambda} parameter for the regularizer (see details). It governs the tradeoff between model fit and complexity. By default, this parameter is chosen by minimizing the sum of the squared leave-one-out errors for KRLS and by minimizing the sum of cross-validation negative log likelihood for KRLogit, with the number of folds set by \code{hyperfolds}. When using logistic loss, \code{lambda} can also be a numeric vector of positive scalars in which case a line search over these values will be used to choose \code{lambda}.
#' @param hyperfolds A positive scalar that sets the number of folds used in selecting \code{lambda} or \code{b} via cross-validation.
#' @param lambdastart A positive scalar that specifices the starting value for a numerical optimization of \code{lambda}. Only is used when jointly optimizing over \code{lambda} and \code{b}
#' @param lambdainterval A numeric vector of length two that specifies the minimum and maxmum \code{lambda} values to look over with \code{optimize}. Both values must be strictly positive. Only used with logistic loss and when \code{lambda} is NULL. Defaults to \code{c(10^-8, 25)}. This is for use with numerical optimization, if you want to do a grid search, instead pass a numerical vector to the \code{lambda} argument.
#' @param L Non-negative scalar that determines the lower bound of the search window for the leave-one-out optimization to find \eqn{\lambda}{lambda} with least squares loss. Default is \code{NULL} which means that the lower bound is found by using an algorithm outlined in \code{\link{lambdaline}}. Ignored with logistic loss.
#' @param U Positive scalar that determines the upper bound of the search window for the leave-one-out optimization to find \eqn{\lambda}{lambda} with least squares loss. Default is \code{NULL} which means that the upper bound is found by using an algorithm outlined in \code{\link{lambdaline}}. Ignored with logistic loss.
#' @param tol Positive scalar that determines the tolerance used in the optimization routine used to find \eqn{\lambda}{lambda} with least squares loss. Default is \code{NULL} which means that convergence is achieved when the difference in the sum of squared leave-one-out errors between the \var{i} and the \var{i+1} iteration is less than \var{N * 10^-3}. Ignored with logistic loss.
#' @param truncate A boolean that defaults to \code{FALSE}. If \code{TRUE} truncates the kernel matrix, keeping as many eigenvectors as needed so that 1-\code{epsilon} of the total variance in the kernel matrix is retained. Alternatively, you can simply specify \code{epsilon} and truncation will be used.
#' @param epsilon Scalar between 0 and 1 that determines the total variance that can be lost in truncation. If not NULL, truncation is automatically set to TRUE. If \code{truncate == TRUE}, default is 0.001.
#' @param lastkeeper Number of columns of \code{U} to keep when \code{truncate == TRUE}. Overrides \code{epsilon}.
#' @param con A list of control arguments passed to optimization for the numerical optimization of the kernel regularized logistic loss function.
#' @param returnopt A boolean that defaults to \code{FALSE}. If \code{TRUE}, returns the result of the \code{optim} method called to optimize the kernel regularized logistic loss function. Returns \code{NULL} with leastsquares loss.
#' @param printlevel A number that is either 0 (default), 1, or 2. 0 Has minimal printing, 1 prints out most diagnostics, and 2 prints out most diagnostics including \code{optim} diagnostics for each fold in the cross-validation selection of hyperparameters.
#' @param warn A number that sets your \code{warn} option. We default to 1 so that warnings print as they occur. You can change this to 2 if you want all warnings to be errors, to 0 if you want all warnings to wait until the top-level call is finished, or to a negative number to ignore them.
#' @param sigma DEPRECATED. Users should now use \code{b}, included for backwards compatability.
#' @details
#' \code{krls} implements the Kernel-based Regularized Least Squares (KRLS) estimator as described in Hainmueller and Hazlett (2014). Please consult this reference for any details.

#' Kernel-based Regularized Least Squares (KRLS) arises as a Tikhonov minimization problem with a squared loss. Assume we have data of the from \eqn{y_i,\,x_i}{y_i, x_i} where \var{i} indexes observations, \eqn{y_i \in R}{y_i in R} is the outcome and \eqn{x_i \in R^D}{x_i in R^D} is a \var{D}-dimensional vector of predictor values. Then KRLS searches over a space of functions \eqn{H} and chooses the best fitting function \eqn{f} according to the rule:
#'
#' \deqn{argmin_{f \in H} \sum_i^N (y_i - f(x_i))^2 + \lambda ||f||_{H^2}}{%
#'       argmin_{f in H} sum_i^N (y_i - f(x_i))^2 + lambda || f ||_H^2}
#'
#' where \eqn{(y_i - f(x_i))^2} is a loss function that computes how `wrong' the function is at each observation \var{i} and \eqn{|| f ||_{H^2}}{|| f ||_H^2} is the regularizer that measures the complexity of the function according to the \eqn{L_2} norm \eqn{||f||^2 = \int f(x)^2 dx}{||f||^2 = int f(x)^2 dx}. \eqn{\lambda}{lambda} is the scalar regularization parameter that governs the tradeoff between model fit and complexity. By default, \eqn{\lambda}{lambda} is chosen by minimizing the sum of the squared leave-one-out errors, but it can also be specified by the user in the \code{lambda} argument to implement other approaches.
#'
#' Under fairly general conditions, the function that minimizes the regularized loss within the hypothesis space established by the choice of a (positive semidefinite) kernel function \eqn{k(x_i,x_j)} is of the form
#'
#' \deqn{f(x_j)= \sum_i^N c_i k(x_i,x_j)}{%
#'       f(x_j)= sum_i^N c_i k(x_i,x_j)}
#'
#' where the kernel function \eqn{k(x_i,x_j)} measures the distance between two observations \eqn{x_i} and \eqn{x_j} and \eqn{c_i} is the choice coefficient for each observation \eqn{i}. Let \eqn{K} be the \eqn{N} by \eqn{N} kernel matrix with all pairwise distances \eqn{K_ij=k(x_i,x_j)} and \eqn{c} be the  \eqn{N} by \eqn{1} vector of choice coefficients for all observations then in matrix notation the space is \eqn{y=Kc}.
#'
#' Accordingly, the \code{krls} function solves the following minimization problem
#'
#' \deqn{argmin_{f \in H} \sum_i^n (y - Kc)'(y-Kc)+ \lambda c'Kc}{%
#'       argmin_{f in H} sum_i^n (y - Kc)'(y-Kc)+ lambda c'Kc}
#'
#' which is convex in \eqn{c} and solved by \eqn{c=(K +\lambda I)^-1 y}{c=(K +lambda I)^-1 y} where \eqn{I} is the identity matrix. Note that this linear solution provides a flexible fitted response surface that typically reduces misspecification bias because it can learn a wide range of nonlinear and or nonadditive functions of the predictors. In an extension, Hazlett and Sonnet consier a logistic loss function, details of which are forthcoming.
#'
#' If \code{vcov=TRUE} is specified, \code{krls} also computes the variance-covariance matrix for the choice coefficients \eqn{c} and fitted values \eqn{y=Kc} based on a variance estimator developed in Hainmueller and Hazlett (2014). Note that both matrices are \var{N} by \var{N} and therefore this results in increased memory and computing time.
#'
#' By default, \code{krls} uses the Gaussian Kernel (\code{whichkernel = "gaussian"}) given by
#'
#' \deqn{k(x_i,x_j)=exp(\frac{-|| x_i - x_j ||^2}{\sigma^2})}{%
#'       k(x_i,x_j)=exp(-|| x_i - x_j ||^2 / sigma^2)}
#'
#' where \eqn{||x_i - x_j||} is the Euclidean distance. The kernel bandwidth \eqn{\sigma^2}{sigma^2} is set to \eqn{D}, the number of dimensions, by default, but the user can also specify other values using the \code{sigma} argument to implement other approaches.
#'
#' If \code{binary=TRUE} is also specified, the function will identify binary predictors and return first differences for these predictors instead of partial derivatives. First differences are computed going from the minimum to the maximum value of each binary predictor. Note that first differences are more appropriate to summarize the effects for binary predictors (see Hainmueller and Hazlett (2014) for details).
#'
#' A few other kernels are also implemented, but derivatives are currently not supported for these: "linear": \eqn{k(x_i,x_j)=x_i'x_j}, "poly1", "poly2", "poly3", "poly4" are polynomial kernels based on  \eqn{k(x_i,x_j)=(x_i'x_j +1)^p} where \eqn{p} is the order.
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
                    bstart = NULL,
                    binterval = c(10^-8, 500*ncol(X)),
                    # Lambda arguments
                    lambda = NULL,
                    hyperfolds = 5,
                    lambdastart = 10^(-6)/length(y),
                    lambdainterval = c(10^-8, 25),
                    L = NULL,
                    U = NULL,
                    tol = NULL,
                    # Truncation arguments
                    truncate = FALSE,
                    epsilon = NULL,
                    lastkeeper = NULL,
                    # Optimization arguments
                    con = list(maxit=500),
                    returnopt = FALSE,
                    printlevel = 0,
                    warn = 1,
                    sigma = NULL, # to provide legacy support for old code,
                                  # simply is interpreted as 'b' if 'b' is NULL;
                                  # ignored otherwise
                    ...) {

  ###----------------------------------------
  ## Input validation and housekeeping
  ###----------------------------------------

  # set R to issue warnings as they occur if default = 1
  options(warn=warn)
  
  ## Prepare the data
  X <- as.matrix(X)
  y <- as.matrix(y)
  y.init <- y

  n <- nrow(X)
  d <- ncol(X)

  if (!is.numeric(X)) {
    stop("X must be numeric")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (var(y)==0) {
    stop("y is a constant")
  }
  if (anyNA(X)){
    stop("X contains missing data")
  }
  if (anyNA(y)){
    stop("y contains missing data")
  }
  if (n != nrow(y)){
    stop("nrow(X) not equal to number of elements in y")
  }

  # column names in case there are none
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", 1:d)
  }

  ## Scale data
  X.init <- X
  X.init.sd <- apply(X.init, 2, sd)
  if (any(X.init.sd == 0)) {
    stop("at least one column (", which(X.init.sd == 0) ,") in X is a constant, please remove the constant(s)")
  }
  X <- scale(X, center = TRUE, scale = X.init.sd)

  if (loss == "leastsquares") {
    y.init.sd <- apply(y.init, 2, sd)
    y.init.mean <- mean(y.init)
    y <- scale(y, center=y.init.mean, scale=y.init.sd)

  } else if (loss == "logistic") {
    if (!all(y %in% c(0,1))) {
      stop("y is not binary data")
    }
  }

  ## Set truncation options
  if (!is.null(epsilon)) {
    if (!is.numeric(epsilon) || epsilon < 0 || epsilon > 1) {
      stop("epsilon must be numeric and between 0 and 1, inclusive")
    } else {
      truncate <- TRUE
    }
  } else if (truncate) {
    epsilon = 0.001 # If truncate but no epsilon specified, use default
  }
  
  ## Warn if not truncating a big dataset
  if (n >= 1000 && !truncate) {
    warning("With n >= 1000 you should consider using truncation for speed. Try setting epsilon to 0.001")
  }
  
  ## todo: it might be faster to have different versions with and without weights
  ## because it saves a bunch of pointless multiplications. However this will
  ## increase the amount of code.
  if (is.null(w)) { 
    w <- rep(1, n)
    weight <- FALSE
  } else if (length(w) != n) {
    stop("w is not the same length as y")
  } else {
    if(loss == "leastsquares" & !truncate) {
      stop("For now, weighted KRLS only works with truncation")
    }
    w <- n * (w / sum(w)) # rescale w to sum to n
    weight <- TRUE
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
                 printlevel=printlevel,
                 returnopt=returnopt,
                 ...)

  ###----------------------------------------
  ## Preparing to search for hyperparameters
  ###----------------------------------------

  ## Initialize hyper-parameter variables
  lambdarange <- NULL
  brange <- NULL

  ## Legacy support for KRLS code with fixed sigma
  if (!is.null(sigma) & is.null(b)) {
    b <- sigma
  }

  ## Default b to be 2*the number of features
  if (is.null(b)) {
    if (is.null(bstart)) {
      b <- 2*d
    } else if (loss == "leastsquares") {
      b <- 2*d
      warning(sprintf("Cannot simultaneously search for both lambda and b with leastsquares loss.
                       Setting b to %s", b))
    }
  } else {
    if(length(b) > 1) {
      brange <- b
      b <- NULL
    } else {
      stopifnot(is.vector(b),
                is.numeric(b),
                b>0)
    }
    if(!is.null(bstart)) warning("bstart ignored when b argument is used.")

  }


  ## todo: unpack truncdat and add checks for if truncate is true, do not have
  ## two seperate processes but instead checks at pertinent steps of process and
  ## try to use similar variable names to economize on code

  ## todo: build distance matrix E first so that we can more quickly compute
  ## K below

  if (!is.null(lambda)) {
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
  if (is.null(lambda) || is.null(b)) {
    chunks <- chunk(sample(n), hyperfolds)

    control <- c(control,
                    list(chunks = chunks,
                         lambdastart = lambdastart,
                         lambdarange = lambdarange,
                         lambdainterval = lambdainterval,
                         brange = brange,
                         bstart = bstart,
                         binterval = binterval,
                         Lbound = L,
                         Ubound = U,
                         tol = tol))
  } else {
    chunks <- NULL
  }

  ###----------------------------------------
  ## searching for hyperparameters
  ###----------------------------------------

  ## If b ix fixed, compute all of the kernels
  if (!is.null(b)) {

    ## Compute kernel matrix and truncated objects
    Kdat <- generateK(X=X,
                      b=b,
                      control=control)

    if (is.null(lambda)) {
      lambda <- lambdasearch(
        y = y,
        X = X,
        Kdat = Kdat,
        control = control,
        b = b
      )
    }
  } else { # if(is.null(b)

    if (is.null(lambda)) {

      hyperOut <- lambdabsearch(y=y,
                                X=X,
                                control=control)
      lambda <- hyperOut$lambda
      b <- hyperOut$b

    } else { # lambda is set

      b <- bsearch(y=y,
                   X=X,
                   control=control,
                   lambda=lambda)
    }

    Kdat <- generateK(X=X,
                      b=b,
                      control=control)

  }

  message("Lambda selected: ", lambda)
  message("b selected: ", b)

  ###----------------------------------------
  ## Estimating choice coefficients (solving)
  ###----------------------------------------

  out <- getDhat(y = y,
                 U = Kdat$U,
                 D = Kdat$D,
                 w = control$w,
                 lambda = lambda,
                 con = con,
                 ctrl = control,
                 printopt = control$printlevel > 0)

  UDinv <- mult_diag(Kdat$U, 1/Kdat$D)

  coefhat <- UDinv %*% out$dhat

  opt <- NULL
  Le <- NULL
  if (loss == "leastsquares") {

    yfitted <- Kdat$K %*% coefhat
    yfitted <- yfitted * y.init.sd + y.init.mean
    Le <- out$Le
    
  } else {

    yfitted <- logistic(K=Kdat$K, coeff=coefhat, beta0 = out$beta0hat)

    if (returnopt) {
      opt <- out$opt
    }
  }

  z <- list(K = Kdat$K,
            U = Kdat$U,
            D = Kdat$D,
            w = control$w,
            lastkeeper = Kdat$lastkeeper,
            truncate = truncate,
            coeffs = coefhat,
            dhat = out$dhat,
            beta0hat = out$beta0hat,
            fitted = yfitted,
            X = X.init,
            y = y.init,
            b = b,
            lambda = lambda,
            kernel = whichkernel,
            loss = loss,
            opt = opt
  )
  
  if (is.numeric(Le)) {
    z[["Looe"]] <- Le * y.init.sd
  }

  class(z) <- "krls2"

  return(z)
}
