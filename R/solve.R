## This function contains functions to solve for coefficients
## Functions:
##            solveForC (used for krlogit; see solve_for_c_ls* in the cpp for krls)

#' @export
getDhat <- function(par = NULL,
                    y,
                    U,
                    D,
                    w,
                    lambda,
                    con = list(),
                    ctrl) {
  
  if(ctrl$loss == "leastsquares") {
    if(ctrl$weight) {
      out <- solve_for_d_ls_w(y=y, U=U, D=D, w=w, lambda=lambda)
    } else {
      out <- solve_for_d_ls(y=y, U=U, D=D, lambda=lambda)
    }
  } else if(ctrl$loss == "logistic") {
    out <- solveForDOptim(par = par,
                   y = y,
                   U = U,
                   D = D,
                   w = w,
                   lambda = lambda,
                   con = con,
                   printout = !ctrl$quiet,
                   returnopt = ctrl$returnopt)
  } else {
    stop("Loss must be either 'leastsquares' or 'logistic'.")
  }
  
  return(out)
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
solveForDOptim <- function(par= NULL,
                      y = NULL,
                      U = NULL,
                      D = NULL,
                      w = NULL,
                      lambda = NULL,
                      con = list(),
                      printout = FALSE,
                      returnopt = TRUE) {
  
  ncoeffs <- ncol(U)
  if(is.null(par)) {
    par <- rep(0, ncoeffs + 1)
  }
  
  opt <- optim(par,
               fn=krlogit_fn_trunc,
               gr=krlogit_gr_trunc,
               U=U,
               D=D,
               y=y,
               w=w,
               lambda=lambda,
               method="BFGS",
               control=con)
  
  dhat <- opt$par[1:ncoeffs]
  beta0hat <- opt$par[ncoeffs+1]
  
  if (!returnopt) opt <- NULL
  if (printout) {
    print(paste("Calls to function:", opt$counts[1], ", calls to gradient:", opt$counts[2]))
    print(paste("Converged? ", as.logical(1 - opt$convergence)))
  }
  
  return(list(dhat = dhat,
              beta0hat = beta0hat,
              opt = opt))
}