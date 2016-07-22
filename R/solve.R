## This function contains functions to solve for coefficients
## Functions:
##            solveForC (used for krlogit; see solve_for_c_ls* in the cpp for krls)

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
                      printout = FALSE,
                      returnopt = TRUE) {
  
  ncoeffs <- ifelse(is.null(D), ncol(K), ncol(Utrunc))
  par <- rep(0, ncoeffs + 1)
  
  if(is.null(D)){
    opt <- optim(par, krlogit.fn, gr=krlogit.gr, K=K, y=y, lambda=lambda,
                 method="BFGS", control=con)
  } else {
    opt <- optim(par, krlogit_fn_trunc, gr=krlogit_gr_trunc, Utrunc=Utrunc, D = D, y=y, lambda=lambda,
                 method="BFGS", control=con)
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
              opt = opt))
}