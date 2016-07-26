## This file contains functions that are now deprecated




## The functions below are not used because their cpp alternatives are preferred
## for speed

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


# Function to multiply a square matrix, X, with a diagonal matrix, diag(d)
# Note that we have replaced this with the cpp counterpart mult_diag
#' @export
multdiag <- function(X,d){	
  R=matrix(NA,nrow=dim(X)[1],ncol=dim(X)[2])		
  for (i in 1:dim(X)[2]){
    R[,i]=X[,i]*d[i]	
  }
  return(R)
}



//' @export
//  [[Rcpp::export]]
Rcpp::List solve_for_c_lst(const arma::vec& y,
                          const arma::mat& K,
                          const double& lambda) {
  
  int nn =  y.n_elem;
  arma::mat Ginv(nn, nn);
  
  Ginv = arma::inv_sympd(K + lambda * arma::eye(nn, nn));

  arma::vec coeffs = Ginv * y;
  arma::vec tempLoss = coeffs / diagvec(Ginv);
  double Le = as_scalar(tempLoss.t() * tempLoss);
  
  return Rcpp::List::create(Rcpp::Named("coeffs") = coeffs,
                            Rcpp::Named("Le") = Le);
}

