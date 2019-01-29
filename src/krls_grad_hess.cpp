#include <RcppArmadillo.h>
#include "krls_helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec krls_gr_trunc(
    const arma::mat& U,
    const arma::vec& D,
    const arma::vec& y,
    const arma::vec& w,
    const arma::vec& fitted,
    const arma::vec& dhat,
    const double& lambda) {

  // arma::vec score = -2 * (mult_diag(U.t(), w)  * (y - fitted) + lambda * (accu(w) / y.n_elem) * (dhat / D));
  arma::vec score = -2 * (mult_diag(U.t(), w)  * (y - fitted) + lambda * (dhat / D));
  
  return score;
}

// Hessian used for sandwich estimator for KRLS
// [[Rcpp::export]]
arma::mat krls_hess_trunc(const arma::mat& U,
                          const arma::vec& D,
                          const arma::vec& w,
                           const double& lambda) {

  arma::mat hess = 2 * (mult_diag(U.t(), w) * U + arma::diagmat(lambda * (1 / D)));

  return hess;

}

// [[Rcpp::export]]
double krlogit_fn_trunc(const arma::vec& par,
                        const arma::mat& U,
                        const arma::vec& D,
                        const arma::vec& y,
                        const arma::vec& w,
                        const double& lambda) {

  arma::vec coef = par.subvec(0, par.n_elem - 2);
  double beta0 = par(par.n_elem-1);
  arma::mat Ud = U * coef;

  double ret = accu(w % (y % log(1 + exp(-(beta0 + Ud))) +
                        (1 - y) % log(1 + exp(beta0 + Ud))) +
                        lambda / y.n_elem * arma::as_scalar(coef.t()  * ((1.0/D) % coef)));

  return ret;

}

// [[Rcpp::export]]
arma::vec krlogit_gr_trunc(const arma::vec& par,
                        const arma::mat& U,
                        const arma::vec& D,
                        const arma::vec& y,
                        const arma::vec& w,
                        const double& lambda) {

  arma::vec coef = par.subvec(0, par.n_elem - 2);
  double beta0 = par(par.n_elem-1);
  arma::vec resid = y - (1 / (1 + exp(-U * coef - beta0)));

  arma::vec ret(par.n_elem);

  ret.subvec(0, par.n_elem - 2) = -U.t() * (w % resid) + 2 * (lambda / D) % coef;
  ret(par.n_elem - 1) = -accu(w % resid);

  return ret;
}

// Generates (exp(-Kc - b)) / (1+exp(-Kc-b))^2 used for the krlogit hessian
// and for for predicted value SEs and first difference SEs
// Should work with either truncated data or full data
// [[Rcpp::export]]
arma::vec partial_logit(const arma::mat& K,
                        const arma::vec& coef,
                        const double& beta0) {
  arma::vec ret = exp(-K * coef - beta0) / pow((1 + exp(-K * coef - beta0)), 2);

  return ret;
}

// [[Rcpp::export]]
arma::mat krlogit_hess_trunc(const arma::vec& par,
                           const arma::mat& U,
                           const arma::vec& D,
                           const arma::vec& y,
                           const arma::vec& w,
                           const double& lambda) {

  arma::vec coef = par.subvec(0, par.n_elem - 2);
  double beta0 = par(par.n_elem-1);
  arma::vec meat = w % partial_logit(U, coef, beta0);

  arma::mat ret(par.n_elem, par.n_elem);

  arma::mat dcdc = mult_diag(U.t(), meat) * U + diagmat(2 * (lambda / D));

  arma::vec dcdb = U.t() * meat;

  double dbdb = accu(meat);

  ret.submat(0, 0, coef.n_elem - 1, coef.n_elem - 1) = dcdc;
  ret.submat(coef.n_elem, 0, coef.n_elem, coef.n_elem - 1) = dcdb.t();
  ret.submat(0, coef.n_elem, coef.n_elem - 1, coef.n_elem) = dcdb;
  ret(coef.n_elem, coef.n_elem) = dbdb;

  return ret;
}