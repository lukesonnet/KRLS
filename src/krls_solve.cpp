#include <RcppArmadillo.h>
#include "krls_helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List solve_for_d_ls(const arma::vec& y,
                          const arma::mat& U,
                          const arma::vec& D,
                          const double& lambda) {

  arma::vec Ginv = 1 / (1 + lambda / D);

  arma::vec dhat = Ginv % (U.t() * y);
  // This is the same as above
  //arma::vec tempLoss = (y - U * coeffs) / diagvec(arma::eye(y.n_elem, y.n_elem) - mult_diag(U, Ginv) * U.t());

  arma::vec tempLoss = (y - U * mult_diag(U, Ginv).t() * y) / diagvec(arma::eye(y.n_elem, y.n_elem) - U * mult_diag(U, Ginv).t());
  double Le = as_scalar(tempLoss.t() * tempLoss);

  return Rcpp::List::create(Rcpp::Named("dhat") = dhat,
                            Rcpp::Named("Le") = Le);
}

// [[Rcpp::export]]
Rcpp::List solve_for_d_ls_w(const arma::vec& y,
                          const arma::mat& U,
                          const arma::vec& D,
                          const arma::vec& w,
                          const double& lambda) {

  arma::mat Uw = mult_diag(U.t(), w);
  arma::mat Ginv = arma::inv_sympd(Uw * U + diagmat(lambda / D));

  arma::vec dhat = Ginv * (Uw * y);
  // This is the same as above
  //arma::vec tempLoss = (y - U * coeffs) / diagvec(arma::eye(y.n_elem, y.n_elem) - mult_diag(U, Ginv) * U.t());
  arma::mat UGinvUw = U * Ginv * Uw;

  arma::vec tempLoss = (y - UGinvUw * y) / diagvec(arma::eye(y.n_elem, y.n_elem) - UGinvUw);
  double Le = as_scalar(tempLoss.t() * tempLoss);

  return Rcpp::List::create(Rcpp::Named("dhat") = dhat,
                            Rcpp::Named("Le") = Le);
}