#include <RcppArmadillo.h>
#include "krls_helpers.h"
using namespace Rcpp;

//----------
// Helper functions
//----------

// [[Rcpp::export]]
arma::mat mult_diag(const arma::mat& x, const arma::vec& d) {
  
  arma::mat out(x.n_rows, x.n_cols);
  for (unsigned j = 0; j < x.n_cols; ++j) {
    out.col(j) = x.col(j) * d(j);
  }
  
  return out;
}

// [[Rcpp::export]]
double trace_mat(const arma::mat& x) {
  return arma::trace(x);
}