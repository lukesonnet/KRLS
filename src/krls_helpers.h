#ifndef KRLS_HELPER_H
#define KRLS_HELPER_H

#include <RcppArmadillo.h>
arma::mat mult_diag(const arma::mat& x, const arma::vec& d);
double trace_mat(const arma::mat& x);

#endif