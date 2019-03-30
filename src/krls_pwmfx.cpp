#include <RcppArmadillo.h>
using namespace Rcpp;

// const ,
// Compute pointwise marginal effects and the var avg pwmfx for one variable
// [[Rcpp::export]]
Rcpp::List pwmfx(const arma::mat& k,
                 const arma::vec& x,
                 const arma::vec& coefhat,
                 const arma::mat& vcovc,
                 const arma::vec& p,
                 const double& b,
                 const bool& computevarderiv)
{
  
  unsigned n = x.n_elem;
  arma::vec deriv(n);
  arma::mat distmat(n, n);
  double val = 0.0;
  
  for (unsigned i = 0; i < n; ++i) {
    val = 0.0;
    for (unsigned i2 = 0; i2 < n; ++i2) {
      distmat(i, i2) = x(i) - x(i2);
      
      val += coefhat(i2) * k(i, i2) * distmat(i, i2);
    }
    
    deriv(i) = - (2 * p(i) / b)  * val;
  }
  
  // We use a matrix of size 1 of value 0 when we do not want
  // to compute the variance of the pwmfx
  // We do this instead of allowing vcovc to be null for memory reasons
  if (vcovc.n_elem == 1 && vcovc(0, 0) == 0) {
    return Rcpp::List::create(
      Rcpp::Named("deriv") = deriv,
      Rcpp::Named("var_deriv") = Rcpp::NumericVector::create(NA_REAL),
      Rcpp::Named("var_avg_deriv") = NA_REAL
    );
  } else {
    arma::mat distk = k % distmat;
    
    arma::mat var_mat = 4 / std::pow(b * n, 2) * p.t() * distk * vcovc * distk.t() * p;
    double var_avg_deriv = var_mat(0, 0);
    arma::mat var_deriv;
    if (computevarderiv) {
      var_deriv = 4 / std::pow(b, 2) * distk.t() * vcovc * distk;
    } else {
      var_deriv = 0.0;
    }
    
    return Rcpp::List::create(
      Rcpp::Named("deriv") = deriv,
      Rcpp::Named("var_deriv") = var_deriv,
      Rcpp::Named("var_avg_deriv") = var_avg_deriv
    );
  }
}
