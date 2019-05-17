#include <RcppArmadillo.h>
using namespace Rcpp;

// Compute pointwise marginal effects and the var avg pwmfx for one variable
// [[Rcpp::export]]
Rcpp::List pwmfx(const arma::mat& Xstar,
                 const arma::mat& X,
                 const unsigned& wrt_column,
                 const arma::mat& K,
                 const arma::vec& coefhat,
                 const arma::mat& Vcovc,
                 const arma::vec& p,
                 const double& b,
                 const bool& computevarderiv,
                 const bool& computederiv2)
{
  unsigned J = Xstar.n_rows;
  unsigned N = X.n_rows;
  arma::vec deriv(J);
  arma::mat distmat(N, J);
  double val = 0.0;
  
  for (unsigned j = 0; j < J; ++j) {
    val = 0.0;
    for (unsigned i = 0; i < N; ++i) {
      distmat(i, j) = X(i, wrt_column) - Xstar(j, wrt_column);
      
      val += coefhat(i) * K(i, j) * distmat(i, j);
    }
    
    deriv(j) = (2 * p(j) / b) * val;
  }
  
  // We use a matrix of size 1 of value 0 when we do not want
  // to compute the variance of the pwmfx
  // We do this instead of allowing vcovc to be null for memory reasons
  if (Vcovc.n_elem == 1 && Vcovc(0, 0) == 0) {
    return Rcpp::List::create(
      Rcpp::Named("deriv") = deriv,
      Rcpp::Named("var_deriv") = Rcpp::NumericVector::create(NA_REAL),
      Rcpp::Named("var_avg_deriv") = NA_REAL
    );
  } else {
    arma::mat distK = K % distmat;
    
    arma::mat var_mat = 4 / std::pow(b * N, 2) * p.t() * distK * Vcovc * distK.t() * p;
    double var_avg_deriv = var_mat(0, 0);
    arma::mat var_deriv;
    if (computevarderiv) {
      var_deriv = 4 / std::pow(b, 2) * distK.t() * Vcovc * distK;
    } else {
      var_deriv = 0.0;
    }
    
    arma::vec avg_second_deriv(X.n_cols);
    if (computederiv2) {
      double val2;
      for (unsigned d = 0; d < X.n_cols; ++d) {
        val2 = 0.0;
        for (unsigned j = 0; j < N; ++j) {
          for (unsigned i = 0; i < N; ++i) {
            val2 += coefhat(i) * K(i, j) * distmat(i, j) * (X(i, d) - X(j, d));
          }
        }
        
        avg_second_deriv(d) = (4 / (std::pow(b, 2))) * val2 / N;
      }
    }
    
    return Rcpp::List::create(
      Rcpp::Named("deriv") = deriv,
      Rcpp::Named("var_deriv") = var_deriv,
      Rcpp::Named("var_avg_deriv") = var_avg_deriv,
      Rcpp::Named("avg_second_deriv") = avg_second_deriv
    );
  }
}