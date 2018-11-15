#include <RcppArmadillo.h>
using namespace Rcpp;

// Compute pointwise marginal effects and the var avg pwmfx for one variable
// [[Rcpp::export]]
Rcpp::List pwmfx(const arma::mat& k,
                const arma::vec& x,
                const arma::vec& coefhat,
                const Rcpp::Nullable<Rcpp::NumericMatrix>& vcovc_mat,
                const arma::vec& p,
                const arma::vec& p2,
                const double& b)
{

  arma::mat vcovc;
  if (vcovc_mat.isNotNull()) {
    vcovc = Rcpp::as<arma::mat>(vcovc_mat);
  }
  double n = x.n_elem;
  arma::vec deriv(n);
  arma::mat var_deriv(n, n);
  double var_avg_deriv;
  arma::mat distmat(n, n);
  double val = 0.0;

  for (unsigned i = 0; i < n; ++i) {
    val = 0.0;
    for (unsigned i2 = 0; i2 < n; ++i2) {
      distmat(i, i2) = x(i) - x(i2);

      val += coefhat(i2) * k(i, i2) * distmat(i, i2);
    }

    deriv(i) = - (p(i) / b)  * val;
  }

  if (vcovc_mat.isNotNull()) {
    arma::mat distk = k % distmat;
    var_avg_deriv = 1 / std::pow(b * n, 2) * accu(p2.t() * distk.t() * vcovc * distk);
    
    var_deriv = 4 / std::pow(b, 2) * distk.t() * vcovc * distk;
  }

  return Rcpp::List::create(
    Rcpp::Named("deriv") = deriv,
    Rcpp::Named("var_deriv") = var_deriv,
    Rcpp::Named("var_avg_deriv") = var_avg_deriv
  );
}
