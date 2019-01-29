#include <RcppArmadillo.h>
using namespace Rcpp;

// Compute pointwise marginal effects and the var avg pwmfx for one variable
// [[Rcpp::export]]
arma::vec pwmfx(const arma::mat& k,
                const arma::vec& x,
                const arma::vec& coefhat,
                const Rcpp::Nullable<Rcpp::NumericMatrix>& vcovc_mat,
                const arma::vec& p,
                const double& b)
{

  arma::mat vcovc;
  if (vcovc_mat.isNotNull()) {
    vcovc = Rcpp::as<arma::mat>(vcovc_mat);
  }
  double n = x.n_elem;
  arma::vec out(n + 1);
  arma::mat distmat(n, n);
  double val = 0.0;

  for (unsigned i = 0; i < n; ++i) {
    val = 0.0;
    for (unsigned i2 = 0; i2 < n; ++i2) {
      distmat(i, i2) = x(i) - x(i2);

      val += coefhat(i2) * k(i, i2) * distmat(i, i2);
    }

    out(i) = - (2 * p(i) / b)  * val;
  }

  if (vcovc_mat.isNotNull()) {
    arma::mat distk = k % distmat;
    arma::mat var1 = 4 / pow(b * n, 2) * p.t() * distk * vcovc * distk.t() * p;
    out(n) = var1(0, 0);
  }

  return out;
}
