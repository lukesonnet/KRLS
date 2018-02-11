#include <RcppArmadillo.h>
using namespace Rcpp;

// ----------
// GAUSSIAN KERNELS
// ----------

// Euclidean distance function

// [[Rcpp::export]]
double euc_dist(const arma::rowvec& x1, const arma::rowvec& x2) {
  double out = 0.0;
  unsigned n = x1.n_elem;

  for (unsigned i = 0; i < n; ++i) {
    out += pow(x1(i) - x2(i), 2);
  }
  return sqrt(out);
}

// Gaussian kernel between two row vectors

// [[Rcpp::export]]
double kern_gauss_1d(const arma::rowvec& x1, const arma::rowvec& x2, const double& b)
{
  return exp(-pow(euc_dist(x1, x2), 2) / (b));
}

// Gaussian kernel matrix for a matrix with itself
// b is sigma, or 2P by default
// [[Rcpp::export]]
arma::mat kern_gauss(const arma::mat& x, const double& b)
{
  unsigned n = x.n_rows;
  double val;
  // Diagonal will remain ones
  arma::mat out(n, n, arma::fill::ones);

  for (unsigned i = 0; i < n; ++i) {

    for (unsigned j = i + 1; j < n; ++j) {
      val = kern_gauss_1d(x.row(i), x.row(j), b);
      out(i, j) = val;
      out(j, i) = val;
    }

  }
  return out;
}

// Kernel matrix for distance between two matrices
// [[Rcpp::export]]
arma::mat new_gauss_kern(const arma::mat& newx, const arma::mat& oldx, const double& b) {
  unsigned n1 = newx.n_rows;
  unsigned n2 = oldx.n_rows;
  arma::mat out(n1, n2);

  for (unsigned i = 0; i < n1; ++i) {
    for (unsigned j = 0; j < n2; ++j) {
      out(i, j) = kern_gauss_1d(newx.row(i), oldx.row(j), b);
    }

  }
  return out;
}