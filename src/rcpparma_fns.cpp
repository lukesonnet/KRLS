// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//----------
// Helper functions
//----------

//' @export
// [[Rcpp::export]]
arma::mat mult_diag(const arma::mat& x, const arma::vec& d) {
  
  arma::mat out(x.n_rows, x.n_cols);
  for (int j = 0; j < x.n_cols; ++j) {
    out.col(j) = x.col(j) * d(j);
  }

  return out;
}


// Eigen decomposition, unused

//' @export
// [[Rcpp::export]]
Rcpp::List eigsym(const arma::mat& x) {
  arma::vec eigval;
  arma::mat eigvec;
  
  eig_sym(eigval, eigvec, x);
  
  return Rcpp::List::create(Rcpp::Named("eigvec") = eigvec,
                            Rcpp::Named("eigval") = eigval);
}

//' @export
// [[Rcpp::export]]
arma::vec krls_gr_trunc(
    const arma::mat& U,
    const arma::vec& D,
    const arma::vec& y,
    const arma::vec& fitted,
    const arma::vec& dhat,
    const double& lambda) {
  
  arma::vec score = -2 * (U.t() * (y - fitted) + lambda * (dhat / D)); 

  return score;
}

// Hessian used for sandwich estimator for KRLS
//' @export
// [[Rcpp::export]]
arma::mat krls_hess_trunc_inv(const arma::mat& U,
                          const arma::vec& D,
                           const double& lambda) {
  
  arma::mat hess = 2 * (U.t() * U + arma::diagmat(lambda * (1 / D)));
  
  return arma::inv_sympd(hess);
  
}

//' @export
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
                        lambda * arma::as_scalar(coef.t()  * ((1.0/D) % coef)));
  
  return ret;
  
}

//' @export
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
  
  ret.subvec(0, par.n_elem - 2) = -U.t() * (w % resid) + 2 * y.n_elem * (lambda / D) % coef;
  ret(par.n_elem - 1) = -accu(w % resid);
  
  return ret;
}

// Generates (exp(-Kc - b)) / (1+exp(-Kc-b))^2 used for the krlogit hessian
// and for for predicted value SEs and first difference SEs
// Should work with either truncated data or full data
//' @export
// [[Rcpp::export]]
arma::vec partial_logit(const arma::mat& K,
                        const arma::vec& coef,
                        const double& beta0) {
  arma::vec ret = exp(-K * coef - beta0) / pow((1 + exp(-K * coef - beta0)), 2);
  
  return ret;
}

//' @export
// [[Rcpp::export]]
arma::mat krlogit_hess_trunc_inv(const arma::vec& par,
                           const arma::mat& U,
                           const arma::vec& D,
                           const arma::vec& y,
                           const arma::vec& w,
                           const double& lambda) {
  
  arma::vec coef = par.subvec(0, par.n_elem - 2);
  double beta0 = par(par.n_elem-1);
  arma::vec meat = w % partial_logit(U, coef, beta0);

  arma::mat ret(par.n_elem, par.n_elem);

  arma::mat dcdc = mult_diag(U.t(), meat) * U + diagmat(2 * y.n_elem * (lambda / D));
  
  arma::vec dcdb = U.t() * meat;
  
  double dbdb = accu(meat);
  
  ret.submat(0, 0, coef.n_elem - 1, coef.n_elem - 1) = dcdc;
  ret.submat(coef.n_elem, 0, coef.n_elem, coef.n_elem - 1) = dcdb.t();
  ret.submat(0, coef.n_elem, coef.n_elem - 1, coef.n_elem) = dcdb;
  ret(coef.n_elem, coef.n_elem) = dbdb;
  
  return arma::inv_sympd(ret);
}

// ----------
// GAUSSIAN KERNELS
// ----------

// Euclidean distance function

//' @export
// [[Rcpp::export]]
double euc_dist(const arma::rowvec& x1, const arma::rowvec& x2) {
  double out = 0.0;
  int n = x1.n_elem;
  
  for (int i = 0; i < n; ++i) {
    out += pow(x1(i) - x2(i), 2);
  }
  return sqrt(out);
}

// Gaussian kernel between two row vectors

//' @export
// [[Rcpp::export]]
double kern_gauss_1d(const arma::rowvec& x1, const arma::rowvec& x2, const double& b)
{
  return exp(-pow(euc_dist(x1, x2), 2) / (b));
}

// Gaussian kernel matrix for a matrix with itself
// b is sigma, or 2P by default

//' @export
// [[Rcpp::export]]
arma::mat kern_gauss(const arma::mat& x, const double& b)
{
  int n = x.n_rows;
  double val;
  // Diagonal will remain ones
  arma::mat out(n, n, arma::fill::ones);
  
  for (int i = 0; i < n; ++i) {

    for (int j = i + 1; j < n; ++j) {
      val = kern_gauss_1d(x.row(i), x.row(j), b);
      out(i, j) = val;
      out(j, i) = val;
    }
    
  }
  return out;
}

// Kernel matrix for distance between two matrices

//' @export
// [[Rcpp::export]]
arma::mat new_gauss_kern(const arma::mat& newx, const arma::mat& oldx, const double& b) {
  int n1 = newx.n_rows;
  int n2 = oldx.n_rows;
  arma::mat out(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      out(i, j) = kern_gauss_1d(newx.row(i), oldx.row(j), b);
    }
    
  }
  return out;
}



//' @export
// [[Rcpp::export]]
double lambda_search(const double& tol,
                     const double& l,
                     const double& u,
                     const arma::vec& y,
                     const arma::vec& eigvals,
                     const arma::mat& eigvecs,
                     const double& eigtrunc) {
  return 42;
}

//' @export
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

// Compute pointwise marginal effects and the var avg pwmfx
//' @export
// [[Rcpp::export]]
arma::mat pwmfx(const arma::mat& k,
                const arma::mat& x,
                const arma::vec& coefhat,
                const arma::mat& vcovc,
                const arma::vec& p,
                const double& b)
{
  //for now until we can easily return list
  double n = x.n_rows;
  arma::mat out(n + 1, x.n_cols);
  arma::mat distmat(n, n);
  arma::mat distk(n, n);
  arma::vec p2 = pow(p, 2);
  double val;
  
  for (int j = 0; j < x.n_cols; ++j) {
    for (int i = 0; i < n; ++i) {
      val = 0;
      for (int i2 = 0; i2 < n; ++i2) {
        distmat(i, i2) = x(i, j) - x(i2, j);
        val += coefhat(i2) * k(i, i2) * distmat(i, i2);
      }
      
      out(i, j) = - (p(i) / b)  * val;
    }
    distk = k % distmat;
    out(n, j) = 1 / pow(b * n, 2) * accu(p2.t() * distk.t() * vcovc * distk);
    
  }
  
  return out;
}
