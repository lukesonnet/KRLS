// No longer used because always use eigen decomposition
//' @export
// [[Rcpp::export]]
arma::vec krls_gr(
               const arma::mat& K,
               const arma::vec& y,
               const arma::vec& fitted,
               const arma::vec& fittedFull,
               const double& lambda) {
  
  arma::vec score = -2 * K * (y - fitted) + 2 * lambda * fittedFull; 
  
  return score;
  
}

//' @export
// [[Rcpp::export]]
arma::mat krls_hess_inv(const arma::mat& K,
                        const double& lambda) {
  
  arma::mat I(K.n_rows, K.n_cols, arma::fill::eye);
  arma::mat hess = 2 * K * (K + I * lambda);
  
  return arma::inv_sympd(hess);
  
}