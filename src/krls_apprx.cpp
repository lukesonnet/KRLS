#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 
using namespace arma;

// training KRR 
// [[Rcpp::export]]
List train_krr(const arma::vec& y_train,
               const arma::mat& R_train,
               const arma::mat& D,
               const double& lambda) {
    // check the condition number of R^TR first
    arma::mat M = R_train.t() * R_train + lambda * D;
    arma::mat dhat;
    double M_cond = cond(M);
    
    if (M_cond > 1e6) {
        cout << "Inverting is ill-conditioned! Added perturbation to the matrix. Try smaller m and/or larger b." << endl;
        arma::mat fixer = 1e-6 * eye( size(D) );
        dhat =  solve( M + fixer,  R_train.t() * y_train );
    }else {
        dhat =  solve( M,  R_train.t() * y_train );
    }
    
    arma::mat yfitted = R_train * dhat;
    //double train_MSE = mean(square(y_train-yh_train));
    
    return List::create(Named("dh") = dhat,
                        Named("M_cond") = M_cond,
                        Named("yfitted") = yfitted);
}

// testing KRR 
// [[Rcpp::export]]
List test_krr(const arma::vec& dhat,
              const arma::vec& y_test,
              const arma::mat& R_test) {
    mat yh_test = R_test * dhat;
    double test_MSE =  mean(square(y_test - yh_test));
    return List::create(Named("yh_test") = yh_test,
                        Named("test_MSE") = test_MSE);
}

// get the average Rsq
// [[Rcpp::export]]
double ave_Rsq(const arma::mat& X,
            const arma::mat& y) {
    double Rsq = 0;
    double mss, rss;
    unsigned m = y.n_cols;
    
    mat yh = X * solve(X.t() * X, X.t()*y);
    mat residuals = y - yh;
    for (unsigned i = 0; i < m; i++) {
        mss = sum(square(yh.col(i)));
        rss = sum(square(residuals.col(i)));
        Rsq += mss/(mss + rss);
    }
    return Rsq/m;
}


