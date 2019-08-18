# KRLS with sketching

ske_KRLS <- function(X, y, b, m, 
                     lambdas =  c(10^(seq(-6, 2, length.out = 20))),
                     folds = 5,
                     n_av, control){
    n <- nrow(X)
    
    # parameters in sketching
    K <- switch(control$whichkernel,
                gaussian=kern_gauss(X, b),
                linear=tcrossprod(X),
                gausslinear=.5*(kern_gauss(X, b)+tcrossprod(X)),
                poly2=(tcrossprod(X)+1)^2,
                poly3=(tcrossprod(X)+1)^3,
                poly4=(tcrossprod(X)+1)^4,
                stop("No valid Kernel specified") )
    S <- matrix(rnorm(m*n), nrow = n, ncol = m)/sqrt(m)
    R <- K %*% S
    D <- t(S) %*% K %*% S
    
    
    # construct R matrix by parallel
    R <- apply(X, 1, function(x){
      onerow <-  switch(control$whichkernel,
                        gaussian=new_gauss_kern(x, X, b),
                        linear=tcrossprod(x, X),
                        gausslinear=.5*(new_gauss_kern(x, X, b)+tcrossprod(x, X)),
                        poly2=(tcrossprod(x, X)+1)^2,
                        poly3=(tcrossprod(x, X)+1)^3,
                        poly4=(tcrossprod(x, X)+1)^4,
                        stop("No valid Kernel specified") )
      apply(S, 2, function(x) crossprod(onerow, x))
    })
    
    # assess the quality of the selected columns
    # calculating the average R square
    K_y <-switch(control$whichkernel,
                 gaussian=new_gauss_kern(X, X[sample(n, n_av), ], b),
                 linear=tcrossprod(X, X[sample(n, n_av), ]),
                 gausslinear=.5*(new_gauss_kern(X, X[sample(n, n_av), ], b)+tcrossprod(X, X[sample(n, n_av), ])),
                 poly2=(tcrossprod(X, X[sample(n, n_av), ])+1)^2,
                 poly3=(tcrossprod(X, X[sample(n, n_av), ])+1)^3,
                 poly4=(tcrossprod(X, X[sample(n, n_av), ])+1)^4,
                 stop("No valid Kernel specified") )
    Rsq <- ave_Rsq(R, K_y)
    
    # cross validation
    lables <- sample(folds, n, replace = T)
    cv_MSE <- rep(0, length(lambdas))
    
    for (i in 1:length(lambdas)){
        lambda <- lambdas[i]
        MSE_vali <- rep(0, folds)
        for (v in 1:folds){
            # training indices in cross validation
            cv_idx <- which(lables != v)
            
            # sketching
            dh_temp <- train_krr(y[cv_idx], R[cv_idx, ], D, lambda)$dh
            MSE_vali[v] <- test_krr(dh_temp, y[-cv_idx], R[-cv_idx, ])$test_MSE
        }
        
        #cross validation MSE nystrom
        cv_MSE[i] <- mean(MSE_vali)
    }
    
    lambda_opt <- lambdas[which.min(cv_MSE)]
    train_result <- train_krr(y, R, D, lambda_opt)
    
    return(list(coeffs = train_result$dh,
                fitted = train_result$yfitted,
                I = S,      
                Rsq = Rsq, 
                lambda = lambda_opt))
}
