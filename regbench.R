# function sum, product and composition -----------------------------------

f_sum <- function(g, alpha= rep(1, length(g))){
    g <- g
    alpha <- alpha
    function(x) sum(alpha*sapply(g, function(gi)gi(x)))
}

f_product <- function(g, alpha = rep(1, length(g))){
    g <- g
    alpha <- alpha
    function(x) prod(alpha*sapply(g, function(gi) gi(x)))
}

f_comp <- function(g){
    g <- g
    function(x){
        y <- g[[1]](x)
        g <- g[-1]
        for (i in seq_along(g)) {
            y <- g[[i]](y)
        }
        return(y)
    }
}


# base functions ----------------------------------------------------------

f_gauss <- function(mu, Vl, idx){
    mu <- mu
    Vl <- Vl
    idx <- idx
    function(x) exp(-0.5 * t(x[idx] - mu) %*% Vl %*% (x[idx] - mu))
}




# composite functions -----------------------------------------------------

# helpers

V <- function(n, min_eigen = 0.1, max_eigen= 2){
    sqrt_d  <- runif(n, min_eigen, max_eigen)
    Ul      <- pracma::randortho(n, type = "orthonormal")
    Vl      <- Ul %*% (diag(n)*sqrt_d^2) %*% t(Ul)
    return(Vl)
}

f_fried <- function(p = 10, ng = 20, base = 1.5, rate = 2){
    p <- p
    ng <- ng
    base <- base
    rate <- rate
    n_int <- pmin(p, floor(base + rexp(ng, rate = rate)))
    g <- vector("list", length= length(ng))
    for(i in 1:ng){
        idx <- sample(1:p, n_int[i])
        mu <- rnorm(n_int[i])
        Vl <- V(n_int[i])
        g[[i]] <- f_gauss(mu = mu, Vl = Vl, idx = idx)
    }
    alpha <- runif(ng, -1, 1)
    fx <- f_sum(g = g, alpha = alpha)
    f <- function(x) {
        x <- x
        if (!is.matrix(x)) stop("x must be a matrix")
        apply(x, 1, fx)
    }
    return(f)
}


# simulators --------------------------------------------------------------

##'@export
##'@import MASS pracma
sim_friedman <- function(p = 10, ng = 2*p, sn_ratio = 1){
    p <- p
    ng <- ng
    sn_ratio <- sn_ratio
    x <- MASS::mvrnorm(n = 1000, mu = rep(0, p), Sigma = diag(p))
    f <- f_fried(p = p, ng = ng)
    fx <- f(x)
    sde <- sd(fx)*sn_ratio
    function(n = 100){
        x <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
        e <- rnorm(n, sd = sde)
        fx <- f(x)
        y <- fx + e
        res <- list(fx = f,
                    data = data.frame(x,y,fx,e))
        return(res)
    }
}
