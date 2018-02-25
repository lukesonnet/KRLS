# KRLS 1.1.0 (Current Github version)

* Temporarily rename package KRLS2
* Use C++ for speed improvements throughout
* Add ability to use logistic loss instead of least squares loss for binary outcomes (paper forthcoming)
* Enable truncation to improve speed with large datsets
* Change `sigma` to `b` for consistency
* Change default value of `b` from `p` to `2*p` where `p` is the number of covariates
* Change first differences algorithm for speed, including correcting variance of average first differences, which was too large by a factor of 2

# KRLS 1.0.0 (Current CRAN version)

* Main functionality as described in [Hainmueller and Hazlett (2014)](https://doi.org/10.1093/pan/mpt019) [[ungated]](https://web.stanford.edu/~jhain/Paper/PA2014a.pdf)
