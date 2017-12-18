## This file contains helper functions
## Functions:
##            chunk
## See in CPP:
##            mult_diag

## Function that splits a vector in to n chunks
## From Stack Overflow answer by mathheadinclouds
## link: https://stackoverflow.com/a/16275428
#' @export
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


## The logistic function that takes values for coeff, b0, and a K or Ktilde
#' @importFrom stats plogis
#' @export
logistic <- function(K, coeff, beta0) {
  plogis(beta0 + K %*% coeff)
}
