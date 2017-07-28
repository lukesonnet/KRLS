## This file contains helper functions
## Functions:
##            chunk
## See in CPP:
##            mult_diag

## Function that splits a vector in to n chunks
#' @export
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


## The logistic function that takes values for coeff, b0, and a K or Ktilde
#' @export
logistic <- function(K, coeff, beta0) {
  
  yhat <- 1 / (1 + exp(-(beta0+K%*%coeff)))
  
  return(yhat)
}
