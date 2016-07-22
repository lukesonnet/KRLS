## This file contains helper functions
## Functions:
##            chunk
## See in CPP:
##            mult_diag

## Function that splits a vector in to n chunks
#' @export
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

