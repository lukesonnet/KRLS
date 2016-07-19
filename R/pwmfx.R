# This file contains functions that generate the pwmfx for krls

#' @export
pwmfxR <- function() {
  rows <- cbind(rep(1:n, each = n), 1:n)
  distances <- X[rows[,1],] - X[ rows[,2],] 
  
  for(k in 1:d){
    print(paste("Computing derivatives, variable", k))
    if(d==1){
      distk <-  matrix(distances,n,n,byrow=TRUE)
    } else {
      distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
    }
    L <-  distk*Kdat$K  #Kfull rather than K here as truncation handled through coefs
    derivatives[,k] <- tau*(-1/sigma)*(L%*%coefhat)
    if(truncate==FALSE) {
      var.avgderivatives = NULL
    } else {
      var.avgderivatives[1,k] <- 1/(sigma * n)^2 * sum(crossprod(tau^2, crossprod(L,vcov.c%*%L)))
    }
  }
}