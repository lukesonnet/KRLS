## This file contains functions that build kernel matrices and truncate them
## Functions:
##            generateK
##            Ktrunc
##            newKernel (although new_gauss_kern in cpp is used in some places)
## Also see in the cpp:
##

## This function generates the kernel matrix and truncated kernel matrix object
#' @export
generateK <- function(X,
                      b,
                      control) {
  K <- NULL
  if(control$whichkernel=="gaussian"){ K <- kern_gauss(X, b)}
  if(control$whichkernel=="linear"){ K <- tcrossprod(X) }
  if(control$whichkernel=="poly2"){ K <- (tcrossprod(X)+1)^2 }
  if(control$whichkernel=="poly3"){ K <- (tcrossprod(X)+1)^3 }
  if(control$whichkernel=="poly4"){ K <- (tcrossprod(X)+1)^4 }
  if(is.null(K)){ stop("No valid Kernel specified") }

  if(control$truncate) {
    truncDat <- Ktrunc(K = K, b=b, lastkeeper=control$lastkeeper, epsilon=control$epsilon,
                       quiet = control$quiet)
    U <- truncDat$Utrunc
    D <- truncDat$eigvals
  } else {
    eigobj <- eigen(K)
    U <- eigobj$vectors
    D <- eigobj$values
  }
  
  lastkeeper <- NULL
  lastkeeper <- if(control$truncate) ncol(U)

  return(list(K = K,
              U = U,
              D = D,
              lastkeeper = lastkeeper))
}


## Function that returns truncated versions of the data if given
## todo: throws warning whenever n < 500 because it uses eigen, should we suppress?
#' @export
Ktrunc <- function(X=NULL, K=NULL, b=NULL, epsilon=NULL, lastkeeper=NULL, quiet = TRUE){
  if(is.null(K)){
    X.sd <- apply(X, 2, sd)
    X.mu <- colMeans(X)
    K <- kern_gauss(scale(X), ifelse(is.null(b), ncol(X), b))
  }

  if(is.null(lastkeeper)){
    #denoms <- c(10, 6, 3, 1) # we could also let people set this to speed up the algorithm
    #CJH: even at 10, with N=5000 it is too big. Let's try this: 
    if (nrow(K) <= 500) numvectorss=nrow(K) else numvectorss=c(50, 100, 200, 300,500, min(c(nrow(K), 1000)))
    enoughvar=FALSE
    j=1
    while (enoughvar==FALSE){
      numvectors=numvectorss[j]
      if(!quiet) print(paste("trying",numvectors,"vectors"))
      eigobj <- NULL
      eigobj <- suppressWarnings({eigs_sym(K, numvectors, which="LM")})
      
      #for now, letting it throw an error in certain failure cases.
      totalvar=sum(eigobj$values)/nrow(K)
      
      if (totalvar>=(1-epsilon)){
        lastkeeper = min(which(cumsum(eigobj$values)/nrow(K) > (1-epsilon)))

        if(!quiet) print(paste("lastkeeper=",lastkeeper))
        enoughvar=TRUE
        
      }
      j=j+1
      # Will eventually add warning here for var < 0.99 at numvector= 300
      # Allow to proceed w/out truncation and SEs or to try full eigen
    } #end while gotenough==FALSE      
  } else { # if !is.null(lastkeeper)
    ## Suppress warning about using all eigenvalues
    eigobj <- suppressWarnings({eigs_sym(K, lastkeeper, which="LM")})
  }
  
  return(
    list(Utrunc=eigobj$vectors[, 1:lastkeeper, drop = F],
         eigvals=eigobj$values[1:lastkeeper, drop = F])
  )
}


## Produce a new Kernel from new data and data that was used to fit model for
## prediction
## Parameters:
##   'X' - the original data matrix, UNSCALED
##   'newData' - the new data matrix (with same features), unscaled
##   'whichkernel' - which kernel was used on the original data
## Values:
##   'newK' - The new Kernel to be used for prediction
#' @export
newKernel <- function(X, newData, whichkernel = "gaussian", b = NULL) {
  
  # Get values of oldData to scale newData
  Xmeans <- colMeans(X)
  Xsd <- apply(X, 2, sd)
  X <- scale(X, center = Xmeans, scale = Xsd)
  
  # scale new data by means and sd of old data
  newData <- scale(newData, center = Xmeans, scale = Xsd)      
  
  # predict based on new kernel matrix
  nn <- nrow(newData)
  
  ## Compute kernel matrix
  K <- NULL
  if(whichkernel=="gaussian"){ 
    if(is.null(b)) {b <- 2 * ncol(X)}
    newK <- new_gauss_kern(newData, X, b)
  }
  
  if(whichkernel=="linear"){ K <- tcrossprod(rbind(newData, X)) }
  if(whichkernel=="poly2"){ K <- (tcrossprod(rbind(newData, X))+1)^2 }
  if(whichkernel=="poly3"){ K <- (tcrossprod(rbind(newData, X))+1)^3 }
  if(whichkernel=="poly4"){ K <- (tcrossprod(rbind(newData, X))+1)^4 }
  if(is.null(K) & is.null(newK)){ stop("No valid Kernel specified") }
  
  if(whichkernel != "gaussian"){
    newK <- matrix(K[1:nn, (nn+1):(nn+nrow(X))],
                   nrow=nrow(newData),
                   byrow=FALSE)
  }
  
  return(newK)
}