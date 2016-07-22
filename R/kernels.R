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
                      sigma,
                      control) {
  K <- NULL
  if(control$whichkernel=="gaussian"){ K <- kern_gauss(X, sigma)}
  if(control$whichkernel=="linear"){ K <- tcrossprod(X) }
  if(control$whichkernel=="poly2"){ K <- (tcrossprod(X)+1)^2 }
  if(control$whichkernel=="poly3"){ K <- (tcrossprod(X)+1)^3 }
  if(control$whichkernel=="poly4"){ K <- (tcrossprod(X)+1)^4 }
  if(is.null(K)){ stop("No valid Kernel specified") }
  
  eigobj <- NULL
  Utrunc <- NULL
  if(control$truncate) {
    full <- FALSE
    if (control$loss == "leastsquares") {
      full <- TRUE
    }
    truncDat <- Ktrunc(K = K, sigma=sigma, lastkeeper=control$lastkeeper, epsilon=control$epsilon,
                       full = full, quiet = control$quiet)
    Utrunc <- truncDat$Utrunc
    eigvals <- truncDat$eigvals
    if (full) {
      eigobj <- truncDat$eigobj
    }
  } else {
    if (control$loss == "leastsquares") {
      eigobj <- eigen(K)
    }
    Utrunc <- NULL
    eigvals <- NULL
  }
  
  lastkeeper <- NULL
  lastkeeper <- if(control$truncate) ncol(Utrunc)
  
  return(list(K = K,
              Utrunc = Utrunc,
              eigvals = eigvals,
              eigobj = eigobj,
              lastkeeper = lastkeeper))
}


## Function that returns truncated versions of the data if given
## todo: throws warning whenever n < 500 because it uses eigen, should we suppress?
#' @export
Ktrunc <- function(X=NULL, K=NULL, sigma=NULL, epsilon=NULL, lastkeeper=NULL, full=FALSE, quiet = TRUE){
  if(is.null(K)){
    X.sd <- apply(X, 2, sd)
    X.mu <- colMeans(X)
    K <- kern_gauss(scale(X), ifelse(is.null(sigma), ncol(X), sigma))
  }
  
  ## Todo: maybe later allow for full return of eigs_sym but not full
  ## eigen decomposition. Might be useful for lambda search for logistic
  if(full){
    eigobj <- eigen(K)
    
    lastkeeper <- ifelse(is.null(lastkeeper),
                         min(which(cumsum(eigobj$values)/nrow(K) > (1-epsilon))),
                         lastkeeper)
    
    Ktilde <- mult_diag(eigobj$vectors[, 1:lastkeeper], eigobj$values)
    
    if(!quiet) print(paste("lastkeeper=",lastkeeper))
    
    return(
      list(Ktilde=Ktilde,
           Utrunc=eigobj$vectors[, 1:lastkeeper, drop = F],
           eigvals=eigobj$values[1:lastkeeper, drop = F],
           eigobj=eigobj)
    )
    
  }
  if(is.null(lastkeeper)){
    #denoms <- c(10, 6, 3, 1) # we could also let people set this to speed up the algorithm
    #CJH: even at 10, with N=5000 it is too big. Let's try this: 
    if (nrow(K) <= 500) numvectorss=nrow(K) else numvectorss=c(50, 100, 200, 300,500,1000)
    enoughvar=FALSE
    j=1
    while (enoughvar==FALSE){
      numvectors=numvectorss[j]
      if(!quiet) print(paste("trying",numvectors,"vectors"))
      eigobj <- NULL
      eigobj <- try({eigs_sym(K, numvectors, which="LM")}, silent = T)
      
      #for now, letting it throw an error in certain failure cases.
      totalvar=sum(eigobj$values)/nrow(K)
      
      if (totalvar>=(1-epsilon)){
        lastkeeper = min(which(cumsum(eigobj$values)/nrow(K) > (1-epsilon)))
        Ktilde <- mult_diag(eigobj$vectors[, 1:lastkeeper, drop = F], eigobj$values)
        
        if(!quiet) print(paste("lastkeeper=",lastkeeper))
        enoughvar=TRUE
        
      }
      j=j+1
      # Will eventually add warning here for var < 0.99 at numvector= 300
      # Allow to proceed w/out truncation and SEs or to try full eigen
    } #end while gotenough==FALSE      
  } else { #end if is.null(lastkeeper)
    ## Suppress warning about using all eigenvalues
    eigobj <- suppressWarnings({eigs_sym(K, lastkeeper, which="LM")})
    Ktilde <- mult_diag(eigobj$vectors, eigobj$values)
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
newKernel <- function(X, newData, whichkernel = "gaussian") {
  
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
    newK <- new_gauss_kern(newData, X, ncol(X))
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