#' Summary function for KRLS2 objects. 
#' 
#' @description Main method for examining results of krls,
#' much as one would use summary.lm for lm objects. Like its lm analog,
#' it provides a summary akin to a regression table: average marginal 
#' effects for each input variables, together with standard errors.
#' 
#' @export
summary.krls2 <- function(object,
                          probs = c(.25, .5, .75),
                          ...) {
  if( class(object)!= "krls2" ) {
    warning("Object not of class 'krls2'")
    UseMethod("summary")
    return(invisible(NULL))
  }
        
  cat("* *********************** *\n")
  cat("Model Summary:\n\n")
  #cat("R2:",object$R2,"\n\n")
      
  d <- ncol(object$X)
  n <- nrow(object$X)
        
  coefficients <- matrix(NA,d,0)
  rownames(coefficients) <- colnames(object$X) 
              
  if(is.null(object$derivatives)){
    object <- inference.krls2(object, ...)
  } 
  # average marginal effects  
  est     <- t(object$avgderivatives)
  se     <- sqrt((object$var.avgderivatives))
  tval   <- est/(se)
  avgcoefficients <- cbind(est, se, tval, 2 * pt(abs(tval),n-d, lower.tail = FALSE))
  colnames(avgcoefficients) <- c("Est", "Std. Error", "t value", "Pr(>|t|)")
        
  #todo: check for underflow of pval
  
  # add stars for binary    
  if(sum(object$binaryindicator)>0){         
    rownames(avgcoefficients)[object$binaryindicator] <- paste(rownames(avgcoefficients)[object$binaryindicator],"*",sep="")
  }
        
  cat("Average Marginal Effects:\n")
  print(avgcoefficients,...)
  if(sum(object$binaryindicator)>0){
    cat("\n(*) average dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
  }
        
  # quantiles of derivatives
  qderiv <- apply(object$derivatives,2,quantile,probs=probs)
  if(sum(object$binaryindicator)>0){         
    colnames(qderiv)[object$binaryindicator] <- paste(colnames(qderiv)[object$binaryindicator],"*",sep="")
  }
  qderiv <- t(qderiv)
        
  cat("\n")
  cat("Quartiles of Marginal Effects:\n")
  print(qderiv,...)
        
  if(sum(object$binaryindicator)>0){
   cat("\n(*) quantiles of dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
  }
                     
  object$sumavgderiv <- avgcoefficients
  object$qderiv <- qderiv

  return(invisible(as.list(object)))
}

