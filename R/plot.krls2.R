#' Plot method for KRLS and KRLogit Model Fits
#' 
#' @description Produces two types of plots. The first type of plot shows
#' histograms for the pointwise partial derivatives to examine the 
#' heterogeneity in the marginal effects of each predictor (\code{which}==1).
#' The second type of plot shows estimates of the conditional expectation
#' functions of \eqn{E[Y|X]} for each predictor (\code{which==2}). For each
#' plot, the predictor of interest varies from its 1st to its 3rd quartile values, while
#' the other predictors are kept at the means (or other values specified in 
#' \code{setx}). For binary varibales the \eqn{E[Y|X]} are predicted at the
#' max and the min value of the predictor (instead of the range from the 1st
#' to the 3rd quantile).
#'
#' @param obj an object of class \code{krls2} that preferably results from a call to 
#'  \code{\link{summary.krls2}}. Also accepts the output of \code{\link{krls}} but
#'  this is slower as it has to recompute pointwise marginal effects if the
#'  histograms of pointwise marginal effects are requested by \code{which}.
#' @param which if a subset of the plots is required, specify a subset of the numbers \code{1:2}.
#' @param main main title for histograms of pointwise partial derivatives.
#' @param setx either one of \code{mean} or \code{median} to hold other
#'  predictors at their mean or median values for the conditional expectation 
#'  plots. Alternativeley the user can specific a numeric vector with predictor
#'  values at which the other predictors should be fixed for the conditional
#'  expectation plots. If specifed in this way there must be one value per
#'  predictor and the order of the values much match the order of the predictor
#'  used in the predictor matrix of the krls fit passed in \code{obj}.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot,
#'  see \code{\link{par}} (\code{ask=.}).
#' @param nvalues scalar that specifies the number of values at which 
#'  conditional expectations should be plotted.
#' @param probs vector with numbers between 0 and 1 that specify the quantiles
#'  that determine the range for of the predictor values for which the 
#'  conditional expectation should be plotted. By default we vary each 
#'  predictor from the 1st quartile to the 3rd quartile value.
#' @param \dots additional arguments to be passed to lower level functions
#' 
#' @details Notice that the histograms for the partial derivatives can only
#'  be plotted if the KRLS object is the result of \code{\link{summary.krls2}}.
#'
#' @examples
#' # non-linear example
#' # set up data
#' N <- 200
#' x1 <- rnorm(N)
#' x2 <- rbinom(N,size=1,prob=.2)
#' y <- x1^3 + .5*x2 + rnorm(N,0,.15)
#' X <- cbind(x1,x2)
#' 
#' # fit model and summarize
#' krlsout <- summary(krls(X=X,y=y))
#' 
#' # plot marginal effects and conditional expectation plots
#' plot(krlsout)
#' 
#' # binary exmaple
#' y <- rbinom(
#'   N,
#'   size = 1,
#'   prob = 1/(1+exp(-y))
#' )
#' 
#' krlogout <- summary(
#'   krls(X=X,y=y,loss='logistic',epsilon=0.01)
#' )
#' 
#' plot(krlogout)
#' 
#' @export
#'
plot.krls2 <-
  function(obj,
           which = c(1:2),
           main = "distributions of pointwise marginal effects",
           setx = "mean",
           ask = prod(par("mfcol")) < nplots,
           nvalues = 50,
           probs = c(.25,.75),
           ...)
  {
    
    if( class(obj)!= "krls2" ){
      warning("obj not of class 'krls2'")
      #UseMethod("summary")
      return(invisible(NULL))
    }
    
    d <- ncol(obj$X)
    n <- nrow(obj$X)
    if(length(probs)!=2){
      stop("length(probs) must be 2")
    }
    
    # check setx
    if(is.numeric(setx)){
      if(length(setx)!=d){
        stop("length(setx) must be equal to number of predictors") 
      }
    } else {
      if(length(setx)!=1){stop("setx must be one of mean or median")}
      if(sum(setx %in% c("mean","median"))<1){stop("setx must be one of mean or median")}
      setx <- apply(obj$X,2,setx)
    }
    
    nplots <- 0
    if(1 %in% which){ nplots <- nplots + 1}
    if(2 %in% which){ nplots <- nplots + d}        
    
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    
    if(is.null(colnames(obj$X))){
      colnames(obj$X) <- paste("x",1:d,sep="") 
    } 
    
    # derivatives
    if(1 %in% which){ # histograms of partial derivatives
      if(is.null(obj$derivatives)){
        warning("running summary.krls2 because object passed does not have PWMFX.
                To save time, pass the result of summary.krls2 instead of krls to this plot method.")
        obj <- inference.krls2(obj)
      } else {
        colnames(obj$derivatives) <- colnames(obj$X)
        
      
        form <-  as.formula(paste("~",paste(colnames(obj$derivatives),collapse="+"),sep=""))
        require(lattice)
        
        print(histogram(form,
                        data=data.frame(obj$derivatives),
                        breaks=NULL,
                        main=main
                        ,...)
        )
        #if(length(which)!=1){readline("Press any key for next plot")}
      }
    }
    
    if(2 %in% which){  # conditional expectation plots
      lengthunique    <- function(x){length(unique(x))}
      # vector with positions of binary variables
      binaryindicator <- which(apply(obj$X,2,lengthunique)==2)
      quantiles <-  apply(obj$X,2,quantile,probs=probs)     
      
      for(i in 1:d){
        
        if(i %in% binaryindicator){ # E[Y|X] for binary Xs
          Xi <- c(min(obj$X[,i]),max(obj$X[,i]))
          Newdata <- matrix(rep(setx,2),ncol=d,byrow=T)
          Newdata[,i] <- Xi
          
        } else {
          # E[Y|X] plots for cont Xs
          Xi <- seq(quantiles[1,i],quantiles[2,i],length.out=nvalues)
          Newdata <- matrix(rep(setx,nvalues),ncol=d,byrow=T)
          Newdata[,i] <- Xi
        }
        pout <- predict(obj,newdata=Newdata,se=TRUE)
        Ylo  <- pout$fit-1.96*pout$se
        Yhi  <- pout$fit+1.96*pout$se
        plot(
          y=pout$fit,x=Xi,
          xlab=colnames(obj$X)[i],
          ylab=c("E[Y|X]"),
          ylim=c(min(Ylo) -.25*c(sqrt(var(pout$fit))),
                 max(Yhi))+.25*c(sqrt(var(pout$fit))),
          pch=19
        )
        arrows(x0=Xi,y0=Ylo,y1=Yhi,length = 0)
        
      }
    }  
    
  }