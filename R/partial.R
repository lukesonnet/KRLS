
#' Plot partial dependence of the  marginal effects
#' 
#' @param obj A \code{krls2} object
#' @param xcol An integer specifying the column of the data we should plot the partial dependence function of (this is the S variable)
#' @param derivcol An integer that specifies what kind of partial dependence plot should be generated. If left blank, the conditional expectation function is plotted; if specified the partial derivative of the CEF with respect to that column in the data is plotted.
#' @param nlocations An integer with the number of locations to plot the function at, defaults to 30 over the range of the relevant variable
#' 
#' @export
partial_dep_mfx_plot <- function(obj, xcol, derivcol = NULL, nlocations = 100) {
  
  if (obj$kernel != 'gaussian' | obj$loss != "leastsquares") {
    stop(
      "Currently only supports partial dependence plots for ",
      "gaussian kernels and least squares loss functions"
    )
  }
  
  if (!is.numeric(derivcol)) {
    stop(
      "Currently only support partial dependence of derivatives of the CEF"
    )
  }
  
  minX <- min(obj$X[, xcol])
  maxX <- max(obj$X[, xcol])
  partial_dat <- data.frame(
    xstar = seq(minX, maxX, length.out = nlocations),
    xstar_deriv = NA,
    xstar_deriv_sd = NA
  )
  
  for (i in seq_along(partial_dat$xstar)) {
    Xstar <- obj$X
    kro <- partial(
      obj, 
      xcol = xcol, 
      derivcol = derivcol, 
      Xs = partial_dat$xstar[i]
    )
    partial_dat$xstar_deriv[i] <- kro$avg_deriv
    partial_dat$xstar_deriv_sd[i] <- sqrt(kro$var_avg_deriv)
  }
  
  partial_dat$xstar_deriv_low <- partial_dat$xstar_deriv - qt(0.975, df = nrow(obj$X) - 1) * partial_dat$xstar_deriv_sd
  partial_dat$xstar_deriv_high <- partial_dat$xstar_deriv + qt(0.975, df = nrow(obj$X) - 1) * partial_dat$xstar_deriv_sd
  ymax <- max(c(partial_dat$xstar_deriv_high, obj$derivatives[, xcol]))
  ymin <- min(c(partial_dat$xstar_deriv_low, obj$derivatives[, xcol]))
  
  with(
    partial_dat,
    plot(
      x = xstar,
      y = xstar_deriv,
      type = "l",
      ylim = c(ymin, ymax),
      ylab = paste0("dy/(d", colnames(obj$X)[derivcol], ")"),
      xlab = colnames(obj$X)[xcol]
    )
  )
  points(x = obj$X[, xcol], y = obj$derivatives[, xcol], col = "red")
  with(partial_dat, segments(x0 = xstar, y0 = xstar_deriv_low, y1 = xstar_deriv_high))
  
  return(invisible(partial_dat))
}

#' @export
partial_dep_mfx <- function(obj, xcol, derivcol, Xs) {
  
  Xstar <- obj$X
  Xstar[, xcol] <- Xs
  
  sdX <- apply(obj$X, 2, sd)
  mX <- apply(obj$X, 2, mean)
  X <- scale(obj$X)
  K <- t(newKernel(obj$X, newData = Xstar, b = obj$b, whichkernel = obj$kernel))
  Xstarscale <- scale(Xstar, center = mX, scale = sdX)
  
  kro <- pwmfx(
    Xstar = Xstarscale, 
    X = X, 
    K = K, 
    wrt_column = derivcol - 1, 
    coefhat = obj$coeffs,
    Vcovc = obj$vcov.c, 
    p = rep(1, nrow(X)),
    b = obj$b,
    computevarderiv = FALSE,
    computederiv2 = FALSE
  )
  
  sdy <- sd(obj$y)
  kro$deriv <- kro$deriv * sdy / sdX[derivcol]
  kro$avg_deriv <- mean(kro$deriv)
  kro$var_avg_deriv <- kro$var_avg_deriv * (sdy/sdX[derivcol])^2
  
  kro
}
