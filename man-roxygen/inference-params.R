#' @param derivative boolean, whether to return the PWMFX. Defaults to TRUE.
#' @param derivativevar boolean, whether to return variance of the individual PWMFX. Defaults to FALSE.
#' @param vcov boolean, whether to return the variance of the AME. Defaults to TRUE.
#' @param robust boolean, whether to use the sandwich estimator for either leastsquares or logistic loss. Defaults to FALSE.
#' @param clusters a vector of any type representing the cluster each observation belongs to.
#' @param return_components boolean, whether to return some of the component pieces of the PWMFX and variance estimates
#' @param lambda_meat boolean, whether to include the normalized lambda in the meat of the robust sandwich estimators. Defaults to TRUE
