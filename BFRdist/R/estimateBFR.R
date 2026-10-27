#' Quantile-Based Estimation of the BFR Distribution Parameters
#'
#' Estimates the two parameters of the BFR distribution using sample quartiles.
#'
#' @details
#' The estimator uses the sample first and third quartiles and the quantile 
#' function of the BFR distribution to obtain estimates of \eqn{\alpha} and \eqn{\beta}.
#' 
#' The BFR quantile function is given by:
#' \deqn{Q_{BFR}(u;\alpha,\beta) = \frac{1}{1+\alpha\left(\tan\left(\frac{\pi u}{2}\right)\right)^{-1/\beta}}}
#' 
#' Rearranging the quantile function (where \eqn{x = Q_{BFR}(u)}) gives:
#' \deqn{\log\left(\frac{1-x}{x}\right) = \log(\alpha) - \frac{1}{\beta}\log\left(\tan\left(\frac{\pi u}{2}\right)\right)}
#' 
#' The first and third sample quartiles are substituted into this relationship to estimate \eqn{\alpha} and \eqn{\beta}.
#' Sample quartiles are calculated using linear interpolation corresponding to quantile type 7.
#'
#' @param x Numeric vector containing observations from the BFR distribution. All observations must lie strictly between 0 and 1.
#'
#' @return A named numeric vector containing:
#' \item{alpha}{The estimated value of \eqn{\alpha}.}
#' \item{beta}{The estimated value of \eqn{\beta}.}
#' \item{q25}{The sample first quartile.}
#' \item{q75}{The sample third quartile.}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' 
#' x <- rbfr(
#'   n = 1000,
#'   alpha = 2,
#'   beta = 3
#' )
#' 
#' estimateBFR(x)
#' }
#'
#' @export
estimateBFR <- function(x) {

  # ------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------

  if (!is.numeric(x)) {
    stop("x must be a numeric vector.")
  }

  if (length(x) < 2) {
    stop("At least two observations are required.")
  }

  if (anyNA(x)) {
    stop("Missing values are not allowed.")
  }

  if (any(!is.finite(x))) {
    stop("All observations must be finite.")
  }

  if (any(x <= 0 | x >= 1)) {
    stop("All observations must be strictly between 0 and 1.")
  }

  # ------------------------------------------------------------
  # Call the C++ estimator
  # ------------------------------------------------------------

  estimateBFR_cpp(x)
}