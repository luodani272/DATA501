#' Quantile-Based Estimation of the BFR Distribution Parameters
#'
#' Estimates the two parameters of the BFR distribution using
#' sample quartiles.
#'
#' The estimator uses the sample first and third quartiles and
#' the quantile function of the BFR distribution to obtain
#' estimates of alpha and beta.
#'
#' @param x Numeric vector containing observations from the BFR
#'   distribution. All observations must lie strictly between 0
#'   and 1.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{alpha}{The estimated value of alpha.}
#'   \item{beta}{The estimated value of beta.}
#'   \item{q25}{The sample first quartile.}
#'   \item{q75}{The sample third quartile.}
#' }
#'
#' @details
#' The BFR quantile function is
#'
#' \deqn{
#' Q(u; alpha, beta) =
#' \left[
#' 1 + alpha
#' \tan\left(\frac{\pi u}{2}\right)^{-1/beta}
#' \right]^{-1}.
#' }
#'
#' Rearranging the quantile function gives
#'
#' \deqn{
#' \log\left(\frac{1-x}{x}\right)
#' =
#' \log(alpha)
#' -
#' \frac{1}{beta}
#' \log\left(\tan\left(\frac{\pi u}{2}\right)\right).
#' }
#'
#' The first and third sample quartiles are substituted into
#' this relationship to estimate alpha and beta.
#'
#' Sample quartiles are calculated using linear interpolation
#' corresponding to quantile type 7.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' x <- rBFR(
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