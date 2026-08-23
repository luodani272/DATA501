#' Predict Method for BFR Model Fits
#'
#' Simulates new predictions based on the fitted BFR model.
#'
#' @param object An object of class "bfr".
#' @param n.ahead Integer. The number of new observations to simulate. Default is 1.
#' @param ... Additional arguments.
#' @return A numeric vector of simulated values from the fitted BFR distribution.
#' @export
predict.bfr <- function(object, n.ahead = 1, ...) {
  if (n.ahead < 1) {
    stop("n.ahead must be at least 1.")
  }
  
  a <- unname(object$estimate["alpha"])
  b <- unname(object$estimate["beta"])
  
  # Simulate new data using the random generation function
  rbfr(n = n.ahead, alpha = a, beta = b)
}