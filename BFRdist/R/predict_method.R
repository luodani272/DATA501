#' Simulate New Predictions from a Fitted BFR Model
#'
#' An S3 method that simulates new observations based on the estimated 
#' shape parameters of a fitted BFR model.
#'
#' @details
#' This function leverages the random generation utility (\code{rbfr}) to 
#' simulate new predictions from the fitted model. Because the BFR 
#' distribution is explicitly designed for double-bounded data such as the 
#' Generalized Radar Vegetation Index (GRVI), all simulated predictions will 
#' strictly fall within the unit interval (0, 1).
#'
#' @param object An object of S3 class \code{bfr}, typically the result of a call to \code{fit_bfr}.
#' @param n.ahead A single positive integer specifying the number of new observations to simulate. Default is 1.
#' @param ... Additional arguments passed to or from other methods. (Currently unused).
#' 
#' @return A numeric vector of length \code{n.ahead} containing simulated values from the fitted BFR distribution.
#' 
#' @seealso \code{\link{fit_bfr}} for fitting the model, and \code{\link{BFRdist}} for the underlying random generation function.
#' 
#' @examples
#' \donttest{
#' # 1. Generate synthetic GRVI-like data
#' set.seed(123)
#' raw_data <- rbfr(n = 100, alpha = 1.5, beta = 2.5)
#' 
#' # 2. Fit the BFR model to the data
#' bfr_model <- fit_bfr(raw_data, method = "mle")
#' 
#' # 3. Simulate 5 new predictions based on the fitted parameters
#' new_predictions <- predict(bfr_model, n.ahead = 5)
#' print(new_predictions)
#' }
#' 
#' @method predict bfr
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