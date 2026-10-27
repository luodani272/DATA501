#' Print Method for BFR Distribution Fits
#'
#' Prints a concise summary of the fitted BFR model, displaying the estimated 
#' shape parameters and the estimation method used.
#'
#' @details
#' This function provides a clean terminal output for objects of S3 class \code{bfr}. 
#' It displays the estimation method—either Quantile Estimation (QE) or Maximum 
#' Likelihood Estimation (MLE)—the number of valid observations used in the fit, 
#' and the final estimated values for the shape parameters \eqn{\alpha} and \eqn{\beta}.
#'
#' @param x An object of class \code{bfr}, typically the result of a call to \code{\link{fit_bfr}}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the original \code{bfr} object.
#'
#' @seealso \code{\link{fit_bfr}} for fitting the BFR distribution.
#'
#' @examples
#' \donttest{
#' # Generate synthetic BFR data
#' set.seed(123)
#' simulated_data <- rbfr(n = 250, alpha = 1.5, beta = 2.5)
#' 
#' # Fit the model using the high-performance Quantile Estimator (QE)
#' model_qe <- fit_bfr(simulated_data, method = "qe")
#' print(model_qe)
#' 
#' # Fit the model using Maximum Likelihood Estimation (MLE)
#' model_mle <- fit_bfr(simulated_data, method = "mle")
#' print(model_mle)
#' }
#'
#' @method print bfr
#' @export
print.bfr <- function(x, ...) {
  cat("\nBFR Distribution Fit\n")
  cat("--------------------\n")
  cat("Estimation Method:", toupper(x$method), "\n\n")
  cat("Parameter Estimates:\n")
  print(x$estimate)
  cat("\n")
  invisible(x)
}

#' Print Summary of a BFR Distribution Fit
#'
#' An S3 method that prints the formatted summary of a fitted BFR distribution.
#'
#' @details
#' This function displays the summary output for objects of class \code{summary.bfr}. 
#' The printed output adapts depending on the estimation method used:
#' \itemize{
#'   \item \strong{Maximum Likelihood Estimation (MLE):} Displays parameter estimates, 
#'     their approximated standard errors, and full goodness-of-fit metrics 
#'     (Log-Likelihood, AIC, and BIC).
#'   \item \strong{Quantile Estimator (QE):} Displays parameter estimates and the 
#'     Log-Likelihood. Standard errors are shown as \code{NA} (since QE is a closed-form 
#'     algebraic estimator), and information criteria (AIC/BIC) are cleanly omitted.
#' }
#'
#' @param x An object of S3 class \code{summary.bfr}, typically the result of calling \code{summary()} on a \code{bfr} object.
#' @param ... Additional arguments passed to or from other methods. (Currently unused).
#'
#' @return Invisibly returns the original \code{summary.bfr} object \code{x}.
#' 
#' @seealso \code{\link{summary.bfr}} for generating the summary object, and \code{\link{fit_bfr}} for model fitting.
#'
#' @examples
#' \donttest{
#' # Generate synthetic data
#' set.seed(123)
#' raw_data <- rbfr(n = 100, alpha = 1.5, beta = 2.5)
#' 
#' # Fit and print summary for the high-speed Quantile Estimator (QE)
#' bfr_qe <- fit_bfr(raw_data, method = "qe")
#' print(summary(bfr_qe))
#' 
#' # Fit and print summary for Maximum Likelihood Estimation (MLE)
#' bfr_mle <- fit_bfr(raw_data, method = "mle")
#' print(summary(bfr_mle))
#' }
#'
#' @method print summary.bfr
#' @export
print.summary.bfr <- function(x, ...) {
  cat("\nBFR Distribution Fit Summary\n")
  cat("----------------------------\n")
  cat("Estimation Method:", toupper(x$method), "\n")
  cat("Observations:", x$n, "\n\n")
  
  # Display parameters and standard errors
  coef_mat <- cbind(Estimate = x$estimate, `Std. Error` = x$se)
  print(coef_mat)
  
  cat("\nGoodness-of-Fit:\n")
  cat("  Log-Likelihood:", round(x$loglik, 4), "\n")
  if (x$method == "mle") {
    cat("  AIC:", round(x$aic, 4), "\n")
    cat("  BIC:", round(x$bic, 4), "\n")
  }
  cat("\n")
  invisible(x)
}