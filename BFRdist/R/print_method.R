#' Print a BFR Distribution Fit
#'
#' An S3 method that prints a clean, concise view of a fitted BFR distribution.
#'
#' @details
#' This function provides a streamlined view of a \code{bfr} object, displaying 
#' only the basic fitted parameter estimates and the estimation method used[cite: 1, 2].
#'
#' @param x An object of S3 class \code{bfr}, typically the result of \code{fit_bfr}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the original \code{bfr} object \code{x}.
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
#' It explicitly formats and outputs the standard errors and goodness-of-fit metrics 
#' (including Log-Likelihood, AIC, and BIC) derived from the fitted model.
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
#' # 1. Generate synthetic GRVI-like data
#' set.seed(123)
#' raw_data <- rbfr(n = 100, alpha = 1.5, beta = 2.5)
#' 
#' # 2. Fit the BFR model to the data
#' bfr_model <- fit_bfr(raw_data, method = "mle")
#' 
#' # 3. Generate and print the summary (implicitly calls print.summary.bfr)
#' model_summary <- summary(bfr_model)
#' print(model_summary)
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