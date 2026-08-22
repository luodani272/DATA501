#' Print Summary of a BFR Distribution Fit
#'
#' Print method for objects of class "summary.bfr". Displays parameter estimates, 
#' standard errors, and goodness-of-fit metrics including Log-Likelihood, AIC, 
#' and BIC (when applicable).
#'
#' @param x An object of class "summary.bfr", typically the result of `summary.bfr()`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the summary object `x`.
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