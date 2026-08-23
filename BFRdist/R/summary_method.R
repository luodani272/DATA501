#' Summary Method for BFR Model Fits
#'
#' Outputs standard errors and goodness-of-fit metrics for a fitted BFR object.
#'
#' @param object An object of class "bfr".
#' @param ... Additional arguments affecting the summary produced.
#' @return An object of class "summary.bfr".
#' @export
summary.bfr <- function(object, ...) {
  x <- object$data
  n <- length(x)
  a <- unname(object$estimate["alpha"])
  b <- unname(object$estimate["beta"])
  
  # Calculate Log-Likelihood
  pdf_vals <- dbfr(x, alpha = a, beta = b)
  # Filter out zeros to prevent log(0) issues
  ll <- sum(log(pdf_vals[pdf_vals > 0])) 
  
  # Goodness-of-Fit Metrics
  k <- 2 # number of parameters
  aic <- 2 * k - 2 * ll
  bic <- k * log(n) - 2 * ll
  
  # Approximate Standard Errors using the Hessian matrix (for MLE only)
  se <- c(alpha = NA_real_, beta = NA_real_)
  if (object$method == "mle") {
    nll <- function(par) {
      -sum(log(dbfr(x, par[1], par[2])))
    }
    # Calculate Hessian numerically
    hess <- stats::optimHess(c(a, b), nll)
    
    # Invert the Hessian to get the covariance matrix
    cov_mat <- tryCatch(solve(hess), error = function(e) matrix(NA, 2, 2))
    se <- sqrt(diag(cov_mat))
    names(se) <- c("alpha", "beta")
  }
  
  # Bundle into a summary object
  res <- list(
    estimate = object$estimate,
    se = se,
    method = object$method,
    loglik = ll,
    aic = aic,
    bic = bic,
    n = n
  )
  class(res) <- "summary.bfr"
  return(res)
}

