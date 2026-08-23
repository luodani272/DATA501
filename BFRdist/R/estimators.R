#' Fit the BFR distribution
#'
#' Estimates parameters of the BFR distribution using Maximum Likelihood (mle) 
#' or Method of Moments (mme).
#'
#' @param x A numeric vector of observations strictly between 0 and 1.
#' @param method A character string specifying the estimation method: "mle" or "mme".
#' @return An object of S3 class \code{bfr}.
#' @export
fit_bfr <- function(x, method = c("mle", "mme")) {
  method <- match.arg(method)
  
  # Use Daniel's C++ quantile estimator for initial starting values
  init_est <- estimateBFR(x) 
  init_alpha <- init_est["alpha"]
  init_beta <- init_est["beta"]
  
  # Maximum Likelihood Estimation
  if (method == "mle") {
    nll <- function(par) {
      a <- par[1]; b <- par[2]
      if (a <= 0 || b <= 0) return(Inf)
      pdf_vals <- dbfr(x, alpha = a, beta = b)
      if (any(pdf_vals <= 0, na.rm = TRUE)) return(Inf)
      -sum(log(pdf_vals))
    }
    opt <- stats::optim(par = c(init_alpha, init_beta), fn = nll, 
                 method = "L-BFGS-B", lower = c(1e-5, 1e-5))
    
    # Method of Moments Estimation
  } else if (method == "mme") {
    m1 <- mean(x)
    m2 <- mean(x^2)
    
    moment_obj <- function(par) {
      a <- par[1]; b <- par[2]
      if (a <= 0 || b <= 0) return(Inf)
      
      e_x <- tryCatch(stats::integrate(function(u) u * dbfr(u, a, b), 0, 1)$value, error = function(e) NA)
      e_x2 <- tryCatch(stats::integrate(function(u) u^2 * dbfr(u, a, b), 0, 1)$value, error = function(e) NA)
      
      if (is.na(e_x) || is.na(e_x2)) return(Inf)
      (e_x - m1)^2 + (e_x2 - m2)^2
    }
    opt <- stats::optim(par = c(init_alpha, init_beta), fn = moment_obj, 
                 method = "L-BFGS-B", lower = c(1e-5, 1e-5))
  }
  
  # Construct and return the S3 object
  out <- list(
    estimate = c(alpha = unname(opt$par[1]), beta = unname(opt$par[2])),
    method = method,
    data = x,
    optim_results = opt
  )
  class(out) <- "bfr"
  return(out)
}