library(testthat)
library(BFRdist)

# --- Data Validation & Error Handling ---

test_that("Core functions reject out-of-bounds data", {
  a <- 1.5; b <- 2.0
  
  # Inputs must be strictly between 0 and 1
  expect_error(dbfr(0, a, b), "strictly greater than 0 and less than 1")
  expect_error(dbfr(1, a, b), "strictly greater than 0 and less than 1")
  expect_error(pbfr(-0.5, a, b), "strictly greater than 0 and less than 1")
  expect_error(qbfr(1.5, a, b), "strictly greater than 0 and less than 1")
})

test_that("Core functions reject NA and Inf", {
  a <- 1.5; b <- 2.0
  
  expect_error(dbfr(NA, a, b), "NA and Inf values are not permitted")
  expect_error(pbfr(Inf, a, b), "NA and Inf values are not permitted")
  expect_error(qbfr(NaN, a, b), "NA and Inf values are not permitted")
})

test_that("Core functions reject invalid shape parameters", {
  expect_error(dbfr(0.5, alpha = -1, beta = 2), "alpha and beta must be positive")
  expect_error(pbfr(0.5, alpha = 2, beta = 0), "alpha and beta must be positive")
})

# --- Mathematical Integrity ---

test_that("Densities are strictly non-negative", {
  x_seq <- seq(0.01, 0.99, length.out = 100)
  densities <- dbfr(x_seq, alpha = 2, beta = 3)
  
  expect_true(all(densities >= 0))
})

test_that("CDF bounds are exactly between 0 and 1", {
  q_seq <- seq(0.01, 0.99, length.out = 100)
  probabilities <- pbfr(q_seq, alpha = 2, beta = 3)
  
  expect_true(all(probabilities >= 0 & probabilities <= 1))
  
  # CDF should be monotonically increasing
  expect_true(all(diff(probabilities) > 0))
})

test_that("CDF and Quantile functions are inverses", {
  p_initial <- c(0.1, 0.5, 0.9)
  q_vals <- qbfr(p_initial, alpha = 1.5, beta = 2.0)
  p_recovered <- pbfr(q_vals, alpha = 1.5, beta = 2.0)
  
  # expect_equal allows for tiny floating-point math differences
  expect_equal(p_recovered, p_initial)
})

test_that("Random generation produces valid bounds", {
  set.seed(42)
  r_vals <- rbfr(100, alpha = 2, beta = 3)
  
  expect_length(r_vals, 100)
  expect_true(all(r_vals > 0 & r_vals < 1))
})

# --- Object-Oriented (S3) Implementation ---

test_that("fit_bfr creates correct S3 objects", {
  set.seed(123)
  sim_data <- rbfr(50, alpha = 1.5, beta = 2.5)
  
  fit_obj <- fit_bfr(sim_data, method = "mle")
  
  # Verify class and structure
  expect_s3_class(fit_obj, "bfr")
  expect_type(fit_obj$estimate, "double")
  expect_named(fit_obj$estimate, c("alpha", "beta"))
  
  # Verify predict method
  pred <- predict(fit_obj, n.ahead = 5)
  expect_length(pred, 5)
  expect_true(all(pred > 0 & pred < 1))
  
  # Verify summary method
  sum_obj <- summary(fit_obj)
  expect_s3_class(sum_obj, "summary.bfr")
  expect_true(is.numeric(sum_obj$loglik))
  expect_true(is.numeric(sum_obj$aic))
})

# --- Order Statistic Validation ---

test_that("Galton skewness and Moors kurtosis handle edge cases", {
  x <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  
  sk <- galton_skewness(x)
  ku <- moors_kurtosis(x)
  
  expect_true(is.numeric(sk))
  expect_true(is.numeric(ku))
  
  # Test zero variance handling
  flat_data <- rep(0.5, 10)
  expect_warning(galton_skewness(flat_data), "Interquartile range is zero")
})

# --- S3 Print and Plot Methods ---

test_that("S3 print methods produce expected console output", {
  set.seed(42)
  fit_obj <- fit_bfr(rbfr(50, 1.5, 2.5), method = "mle")
  
  # expect_output captures the console text to verify it prints correctly
  expect_output(print(fit_obj), "BFR Distribution Fit")
  expect_output(print(summary(fit_obj)), "BFR Distribution Fit Summary")
  expect_output(print(summary(fit_obj)), "Log-Likelihood")
})

test_that("S3 plot method executes without errors", {
  set.seed(42)
  fit_obj <- fit_bfr(rbfr(50, 1.5, 2.5), method = "mle")
  
  # Use a very small B for the test so it runs instantly
  # expect_invisible verifies the function completes and returns its data invisibly
  expect_invisible(plot(fit_obj, B = 5, show_observed = TRUE))
})

# --- Alternative Methods and Edge Cases ---

test_that("Method of Moments (MME) fitting executes correctly", {
  set.seed(42)
  fit_obj <- fit_bfr(rbfr(50, 1.5, 2.5), method = "mme")
  
  expect_s3_class(fit_obj, "bfr")
  expect_equal(fit_obj$method, "mme")
  expect_true(is.numeric(fit_obj$estimate["alpha"]))
})

test_that("Predict method correctly rejects invalid inputs", {
  set.seed(42)
  fit_obj <- fit_bfr(rbfr(50, 1.5, 2.5), method = "mle")
  
  # Base R's S3 dispatch handles class checking automatically
  expect_error(predict(fit_obj, n.ahead = 0), "n.ahead must be at least 1")
})

test_that("Summary statistics handle NA values when requested", {
  x_with_na <- c(0.1, 0.2, 0.3, 0.4, 0.5, NA)
  
  # Should successfully calculate by ignoring the NA
  expect_true(is.numeric(galton_skewness(x_with_na, na.rm = TRUE)))
  expect_true(is.numeric(moors_kurtosis(x_with_na, na.rm = TRUE)))
})

# --- Direct C++ Wrapper Validation ---

test_that("estimateBFR directly validates inputs before C++ execution", {
  # Non-numeric data
  expect_error(estimateBFR(c("a", "b")), "must be a numeric vector")
  
  # Insufficient length
  expect_error(estimateBFR(0.5), "At least two observations are required")
  
  # Missing values (NA)
  expect_error(estimateBFR(c(0.5, NA, 0.7)), "Missing values are not allowed")
  
  # Infinite/NaN values
  expect_error(estimateBFR(c(0.5, Inf, 0.7)), "All observations must be finite")
  expect_error(estimateBFR(c(0.5, NaN, 0.7)), "Missing values are not allowed")
  
  # Out-of-bounds values (<= 0 or >= 1)
  expect_error(estimateBFR(c(0, 0.5)), "strictly between 0 and 1")
  expect_error(estimateBFR(c(0.5, 1)), "strictly between 0 and 1")
  expect_error(estimateBFR(c(-0.5, 0.5)), "strictly between 0 and 1")
})

test_that("MLE recovers known BFR parameters", {

  set.seed(123)

  alpha_true <- 2
  beta_true <- 3

  # Generate observations from a known BFR distribution
  x <- rbfr(
    n = 5000,
    alpha = alpha_true,
    beta = beta_true
  )

  # Fit using maximum likelihood
  fit <- fit_bfr(x, method = "mle")

  # Check that the result has the expected structure
  expect_s3_class(fit, "bfr")

  # Check that MLE was actually used
  expect_equal(fit$method, "mle")

  # Check that both parameters were estimated
  expect_named(
    fit$estimate,
    c("alpha", "beta")
  )

  # Check that estimates are positive
  expect_gt(fit$estimate["alpha"], 0)
  expect_gt(fit$estimate["beta"], 0)

  # Check that estimates are reasonably close to the true values
  expect_equal(
    fit$estimate[["alpha"]],
    alpha_true,
    tolerance = 0.15
  )

  expect_equal(
    fit$estimate[["beta"]],
    beta_true,
    tolerance = 0.15
  )
})


test_that("MLE optimisation completes successfully", {

  set.seed(123)

  x <- rbfr(
    n = 1000,
    alpha = 2,
    beta = 3
  )

  fit <- fit_bfr(x, method = "mle")

  # optim() convergence code of 0 means successful convergence
  expect_equal(fit$optim_results$convergence, 0)

  # Negative log-likelihood should be finite
  expect_true(is.finite(fit$optim_results$value))
})


test_that("MLE rejects observations outside the BFR support", {
  x <- c(0.1, 0.2, 0.5, 1.0, 0.8)
  expect_error(fit_bfr(x, method = "mle"), "strictly greater than 0 and less than 1")
})

test_that("MLE rejects observations below zero", {
  x <- c(0.1, 0.2, 0.5, -0.3, 0.8)
  expect_error(fit_bfr(x, method = "mle"), "strictly greater than 0 and less than 1")
})