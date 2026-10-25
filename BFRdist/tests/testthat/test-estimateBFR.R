library(testthat)
library(BFRdist)

# --- Consolidated Synthetic Data Setup ---
# Defining these once at the top speeds up test execution
set.seed(42)
sim_data_small <- rbfr(50, alpha = 1.5, beta = 2.5)
sim_data_large <- rbfr(1000, alpha = 2.0, beta = 3.0)
flat_data <- rep(0.5, 10)


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
  expect_length(sim_data_small, 50)
  expect_true(all(sim_data_small > 0 & sim_data_small < 1))
})


# --- Object-Oriented (S3) Implementation ---

test_that("fit_bfr creates correct S3 objects and summary outputs", {
  fit_obj <- fit_bfr(sim_data_small, method = "mle")
  
  # Verify class and structure
  expect_s3_class(fit_obj, "bfr")
  expect_type(fit_obj$estimate, "double")
  expect_named(fit_obj$estimate, c("alpha", "beta"))
  
  # Verify predict method
  pred <- predict(fit_obj, n.ahead = 5)
  expect_length(pred, 5)
  expect_true(all(pred > 0 & pred < 1))
  
  # Verify summary method and deep structure
  sum_obj <- summary(fit_obj)
  expect_s3_class(sum_obj, "summary.bfr")
  expect_true(is.numeric(sum_obj$loglik))
  expect_true(is.numeric(sum_obj$aic))
  expect_true(is.numeric(sum_obj$se["alpha"]))
  expect_true(is.numeric(sum_obj$se["beta"]))
})

test_that("S3 print methods produce expected console output", {
  fit_obj <- fit_bfr(sim_data_small, method = "mle")
  
  # expect_output captures the console text to verify it prints correctly
  expect_output(print(fit_obj), "BFR Distribution Fit")
  expect_output(print(summary(fit_obj)), "BFR Distribution Fit Summary")
  expect_output(print(summary(fit_obj)), "Log-Likelihood")
})

test_that("S3 plot method executes without errors", {
  fit_obj <- fit_bfr(sim_data_small, method = "mle")
  
  # Use a very small B for the test so it runs instantly
  # expect_invisible verifies the function completes and returns its data invisibly
  expect_invisible(plot(fit_obj, B = 5, show_observed = TRUE))
})


# --- Alternative Methods and Edge Cases ---

test_that("fit_bfr successfully fits using Method of Moments (MME)", {
  fit_obj <- fit_bfr(sim_data_small, method = "mme")
  
  expect_s3_class(fit_obj, "bfr")
  expect_equal(fit_obj$method, "mme")
  expect_true(is.numeric(fit_obj$estimate["alpha"]))
})

test_that("fit_bfr rejects unsupported methods", {
  expect_error(fit_bfr(sim_data_small, method = "bayesian"))
})

test_that("fit_bfr rejects zero-variance data", {
  expect_error(fit_bfr(flat_data, method = "mle"))
})

test_that("Predict method correctly rejects invalid inputs", {
  fit_obj <- fit_bfr(sim_data_small, method = "mle")
  
  # Base R's S3 dispatch handles class checking automatically
  expect_error(predict(fit_obj, n.ahead = 0), "n.ahead must be at least 1")
})


# --- Order Statistic Validation ---

test_that("Galton skewness and Moors kurtosis handle edge cases", {
  sk <- galton_skewness(sim_data_small)
  ku <- moors_kurtosis(sim_data_small)
  
  expect_true(is.numeric(sk))
  expect_true(is.numeric(ku))
  
  # Test zero variance handling
  expect_warning(galton_skewness(flat_data), "Interquartile range is zero")
})

test_that("Summary statistics handle NA values when requested", {
  x_with_na <- c(sim_data_small, NA)
  
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


# --- Empirical Data & MLE Validation ---

test_that("fit_bfr rejects unfiltered real-world data with exact boundary values", {
  data("grvi_data", package = "BFRdist")
  
  # Spring Wheat GRVI on August 27 contains exact 0.0 values
  wheat_aug27 <- grvi_data$GRVI_27Aug[grvi_data$Crop == "Spring Wheat"]
  
  expect_error(fit_bfr(wheat_aug27, method = "mle"), "strictly greater than 0 and less than 1")
})

test_that("fit_bfr works accurately on cleaned real-world GRVI data", {
  data("grvi_data", package = "BFRdist")
  
  # Prepare the data exactly as documented
  canola_july <- grvi_data$GRVI_24_july[grvi_data$Crop == "Canola"]
  canola_clean <- canola_july[canola_july > 0 & canola_july < 1]
  
  # Fit the model
  model_mle <- fit_bfr(canola_clean, method = "mle")
  
  # Test that it returns the correct S3 class and doesn't throw errors
  expect_s3_class(model_mle, "bfr")
  expect_true(is.numeric(model_mle$estimate))
  expect_length(model_mle$estimate, 2)
})

test_that("MLE recovers known BFR parameters", {
  # Use the large pre-defined dataset to ensure asymptotic convergence
  fit <- fit_bfr(sim_data_large, method = "mle")
  
  # Check that the result has the expected structure
  expect_s3_class(fit, "bfr")
  expect_equal(fit$method, "mle")
  expect_named(fit$estimate, c("alpha", "beta"))
  
  # Check that estimates are positive
  expect_gt(fit$estimate["alpha"], 0)
  expect_gt(fit$estimate["beta"], 0)
  
  # Check that estimates are reasonably close to the true values (alpha=2, beta=3)
  expect_equal(fit$estimate[["alpha"]], 2.0, tolerance = 0.15)
  expect_equal(fit$estimate[["beta"]], 3.0, tolerance = 0.15)
})

test_that("MLE optimisation completes successfully", {
  fit <- fit_bfr(sim_data_large, method = "mle")
  
  # optim() convergence code of 0 means successful convergence
  expect_equal(fit$optim_results$convergence, 0)
  
  # Negative log-likelihood should be finite
  expect_true(is.finite(fit$optim_results$value))
})