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
    fit$estimate["alpha"],
    alpha_true,
    tolerance = 0.15
  )

  expect_equal(
    fit$estimate["beta"],
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
  expect_error(fit_bfr(x, method = "mle"), "strictly between 0 and 1")
})

test_that("MLE rejects observations below zero", {
  x <- c(0.1, 0.2, 0.5, -0.3, 0.8)
  expect_error(fit_bfr(x, method = "mle"), "strictly between 0 and 1")
})