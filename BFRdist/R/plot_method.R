#' SK_G x KU_M Bootstrap Plot
#'
#' Generates a parametric bootstrap distribution of Galton skewness (SK_G)
#' and Moors kurtosis (KU_M) for the BFR distribution and plots the
#' bootstrap replications together with confidence ellipses.
#'
#' @param data A numeric vector of observations from the BFR distribution.
#' @param alpha The alpha parameter of the BFR distribution.
#' @param beta The beta parameter of the BFR distribution.
#' @param B The number of bootstrap replications. Default is 1000.
#' @param conf A numeric vector containing the confidence levels for the
#'   confidence ellipses. Default is c(0.90, 0.95).
#' @param seed An optional numeric value used to set the random seed.
#'   Default is NULL.
#' @param show_observed Logical. If TRUE, the SK_G and KU_M values calculated
#'   from the original data are displayed on the plot. Default is TRUE.
#'
#' @return Invisibly returns a list containing the bootstrap statistics,
#'   observed statistics, bootstrap centre, confidence levels, parameters,
#'   and bootstrap sample information.
#'
#' @examples
#' \dontrun{
#' x <- rbfr(100, alpha = 2, beta = 3)
#'
#' plot_method(
#'   data = x,
#'   alpha = 2,
#'   beta = 3,
#'   B = 1000,
#'   seed = 123
#' )
#' }
#'
#' @method plot bfr
#' @export
plot.bfr <- function(x, B = 1000, conf = c(0.90, 0.95), seed = NULL, show_observed = TRUE, ...) {
  
  # Extract data and parameters from the fitted bfr list object
  data <- x$data
  alpha <- x$estimate["alpha"] 
  beta <- x$estimate["beta"]
  
  ellipse_colours <- c("steelblue", "firebrick")
  # Input validation

  if (!is.numeric(data)) {
    stop("'data' must be a numeric vector.")
  }

  if (length(data) < 2) {
    stop("'data' must contain at least two observations.")
  }

  if (any(!is.finite(data))) {
    stop("'data' must contain only finite values.")
  }

  if (any(data <= 0 | data >= 1)) {
    stop("'data' must contain observations strictly between 0 and 1.")
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1 ||
      !is.finite(alpha) ||
      alpha <= 0) {
    stop("'alpha' must be a single positive number.")
  }

  if (!is.numeric(beta) ||
      length(beta) != 1 ||
      !is.finite(beta) ||
      beta <= 0) {
    stop("'beta' must be a single positive number.")
  }

  if (!is.numeric(B) ||
      length(B) != 1 ||
      !is.finite(B) ||
      B < 2 ||
      B != floor(B)) {
    stop("'B' must be an integer greater than or equal to 2.")
  }

  if (!is.numeric(conf) ||
      length(conf) < 1 ||
      any(!is.finite(conf)) ||
      any(conf <= 0 | conf >= 1)) {
    stop("'conf' must contain values strictly between 0 and 1.")
  }

  if (!is.null(seed)) {

    if (!is.numeric(seed) ||
        length(seed) != 1 ||
        !is.finite(seed)) {
      stop("'seed' must be NULL or a single finite numeric value.")
    }

    set.seed(seed)
  }

  # Remove duplicate confidence levels
  conf <- unique(conf)

  # Generate parametric bootstrap samples

  n <- length(data)

  bootstrap_stats <- matrix(
    NA_real_,
    nrow = B,
    ncol = 2
  )

  colnames(bootstrap_stats) <- c(
    "SK_G",
    "KU_M"
  )

  for (i in seq_len(B)) {

    # Generate a bootstrap sample from the fitted BFR distribution.
    #
    # rbfr() is defined in core_functions.R.
    bootstrap_sample <- rbfr(
      n = n,
      alpha = alpha,
      beta = beta
    )

    # Calculate Galton skewness.
    bootstrap_stats[i, "SK_G"] <- tryCatch(
      galton_skewness(
        bootstrap_sample,
        na.rm = TRUE
      ),
      error = function(e) NA_real_
    )

    # Calculate Moors kurtosis.
    bootstrap_stats[i, "KU_M"] <- tryCatch(
      moors_kurtosis(
        bootstrap_sample,
        na.rm = TRUE
      ),
      error = function(e) NA_real_
    )
  }

  bootstrap_stats <- as.data.frame(bootstrap_stats)

  # Remove unsuccessful replications.
  bootstrap_stats <- bootstrap_stats[
    is.finite(bootstrap_stats$SK_G) &
      is.finite(bootstrap_stats$KU_M),
    ,
    drop = FALSE
  ]

  if (nrow(bootstrap_stats) < 3) {
    stop(
      "Fewer than three valid bootstrap replications were obtained. ",
      "Try increasing 'B'."
    )
  }
 
  # Calculate statistics from the original data
 
  observed_sk <- galton_skewness(
    data,
    na.rm = TRUE
  )

  observed_ku <- moors_kurtosis(
    data,
    na.rm = TRUE
  )

  # Calculate bootstrap centre
 
  bootstrap_mean <- c(
    mean(bootstrap_stats$SK_G),
    mean(bootstrap_stats$KU_M)
  )

  names(bootstrap_mean) <- c(
    "SK_G",
    "KU_M"
  )
 
  # Calculate covariance matrix

  covariance_matrix <- stats::cov(
    bootstrap_stats[, c("SK_G", "KU_M")]
  )

  covariance_valid <- (
    all(is.finite(covariance_matrix)) &&
    det(covariance_matrix) > 0
  )

  # Determine plotting limits

  plot_x <- bootstrap_stats$SK_G
  plot_y <- bootstrap_stats$KU_M

  if (show_observed &&
      is.finite(observed_sk) &&
      is.finite(observed_ku)) {

    plot_x <- c(plot_x, observed_sk)
    plot_y <- c(plot_y, observed_ku)
  }

  x_range <- range(plot_x, finite = TRUE)
  y_range <- range(plot_y, finite = TRUE)

  x_padding <- diff(x_range) * 0.10
  y_padding <- diff(y_range) * 0.10

  # Prevent zero-width plotting ranges.

  if (x_padding == 0) {
    x_padding <- 0.1
  }

  if (y_padding == 0) {
    y_padding <- 0.1
  }

  x_limits <- c(
    x_range[1] - x_padding,
    x_range[2] + x_padding
  )

  y_limits <- c(
    y_range[1] - y_padding,
    y_range[2] + y_padding
  ) 

  # Create the SK_G x KU_M plot
 
  graphics::plot(
    bootstrap_stats$SK_G,
    bootstrap_stats$KU_M,
    xlim = x_limits,
    ylim = y_limits,
    pch = 16,
    cex = 0.65,
    xlab = expression(SK[G]),
    ylab = expression(KU[M]),
    main = "BFR Bootstrap SK_G x KU_M Plot"
  )
 
# Draw confidence ellipses

if (covariance_valid) {

  eigen_result <- eigen(
    covariance_matrix,
    symmetric = TRUE
  )

  eigen_values <- eigen_result$values
  eigen_vectors <- eigen_result$vectors

  theta <- seq(
    0,
    2 * pi,
    length.out = 361
  )

  unit_circle <- rbind(
    cos(theta),
    sin(theta)
  )

  # Sort confidence levels from lowest to highest
  sorted_conf <- sort(conf)

  # Make sure there is one colour for each ellipse
  if (length(ellipse_colours) < length(sorted_conf)) {
    ellipse_colours <- rep(
      ellipse_colours,
      length.out = length(sorted_conf)
    )
  }

  # Draw each confidence ellipse
  for (i in seq_along(sorted_conf)) {

    level <- sorted_conf[i]

    # Radius corresponding to the confidence level
    # of a bivariate normal distribution.
    radius <- sqrt(
      stats::qchisq(
        level,
        df = 2
      )
    )

    # Calculate ellipse
    ellipse <- eigen_vectors %*%
      diag(sqrt(pmax(eigen_values, 0))) %*%
      unit_circle * radius

    # Move ellipse to bootstrap centre
    ellipse[1, ] <- (
      ellipse[1, ] +
        bootstrap_mean["SK_G"]
    )

    ellipse[2, ] <- (
      ellipse[2, ] +
        bootstrap_mean["KU_M"]
    )

    # Draw coloured ellipse
    graphics::lines(
      ellipse[1, ],
      ellipse[2, ],
      col = ellipse_colours[i],
      lwd = 3
    )
  }

} else {

  warning(
    "The bootstrap covariance matrix is singular or invalid. ",
    "Confidence ellipses could not be drawn."
  )
}
 
  # Plot bootstrap centre

  graphics::points(
    bootstrap_mean["SK_G"],
    bootstrap_mean["KU_M"],
    pch = 4,
    cex = 1.5,
    lwd = 2
  ) 
  # Plot observed statistics

  if (show_observed &&
      is.finite(observed_sk) &&
      is.finite(observed_ku)) {

    graphics::points(
      observed_sk,
      observed_ku,
      pch = 8,
      cex = 1.5,
      lwd = 2
    )
  }

# Legend

legend_labels <- c(
  "Bootstrap replications",
  "Bootstrap centre"
)

legend_pch <- c(
  16,
  4
)

legend_lty <- c(
  NA,
  NA
)

legend_col <- c(
  "black",
  "black"
)

if (show_observed &&
    is.finite(observed_sk) &&
    is.finite(observed_ku)) {

  legend_labels <- c(
    legend_labels,
    "Observed statistic"
  )

  legend_pch <- c(
    legend_pch,
    8
  )

  legend_lty <- c(
    legend_lty,
    NA
  )

  legend_col <- c(
    legend_col,
    "black"
  )
}

sorted_conf <- sort(conf)

for (i in seq_along(sorted_conf)) {

  level <- sorted_conf[i]

  legend_labels <- c(
    legend_labels,
    paste0(
      level * 100,
      "% confidence ellipse"
    )
  )

  legend_pch <- c(
    legend_pch,
    NA
  )

  legend_lty <- c(
    legend_lty,
    1
  )

  legend_col <- c(
    legend_col,
    ellipse_colours[i]
  )
}

graphics::legend(
  "topright",
  legend = legend_labels,
  col = legend_col,
  pch = legend_pch,
  lty = legend_lty,
  lwd = 3,
  bty = "n"
)
 
  # Return results invisibly

  invisible(
    list(
      bootstrap = bootstrap_stats,
      observed = c(
        SK_G = observed_sk,
        KU_M = observed_ku
      ),
      centre = bootstrap_mean,
      confidence_levels = conf,
      alpha = alpha,
      beta = beta,
      sample_size = n,
      requested_bootstrap_replicates = B,
      valid_bootstrap_replicates = nrow(bootstrap_stats)
    )
  )
}