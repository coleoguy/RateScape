#' Plot Rate-Painted Phylogenetic Trees
#'
#' Produces a phylogenetic tree visualization with branches colored by their
#' estimated rate scalars, creating a "painted phylogeny" showing lineage-specific
#' rate heterogeneity.
#'
#' @param x An object of class "ratescapeFit" or "ratescapeML" (from fitRateScape
#'   or rateCategories).
#' @param y Unused (for compatibility with plot generic).
#' @param quantile_threshold Numeric. Posterior probability threshold for identifying
#'   rate-shifted branches (default 0.5). Only used for Bayesian fits.
#' @param color_palette Character. Color scheme: "diverging" (blue-grey-red, default),
#'   "sequential" (white-orange), or "viridis".
#' @param label_shifted Logical. If TRUE, label branches with rate shifts.
#'   Default is FALSE.
#' @param ... Additional arguments passed to ape::plot.phylo.
#'
#' @details
#'
#' **For Bayesian fits:**
#'   - Branch colors represent posterior mean r_i estimates.
#'   - Blue branches: r_i < 1 (slower evolution).
#'   - Red branches: r_i > 1 (faster evolution).
#'   - Grey branches: posterior probability of shift < threshold (consistent with background rate).
#'
#' **For ML fits:**
#'   - Branch colors represent maximum likelihood estimates.
#'   - Color intensity reflects the magnitude of rate deviation from 1.
#'
#' **Diverging palette:**
#'   Centered at 1 (white/grey), with blue for r < 1 and red for r > 1.
#'   This makes the background rate (r = 1) clearly visible.
#'
#' @return Invisibly returns the plot object (for chaining with other graphics commands).
#'
#' @examples
#' \dontrun{
#'   fit <- fitRateScape(tree, data, Q, lambda_sigma = 1.0)
#'   plot(fit, main = "Rate-Painted Tree")
#'
#'   fit_ml <- rateCategories(tree, data, Q)
#'   plot(fit_ml, color_palette = "sequential")
#' }
#'
#' @export
plotRateTree <- function(
    x,
    y = NULL,
    quantile_threshold = 0.5,
    color_palette = "diverging",
    label_shifted = FALSE,
    ...) {

  if (inherits(x, "ratescapeFit")) {
    return(plot_bayesian_tree(x, quantile_threshold, color_palette, label_shifted, ...))
  } else if (inherits(x, "ratescapeML")) {
    return(plot_ml_tree(x, color_palette, ...))
  } else {
    stop("x must be of class 'ratescapeFit' or 'ratescapeML'")
  }
}


#' Plot Bayesian painted tree
#'
#' @keywords internal
plot_bayesian_tree <- function(x, quantile_threshold, color_palette, label_shifted, ...) {

  tree <- x$tree
  r_samples <- x$mcmc_samples$r
  z_samples <- x$mcmc_samples$z

  # Compute posterior summaries
  r_mean <- colMeans(r_samples)
  prob_shifted <- colMeans(z_samples == 0)

  # For display: mark branches with low probability of shift as background (r ≈ 1)
  r_display <- r_mean
  r_display[prob_shifted < quantile_threshold] <- 1.0

  # Create color palette
  colors <- get_rate_colors(r_display, palette = color_palette)

  # Plot tree with colored edges
  plot(tree, edge.color = colors, edge.width = 1.5, ...)

  # Add legend
  add_rate_legend(r_display, palette = color_palette)

  invisible(NULL)
}


#' Plot ML painted tree
#'
#' @keywords internal
plot_ml_tree <- function(x, color_palette, ...) {

  tree <- x$tree
  best_fit <- x$best_fit

  # Extract rate estimates from best fit
  r_estimates <- best_fit$rates

  # Create color palette
  colors <- get_rate_colors(r_estimates, palette = color_palette)

  # Plot tree
  plot(tree, edge.color = colors, edge.width = 1.5, ...)

  # Add legend
  add_rate_legend(r_estimates, palette = color_palette)

  invisible(NULL)
}


#' Generate colors for rate values
#'
#' @keywords internal
get_rate_colors <- function(rate_values, palette = "diverging") {

  n_colors <- 100
  r_range <- range(rate_values, na.rm = TRUE)

  if (palette == "diverging") {
    # Blue-grey-red centered at 1
    colors_palette <- colorspace::diverging_hcl(n_colors, palette = "RdBu")
  } else if (palette == "sequential") {
    # White to orange gradient
    colors_palette <- colorspace::sequential_hcl(n_colors, palette = "OrYe")
  } else if (palette == "viridis") {
    # Viridis-like gradient
    colors_palette <- grDevices::hcl.colors(n_colors, palette = "viridis")
  } else {
    stop("Unknown palette: ", palette)
  }

  # Map rate values to color indices
  # Normalize rates to [0, 1] centered at 1
  if (palette == "diverging") {
    # Symmetric around 1: rates from 1/max to max map to [0, 1]
    r_max <- max(abs(log(rate_values)))
    color_idx <- pmin(n_colors, pmax(1, round(50 + 50 * log(rate_values) / r_max)))
  } else {
    # Sequential: map from min to max
    color_idx <- pmin(n_colors, pmax(1, round((rate_values - r_range[1]) /
                                               (r_range[2] - r_range[1]) * n_colors)))
  }

  return(colors_palette[color_idx])
}


#' Add legend to rate tree plot
#'
#' @keywords internal
add_rate_legend <- function(rate_values, palette = "diverging") {

  r_range <- range(rate_values, na.rm = TRUE)

  if (palette == "diverging") {
    legend_text <- sprintf(
      "Rates: blue = slow (r=%.1f), red = fast (r=%.1f), background = r=1",
      r_range[1], r_range[2]
    )
  } else {
    legend_text <- sprintf(
      "Rates: min = %.2f, max = %.2f",
      r_range[1], r_range[2]
    )
  }

  mtext(legend_text, side = 1, line = -1, cex = 0.8, outer = FALSE)
}


#' Generic plot method for RateScape fits
#'
#' @export
plot.ratescapeFit <- function(x, y = NULL, ...) {
  plotRateTree(x, ...)
}


#' Generic plot method for ML fits
#'
#' @export
plot.ratescapeML <- function(x, y = NULL, ...) {
  plotRateTree(x, ...)
}
