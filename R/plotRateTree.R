#' Plot Painted Phylogeny with Branch-Specific Rates
#'
#' Produces a phylogeny with branches colored by their estimated rate scalar,
#' from blue (decelerated) through grey (background) to red (accelerated).
#' Optionally highlights branches with strong support for rate shifts.
#'
#' @param x A fitted ratescape or ratescape_ml object
#' @param type Character: "fan", "phylogram", or "cladogram". Default "phylogram".
#' @param color_by Character: "mean" (posterior mean rate), "median" (posterior
#'   median rate), or "shift_prob" (posterior probability of shift). Default "mean".
#' @param palette A vector of 3 colors for (slow, background, fast).
#'   Default c("blue", "grey80", "red").
#' @param n_colors Number of colors in the gradient. Default 100.
#' @param rate_range Numeric vector of length 2 giving the range of rates to
#'   map to the color scale. Default NULL (auto from data).
#' @param highlight_bf Bayes factor threshold for highlighting shifted branches.
#'   Default NULL (no highlighting). Set to e.g. 10 for strong evidence.
#' @param highlight_lwd Line width multiplier for highlighted branches. Default 3.
#' @param show_legend Logical. Show color legend? Default TRUE.
#' @param cex Tip label size. Default 0.6.
#' @param ... Additional arguments passed to ape::plot.phylo
#' @export
plotRateTree <- function(x, type = "phylogram",
                         color_by = "mean",
                         palette = c("blue", "grey80", "red"),
                         n_colors = 100,
                         rate_range = NULL,
                         highlight_bf = NULL,
                         highlight_lwd = 3,
                         show_legend = TRUE,
                         cex = 0.6, ...) {

  if (!inherits(x, "ratescape")) stop("x must be a ratescape object")

  tree <- x$tree

  # Get rates to color by
  if (inherits(x, "ratescape_ml")) {
    rates <- x$branch_rates
    has_bf <- FALSE
  } else {
    if (color_by == "mean") {
      rates <- x$rate_means
    } else if (color_by == "median") {
      rates <- x$rate_medians
    } else if (color_by == "shift_prob") {
      rates <- x$shift_probs
    } else {
      stop("color_by must be 'mean', 'median', or 'shift_prob'")
    }
    has_bf <- TRUE
  }

  # Build color gradient
  color_func <- grDevices::colorRampPalette(palette)
  colors <- color_func(n_colors)

  # Map rates to colors
  if (is.null(rate_range)) {
    if (color_by == "shift_prob") {
      rate_range <- c(0, 1)
    } else {
      # Symmetric around 1
      max_dev <- max(abs(log(rates)))
      rate_range <- exp(c(-max_dev, max_dev))
    }
  }

  # Map rates to color index
  if (color_by == "shift_prob") {
    # Linear mapping 0 to 1
    color_idx <- round((rates - rate_range[1]) /
                          (rate_range[2] - rate_range[1]) * (n_colors - 1)) + 1
  } else {
    # Log-scale mapping centered at 1
    log_rates <- log(rates)
    log_range <- log(rate_range)
    color_idx <- round((log_rates - log_range[1]) /
                          (log_range[2] - log_range[1]) * (n_colors - 1)) + 1
  }
  color_idx <- pmax(1, pmin(n_colors, color_idx))
  edge_colors <- colors[color_idx]

  # Edge widths
  edge_widths <- rep(2, ape::Nedge(tree))

  # Highlight branches with strong BF

  if (!is.null(highlight_bf) && has_bf) {
    strong <- x$bayes_factors > highlight_bf
    edge_widths[strong] <- highlight_lwd
  }

  # Plot
  ape::plot.phylo(tree, type = type,
                  edge.color = edge_colors,
                  edge.width = edge_widths,
                  cex = cex, ...)

  # Legend
  if (show_legend) {
    if (color_by == "shift_prob") {
      legend_labels <- c("0", "0.5", "1")
      legend_title <- "P(rate shift)"
    } else {
      legend_labels <- c(
        sprintf("%.2f", rate_range[1]),
        "1.00",
        sprintf("%.2f", rate_range[2])
      )
      legend_title <- "Rate scalar"
    }

    # Simple color bar via legend
    legend_colors <- colors[c(1, n_colors / 2, n_colors)]
    graphics::legend("bottomleft",
                     legend = legend_labels,
                     col = legend_colors,
                     lwd = 4,
                     title = legend_title,
                     bty = "n",
                     cex = 0.7)
  }
}
