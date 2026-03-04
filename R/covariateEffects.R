#' Summarize Covariate Effects on Branch-Specific Rates
#'
#' Extracts and summarizes the posterior distribution of regression coefficients
#' from a covariate-driven RateScape fit. Reports posterior means, credible
#' intervals, and the probability that each coefficient is non-zero (positive
#' or negative).
#'
#' @param fit An object of class "ratescapeFit" with covariates.
#' @param burnin_frac Additional burn-in fraction to discard. Default 0.
#' @param prob Credible interval probability. Default 0.95.
#'
#' @details
#'
#' When covariates are included in fitRateScape, the slab component of the
#' spike-and-slab model becomes:
#'
#'   log(r_i) = X_i * beta + epsilon_i,  epsilon_i ~ N(0, sigma^2)
#'
#' where X is the standardized covariate matrix and beta are the regression
#' coefficients. Positive beta_j means that higher values of covariate j are
#' associated with faster rates of discrete character evolution on that branch.
#'
#' **Coefficients are on the standardized scale.** To convert to the original
#' covariate scale, divide by the standard deviation of the covariate and
#' multiply by the original SD:
#'
#'   beta_original = beta_standardized / sd(covariate)
#'
#' The function automatically provides both scales.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{covariate}{Covariate name}
#'     \item{mean_std}{Posterior mean on standardized scale}
#'     \item{sd_std}{Posterior SD on standardized scale}
#'     \item{lower_std}{Lower CI bound (standardized)}
#'     \item{upper_std}{Upper CI bound (standardized)}
#'     \item{mean_orig}{Posterior mean on original scale}
#'     \item{prob_positive}{Posterior probability that beta > 0}
#'     \item{prob_negative}{Posterior probability that beta < 0}
#'     \item{significant}{Whether 0 is outside the credible interval}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Fit with body size as a covariate
#'   body_size <- log(species_masses[match(tree$tip.label, names(species_masses))])
#'   # Branch covariates: average of parent and child tip values (or custom)
#'   edge_covars <- matrix(body_size_by_edge, ncol = 1)
#'   colnames(edge_covars) <- "log_body_size"
#'
#'   fit <- fitRateScape(tree, data, Q, lambda_sigma = 1.0,
#'                        covariates = edge_covars)
#'   covariateEffects(fit)
#' }
#'
#' @export
covariateEffects <- function(fit, burnin_frac = 0, prob = 0.95) {

  if (!inherits(fit, "ratescapeFit")) {
    stop("fit must be of class 'ratescapeFit'")
  }
  if (!fit$use_covariates) {
    stop("This fit does not include covariates. Re-run fitRateScape with covariates argument.")
  }

  beta_samples <- fit$mcmc_samples$beta
  n_samples <- nrow(beta_samples)
  burn_idx <- round(burnin_frac * n_samples)
  keep_idx <- (burn_idx + 1):n_samples
  beta_kept <- beta_samples[keep_idx, , drop = FALSE]

  covar_names <- fit$covariate_info$names
  raw_sds <- fit$covariate_info$raw_sds
  n_covars <- length(covar_names)

  alpha <- (1 - prob) / 2

  results <- data.frame(
    covariate = covar_names,
    mean_std = numeric(n_covars),
    sd_std = numeric(n_covars),
    lower_std = numeric(n_covars),
    upper_std = numeric(n_covars),
    mean_orig = numeric(n_covars),
    prob_positive = numeric(n_covars),
    prob_negative = numeric(n_covars),
    significant = logical(n_covars),
    stringsAsFactors = FALSE
  )

  for (j in 1:n_covars) {
    samples_j <- beta_kept[, j]
    results$mean_std[j] <- mean(samples_j)
    results$sd_std[j] <- sd(samples_j)
    results$lower_std[j] <- quantile(samples_j, alpha)
    results$upper_std[j] <- quantile(samples_j, 1 - alpha)
    results$mean_orig[j] <- mean(samples_j) / raw_sds[j]
    results$prob_positive[j] <- mean(samples_j > 0)
    results$prob_negative[j] <- mean(samples_j < 0)
    results$significant[j] <- results$lower_std[j] > 0 || results$upper_std[j] < 0
  }

  cat("RateScape Covariate Effects\n")
  cat("===========================\n\n")
  cat(sprintf("Posterior samples: %d (after %.0f%% burn-in)\n",
              length(keep_idx), burnin_frac * 100))
  cat(sprintf("Credible interval: %.0f%%\n\n", prob * 100))

  for (j in 1:n_covars) {
    sig_mark <- if (results$significant[j]) " *" else ""
    cat(sprintf("  %s%s:\n", results$covariate[j], sig_mark))
    cat(sprintf("    Standardized: %.3f [%.3f, %.3f]\n",
                results$mean_std[j], results$lower_std[j], results$upper_std[j]))
    cat(sprintf("    Original scale: %.4f per unit\n", results$mean_orig[j]))
    cat(sprintf("    P(beta > 0) = %.3f, P(beta < 0) = %.3f\n",
                results$prob_positive[j], results$prob_negative[j]))
  }

  if (any(results$significant)) {
    cat("\n  * Coefficient significantly different from zero\n")
  }

  invisible(results)
}


#' Create branch-level covariates from tip-level data
#'
#' Helper function that converts tip-level continuous traits into edge-level
#' covariates suitable for the covariates argument of fitRateScape. For each
#' edge, the covariate is computed as the mean of the values at the tips
#' descending from that edge's child node.
#'
#' @param tree phylo object.
#' @param tip_values Named numeric vector of values at tips. Names must match
#'   tree tip labels.
#' @param method How to compute edge-level values from tip descendant values.
#'   "mean" (default), "median", "max", "min", or "contrast" (absolute value
#'   of the difference between the two child-clade means, only for internal
#'   edges).
#'
#' @return A numeric vector with one value per edge, in the same order as
#'   tree$edge.
#'
#' @examples
#' \dontrun{
#'   # Body size by edge
#'   body_sizes <- c(sp1 = 2.3, sp2 = 4.1, sp3 = 1.8, ...)
#'   edge_sizes <- tipToEdgeCovariate(tree, body_sizes)
#'   covars <- matrix(edge_sizes, ncol = 1)
#'   colnames(covars) <- "body_size"
#' }
#'
#' @export
tipToEdgeCovariate <- function(tree, tip_values, method = "mean") {
  method <- match.arg(method, c("mean", "median", "max", "min", "contrast"))

  ntips <- length(tree$tip.label)
  nedges <- nrow(tree$edge)

  if (length(tip_values) != ntips) {
    stop("tip_values must have length equal to number of tips")
  }

  # Match to tree order
  if (!is.null(names(tip_values))) {
    mi <- match(tree$tip.label, names(tip_values))
    if (any(is.na(mi))) stop("Some tip labels not found in tip_values names")
    tip_values <- tip_values[mi]
  }

  # For each edge, get descendant tips
  edge_values <- numeric(nedges)
  for (e in 1:nedges) {
    child <- tree$edge[e, 2]
    if (child <= ntips) {
      # Terminal edge: just use the tip value
      edge_values[e] <- tip_values[child]
    } else {
      # Internal edge: aggregate descendant tip values
      desc <- get_descendants(child, tree, ntips)
      vals <- tip_values[desc]
      edge_values[e] <- switch(method,
        "mean" = mean(vals),
        "median" = median(vals),
        "max" = max(vals),
        "min" = min(vals),
        "contrast" = {
          # Difference between the two child clades
          children_edges <- which(tree$edge[, 1] == child)
          if (length(children_edges) == 2) {
            c1 <- tree$edge[children_edges[1], 2]
            c2 <- tree$edge[children_edges[2], 2]
            d1 <- if (c1 <= ntips) c1 else get_descendants(c1, tree, ntips)
            d2 <- if (c2 <= ntips) c2 else get_descendants(c2, tree, ntips)
            abs(mean(tip_values[d1]) - mean(tip_values[d2]))
          } else {
            0
          }
        }
      )
    }
  }

  edge_values
}
