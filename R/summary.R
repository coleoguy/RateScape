#' Summarize Branch-Specific Rate Estimates
#'
#' Extracts a data frame of branch-specific rate estimates, credible intervals,
#' posterior shift probabilities, and Bayes factors from a fitted RateScape model.
#'
#' @param object A fitted ratescape object
#' @param prob Credible interval probability. Default 0.95.
#' @param ... Additional arguments (ignored)
#' @return A data frame with one row per edge and columns:
#'   edge, parent, child, rate_mean, rate_median, rate_lower, rate_upper,
#'   shift_prob, bayes_factor
#' @export
summarizeRates <- function(object, prob = 0.95, ...) {
  if (!inherits(object, "ratescape")) stop("object must be a ratescape object")

  tree <- object$tree
  nedge <- ape::Nedge(tree)

  if (inherits(object, "ratescape_ml")) {
    # ML approach
    df <- data.frame(
      edge = 1:nedge,
      parent = tree$edge[, 1],
      child = tree$edge[, 2],
      rate_mean = object$branch_rates,
      rate_category_probs = I(lapply(1:nedge, function(e) {
        setNames(object$branch_posteriors[e, ],
                 paste0("cat", 1:object$ncat))
      })),
      stringsAsFactors = FALSE
    )

    # Add tip labels for terminal edges
    ntip <- ape::Ntip(tree)
    df$tip_label <- ifelse(df$child <= ntip,
                           tree$tip.label[df$child], NA)

  } else {
    # Bayesian approach
    alpha <- (1 - prob) / 2

    rate_lower <- apply(object$rate_samples, 2,
                        stats::quantile, probs = alpha)
    rate_upper <- apply(object$rate_samples, 2,
                        stats::quantile, probs = 1 - alpha)

    df <- data.frame(
      edge = 1:nedge,
      parent = tree$edge[, 1],
      child = tree$edge[, 2],
      rate_mean = object$rate_means,
      rate_median = object$rate_medians,
      rate_lower = rate_lower,
      rate_upper = rate_upper,
      shift_prob = object$shift_probs,
      bayes_factor = object$bayes_factors,
      stringsAsFactors = FALSE
    )

    # Add tip labels for terminal edges
    ntip <- ape::Ntip(tree)
    df$tip_label <- ifelse(df$child <= ntip,
                           tree$tip.label[df$child], NA)
  }

  # Sort by Bayes factor (or rate) descending
  if ("bayes_factor" %in% names(df)) {
    df <- df[order(-df$bayes_factor), ]
  } else {
    df <- df[order(-abs(log(df$rate_mean))), ]
  }

  return(df)
}


#' Compare Homogeneous vs. Rate-Variable Models
#'
#' Computes support for rate heterogeneity by comparing the homogeneous
#' (all rates = 1) model against the fitted rate-variable model.
#'
#' @param object A fitted ratescape object
#' @param ... Additional arguments (ignored)
#' @return A list with model comparison statistics
#' @export
compareModels <- function(object, ...) {
  if (!inherits(object, "ratescape")) stop("object must be a ratescape object")

  tree <- object$tree
  nedge <- ape::Nedge(tree)

  if (inherits(object, "ratescape_ml")) {
    # ML comparison
    # Null model: homogeneous rates
    tip_states <- as.integer(object$data)
    names(tip_states) <- names(object$data)
    rates_null <- rep(1.0, nedge)
    loglik_null <- pruning_likelihood_scaled(tree, tip_states, object$Q,
                                             rates_null, "equal")

    result <- list(
      method = "ML",
      loglik_null = loglik_null,
      loglik_alt = object$loglik,
      delta_AIC = (-2 * loglik_null + 2 * 0) - object$AIC,
      LRT_stat = 2 * (object$loglik - loglik_null),
      n_params_null = 0,
      n_params_alt = 1  # alpha
    )

  } else {
    # Bayesian comparison via spike-and-slab
    # The posterior on pi directly tells us about rate heterogeneity
    # If pi is near 1, homogeneous model is supported

    # Posterior probability of any heterogeneity
    p_any_shift <- mean(apply(object$z_samples, 1, function(z) any(z == 1)))

    # Savage-Dickey approximation for BF of heterogeneous vs homogeneous
    # BF = p(pi=1 | prior) / p(pi=1 | posterior)
    # Under Beta(a,b) prior, p(pi=1) = 0 (continuous), so use approximate
    # near-boundary density

    result <- list(
      method = "Bayesian",
      pi_posterior_mean = object$pi_mean,
      pi_95CI = stats::quantile(object$pi_samples, c(0.025, 0.975)),
      prob_any_shift = p_any_shift,
      n_shifted_mean = mean(rowSums(object$z_samples)),
      sigma2_posterior_mean = object$sigma2_mean,
      mean_loglik = mean(object$llik_samples),
      loglik_95CI = stats::quantile(object$llik_samples, c(0.025, 0.975))
    )
  }

  class(result) <- "ratescape_comparison"
  return(result)
}


#' Print model comparison results
#' @param x A ratescape_comparison object
#' @param ... Additional arguments (ignored)
#' @export
print.ratescape_comparison <- function(x, ...) {
  cat("RateScape Model Comparison\n")
  cat("==========================\n\n")

  if (x$method == "ML") {
    cat("Method: Maximum Likelihood (discretized gamma)\n")
    cat(sprintf("  Homogeneous logLik: %.2f\n", x$loglik_null))
    cat(sprintf("  Rate-variable logLik: %.2f\n", x$loglik_alt))
    cat(sprintf("  LRT statistic: %.2f\n", x$LRT_stat))
    cat(sprintf("  Delta AIC (null - alt): %.2f\n", x$delta_AIC))
    if (x$delta_AIC > 0) {
      cat("  => Rate-variable model preferred\n")
    } else {
      cat("  => Homogeneous model preferred\n")
    }
  } else {
    cat("Method: Bayesian (spike-and-slab)\n")
    cat(sprintf("  Posterior P(background) [pi]: %.3f (%.3f, %.3f)\n",
                x$pi_posterior_mean, x$pi_95CI[1], x$pi_95CI[2]))
    cat(sprintf("  P(any rate shift): %.3f\n", x$prob_any_shift))
    cat(sprintf("  Mean number of shifted branches: %.1f\n", x$n_shifted_mean))
    cat(sprintf("  Posterior sigma2: %.3f\n", x$sigma2_posterior_mean))
    cat(sprintf("  Mean logLik: %.2f (%.2f, %.2f)\n",
                x$mean_loglik, x$loglik_95CI[1], x$loglik_95CI[2]))
  }

  invisible(x)
}
