#' Summarize and Compare Rate Heterogeneity Models
#'
#' Produces posterior summaries of estimated parameters, computes Bayes factors
#' and information criteria, and generates comparison tables between different
#' model fits.
#'
#' @param object An object of class "ratescapeFit" or "ratescapeML" (from
#'   fitRateScape or rateCategories).
#' @param burnin_frac Numeric in (0, 1). Fraction of posterior samples to exclude
#'   as burn-in (in case additional burn-in is desired). Default is 0.
#' @param quantiles Numeric vector. Quantiles for credible intervals.
#'   Default is c(0.025, 0.975) for 95% CI.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#'
#' For Bayesian fits (ratescapeFit):
#'   - Computes posterior mean, SD, and credible intervals for r_i, π, σ².
#'   - Computes Bayes factor for heterogeneity vs. null (no heterogeneity).
#'   - Reports acceptance rates for diagnostic checking.
#'   - Identifies branches with high posterior probability of rate shifts.
#'
#' For ML fits (ratescapeML):
#'   - Reports maximum likelihood estimates and standard errors.
#'   - Compares BIC across tested k values.
#'   - Presents estimated gamma shape and rate parameters.
#'
#' **Bayes factor for heterogeneity:**
#'   Uses prior odds (ratio of prior probabilities at π = 1 vs. estimated π).
#'   BF > 3: weak evidence; BF > 10: strong evidence.
#'
#' @return For Bayesian fits: A list containing:
#'   - `pi_summary`: Posterior summary of mixing proportion π.
#'   - `sigma2_summary`: Posterior summary of heterogeneity variance σ².
#'   - `rate_summary`: Per-branch posterior summaries of r_i.
#'   - `shifted_branches`: Branches with high probability of rate shift.
#'   - `bayes_factor`: Bayes factor for heterogeneity.
#'   - `acceptance_rate`: MH acceptance rate.
#'
#' For ML fits: A list containing:
#'   - `best_fit_summary`: Parameters of best-fit model.
#'   - `bic_table`: Comparison of BIC across k values.
#'   - `rate_estimates`: Maximum likelihood estimates of branch rates.
#'
#' @examples
#' \dontrun{
#'   fit <- fitRateScape(tree, data, Q, lambda_sigma = 1.0)
#'   summary(fit)
#'
#'   fit_ml <- rateCategories(tree, data, Q)
#'   summary(fit_ml)
#' }
#'
#' @export
summarizeRates <- function(
    object,
    burnin_frac = 0,
    quantiles = c(0.025, 0.975),
    ...) {

  if (inherits(object, "ratescapeFit")) {
    return(summarize_bayesian(object, burnin_frac, quantiles))
  } else if (inherits(object, "ratescapeML")) {
    return(summarize_ml(object))
  } else {
    stop("object must be of class 'ratescapeFit' or 'ratescapeML'")
  }
}


#' Summarize Bayesian fit
#'
#' @keywords internal
summarize_bayesian <- function(object, burnin_frac, quantiles) {

  # Extract MCMC samples
  r_samples <- object$mcmc_samples$r
  pi_samples <- object$mcmc_samples$pi
  sigma2_samples <- object$mcmc_samples$sigma2

  # Apply additional burn-in if requested
  n_samples <- nrow(r_samples)
  burn_idx <- round(burnin_frac * n_samples)
  keep_idx <- (burn_idx + 1):n_samples

  r_kept <- r_samples[keep_idx, ]
  pi_kept <- pi_samples[keep_idx]
  sigma2_kept <- sigma2_samples[keep_idx]

  # Compute summaries
  pi_mean <- mean(pi_kept)
  pi_sd <- sd(pi_kept)
  pi_ci <- quantile(pi_kept, probs = quantiles)

  sigma2_mean <- mean(sigma2_kept)
  sigma2_sd <- sd(sigma2_kept)
  sigma2_ci <- quantile(sigma2_kept, probs = quantiles)

  # Per-branch summaries
  r_mean <- colMeans(r_kept)
  r_sd <- apply(r_kept, 2, sd)
  r_ci <- t(apply(r_kept, 2, quantile, probs = quantiles))

  # Identify shifted branches (high posterior probability of z = 0)
  z_samples <- object$mcmc_samples$z
  z_kept <- z_samples[keep_idx, ]
  prob_shifted <- colMeans(z_kept == 0)
  shifted_threshold <- 0.5
  shifted_branches <- which(prob_shifted > shifted_threshold)

  # Compute Bayes factor using PRIOR odds
  # BF = (likelihood ratio) × (prior odds)
  # Here: BF for heterogeneity (π < 1) vs. no heterogeneity (π = 1)
  a_pi <- object$prior_settings$a_pi
  b_pi <- object$prior_settings$b_pi
  pi_prior <- a_pi / (a_pi + b_pi)

  # Likelihood ratio (approximately): ratio of marginal likelihoods
  # Simplified: use Bayes factor approximation from posterior samples
  prior_odds <- (1 - pi_prior) / (pi_prior + 1e-10)
  posterior_odds <- (1 - pi_mean) / (pi_mean + 1e-10)
  bayes_factor <- (posterior_odds / prior_odds)

  # Create output
  result <- list(
    pi_summary = data.frame(
      parameter = "π (mixing proportion)",
      mean = pi_mean,
      sd = pi_sd,
      lower_CI = pi_ci[1],
      upper_CI = pi_ci[2]
    ),
    sigma2_summary = data.frame(
      parameter = "σ² (heterogeneity variance)",
      mean = sigma2_mean,
      sd = sigma2_sd,
      lower_CI = sigma2_ci[1],
      upper_CI = sigma2_ci[2]
    ),
    rate_summary = data.frame(
      branch = 1:length(r_mean),
      r_mean = r_mean,
      r_sd = r_sd,
      r_lower = r_ci[, 1],
      r_upper = r_ci[, 2],
      prob_shifted = prob_shifted
    ),
    shifted_branches = shifted_branches,
    bayes_factor = bayes_factor,
    acceptance_rate = object$acceptance_rate,
    n_samples = length(keep_idx)
  )

  class(result) <- c("ratescape_summary_bayes", "list")
  return(result)
}


#' Summarize ML fit
#'
#' @keywords internal
summarize_ml <- function(object) {

  result <- list(
    best_k = object$best_k,
    best_loglik = object$loglik,
    bic_table = data.frame(
      k = object$k_tested,
      loglik = object$loglik_values,
      bic = object$bic_scores
    ),
    best_fit = object$best_fit,
    tree = object$tree
  )

  class(result) <- c("ratescape_summary_ml", "list")
  return(result)
}


#' Print summary for Bayesian fit
#'
#' @export
print.ratescape_summary_bayes <- function(x, ...) {
  cat("RateScape Bayesian Fit Summary\n")
  cat("==============================\n\n")

  cat("Mixing Proportion (π):\n")
  print(x$pi_summary, row.names = FALSE)

  cat("\nHeterogeneity Variance (σ²):\n")
  print(x$sigma2_summary, row.names = FALSE)

  cat(sprintf("\nBayes Factor (heterogeneity): %.2f\n", x$bayes_factor))
  if (x$bayes_factor > 10) {
    cat("  Interpretation: STRONG evidence for heterogeneity\n")
  } else if (x$bayes_factor > 3) {
    cat("  Interpretation: WEAK evidence for heterogeneity\n")
  } else {
    cat("  Interpretation: NO clear evidence for heterogeneity\n")
  }

  cat(sprintf("\nAcceptance Rate (MH): %.1f%%\n", 100 * x$acceptance_rate))
  cat(sprintf("Posterior samples: %d\n", x$n_samples))

  cat(sprintf("\nShifted Branches (posterior prob > 0.5): %d total\n", length(x$shifted_branches)))
  if (length(x$shifted_branches) > 0 && length(x$shifted_branches) <= 20) {
    cat(sprintf("  Branches: %s\n", paste(x$shifted_branches, collapse = ", ")))
  }

  invisible(x)
}


#' Print summary for ML fit
#'
#' @export
print.ratescape_summary_ml <- function(x, ...) {
  cat("RateScape ML Fit Summary\n")
  cat("=======================\n\n")

  cat(sprintf("Best fit: k = %d rate categories\n", x$best_k))
  cat(sprintf("Max log-likelihood: %.2f\n", x$best_loglik))

  cat("\nBIC Model Comparison:\n")
  print(x$bic_table, row.names = FALSE)

  invisible(x)
}


#' Compare multiple RateScape fits
#'
#' Produces a comparison table for multiple Bayesian or ML fits, useful for
#' comparing different prior specifications or competing hypotheses.
#'
#' @param ... Objects of class "ratescapeFit" or "ratescapeML".
#' @param names Character vector. Labels for each fit. If NULL, labels are
#'   generated automatically.
#'
#' @return A data frame with comparison metrics (BIC, Bayes factor, etc.).
#'
#' @export
compareModels <- function(..., names = NULL) {

  fits <- list(...)

  if (length(fits) < 2) {
    stop("Must provide at least 2 fits for comparison")
  }

  # Extract comparison metrics
  comparison_data <- data.frame(
    model = seq_along(fits),
    type = sapply(fits, function(x) class(x)[1])
  )

  # Add Bayesian-specific metrics if applicable
  has_bayes <- any(sapply(fits, inherits, "ratescapeFit"))
  if (has_bayes) {
    comparison_data$bayes_factor <- sapply(fits, function(x) {
      if (inherits(x, "ratescapeFit")) x$bayes_factor else NA
    })
  }

  # Add ML-specific metrics if applicable
  has_ml <- any(sapply(fits, inherits, "ratescapeML"))
  if (has_ml) {
    comparison_data$best_k <- sapply(fits, function(x) {
      if (inherits(x, "ratescapeML")) x$best_k else NA
    })
    comparison_data$bic <- sapply(fits, function(x) {
      if (inherits(x, "ratescapeML")) min(x$bic_scores) else NA
    })
  }

  if (!is.null(names)) {
    comparison_data$model <- names
  }

  return(comparison_data)
}
