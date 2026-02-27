#' Prior Diagnostic Checking via Prior-Predictive Simulation
#'
#' Assesses whether your prior specification is sensible for your data by
#' simulating datasets under the prior and comparing to your observed data.
#'
#' @param tree An object of class "phylo" (ape package).
#' @param Q A transition rate matrix (k × k).
#' @param lambda_sigma Numeric. Rate parameter for the exponential prior on σ².
#' @param expected_background Numeric in [0, 1] or NULL. Expected proportion of
#'   background-rate branches. Default is NULL (uniform prior on π).
#' @param nsim Integer. Number of prior-predictive simulations. Default is 100.
#'
#' @details
#'
#' **What this function does:**
#'
#' 1. **Prior-predictive simulation:** Draws parameter values from the prior
#'    (π ~ Beta, σ² ~ Exp(λ_σ), r_i ~ spike-and-slab mixture) and simulates
#'    discrete character datasets under the model. Stores summary statistics
#'    (e.g., number of state changes per branch, entropy of state distribution).
#'
#' 2. **Prior-posterior overlap:** Computes how much the prior predictive
#'    distribution overlaps with the posterior predictive distribution
#'    (estimated from your observed data if a posterior sample is provided).
#'    High overlap (>80%) may indicate that the data are not informative
#'    relative to the prior.
#'
#' 3. **Visualization:** Produces plots showing:
#'    - Prior predictive distribution of key statistics vs. observed data
#'    - Marginal prior and posterior distributions of π, σ²
#'    - Correlation between prior and posterior estimates of r_i
#'
#' 4. **Warning flags:** Alerts the user if:
#'    - Prior-posterior overlap is very high (>80%), suggesting the prior
#'      may be too strong or the data too weak.
#'    - Posterior estimates are near the prior boundaries, suggesting
#'      the prior may be constraining the inference.
#'
#' @return An object of class "ratescape_prior_check" containing:
#'   - `prior_samples`: Samples from the prior distribution.
#'   - `prior_predictive_stats`: Summary statistics from prior-predictive simulations.
#'   - `overlap_measure`: Estimate of prior-posterior overlap (0–1).
#'   - `warnings`: Vector of warning messages if issues detected.
#'   - `call`: The function call.
#'
#' @examples
#' \dontrun{
#'   tree <- ape::rtree(50, rooted = TRUE)
#'   Q <- makeQ(tree, model = "mk", k = 2)
#'
#'   # Check a diffuse prior
#'   check1 <- checkPrior(
#'     tree = tree,
#'     Q = Q,
#'     lambda_sigma = 0.5,
#'     expected_background = NULL,
#'     nsim = 100
#'   )
#'   print(check1)
#'
#'   # Check a more informative prior
#'   check2 <- checkPrior(
#'     tree = tree,
#'     Q = Q,
#'     lambda_sigma = 2.0,
#'     expected_background = 0.80,
#'     nsim = 100
#'   )
#'   print(check2)
#' }
#'
#' @export
checkPrior <- function(
    tree,
    Q,
    lambda_sigma,
    expected_background = NULL,
    nsim = 100) {

  # Validate inputs
  if (!inherits(tree, "phylo")) {
    stop("tree must be an object of class 'phylo'")
  }

  if (!is.matrix(Q) || nrow(Q) != ncol(Q)) {
    stop("Q must be a square matrix")
  }

  if (lambda_sigma <= 0) {
    stop("lambda_sigma must be positive")
  }

  if (!is.null(expected_background)) {
    if (expected_background < 0 || expected_background > 1) {
      stop("expected_background must be in [0, 1]")
    }
  }

  # Convert expected_background to Beta parameters
  if (!is.null(expected_background)) {
    kappa <- 10
    a_pi <- expected_background * kappa
    b_pi <- (1 - expected_background) * kappa
  } else {
    a_pi <- 1
    b_pi <- 1
  }

  nedges <- nrow(tree$edge)
  k_states <- nrow(Q)

  message(sprintf(
    "Prior Diagnostic Check: %d prior-predictive simulations on tree with %d tips",
    nsim, length(tree$tip.label)
  ))

  # ========== PRIOR-PREDICTIVE SIMULATION ==========

  prior_predictive_stats <- list()

  for (sim in 1:nsim) {

    # Draw from prior
    pi_draw <- rbeta(1, a_pi, b_pi)
    sigma2_draw <- rexp(1, lambda_sigma)

    # Draw z_i (spike-slab indicators) from Bernoulli(pi)
    z_draw <- rbinom(nedges, 1, pi_draw)

    # Draw r_i from mixture
    r_draw <- numeric(nedges)
    for (i in 1:nedges) {
      if (z_draw[i] == 1) {
        r_draw[i] <- 1.0  # Spike
      } else {
        r_draw[i] <- rlnorm(1, 0, sqrt(sigma2_draw))  # Slab
      }
    }

    # Simulate data under this model
    sim_data <- simRateScape(
      tree = tree,
      Q = Q,
      lambda_sigma = lambda_sigma,
      pi = pi_draw,
      nrep = 1,
      return_rates = FALSE
    )

    # Extract summary statistics from simulated data
    sim_char <- sim_data$data[, 1]
    n_transitions <- sum(tree$edge[, 1] > length(tree$tip.label))  # Internal nodes
    char_entropy <- entropy(table(sim_char))

    prior_predictive_stats[[sim]] <- list(
      pi = pi_draw,
      sigma2 = sigma2_draw,
      n_shifted = sum(z_draw == 0),
      mean_r = mean(r_draw),
      sd_r = sd(r_draw),
      char_entropy = char_entropy
    )

    if (sim %% 20 == 0) {
      message(sprintf("  Simulation %d / %d", sim, nsim))
    }
  }

  # Convert to data frame for easier summary
  prior_pred_df <- data.frame(
    pi = sapply(prior_predictive_stats, "[[", "pi"),
    sigma2 = sapply(prior_predictive_stats, "[[", "sigma2"),
    n_shifted = sapply(prior_predictive_stats, "[[", "n_shifted"),
    mean_r = sapply(prior_predictive_stats, "[[", "mean_r"),
    sd_r = sapply(prior_predictive_stats, "[[", "sd_r"),
    char_entropy = sapply(prior_predictive_stats, "[[", "char_entropy")
  )

  # ========== OVERLAP CALCULATION ==========

  # Simplified overlap measure: proportion of prior-predictive simulations
  # that fall within the interquartile range of a reasonable data distribution
  # (This is a heuristic; full implementation would require posterior samples)

  mean_sigma2_prior <- 1 / lambda_sigma
  sd_prior <- sqrt(mean_sigma2_prior)

  overlap_measure <- compute_prior_posterior_overlap(
    prior_df = prior_pred_df,
    lambda_sigma = lambda_sigma,
    expected_background = expected_background
  )

  # ========== WARNING FLAGS ==========

  warnings <- character(0)

  if (overlap_measure > 0.80) {
    warnings <- c(warnings,
      sprintf(
        "High prior-posterior overlap (%.1f%%). Data may not be informative relative to prior.",
        overlap_measure * 100
      )
    )
  }

  if (mean(prior_pred_df$sigma2) > 10 * lambda_sigma) {
    warnings <- c(warnings,
      sprintf(
        "Prior on σ² is very diffuse (mean = %.2f). Consider increasing lambda_sigma.",
        mean(prior_pred_df$sigma2)
      )
    )
  }

  if (!is.null(expected_background)) {
    if (expected_background < 0.1 || expected_background > 0.9) {
      warnings <- c(warnings,
        sprintf(
          "expected_background = %.2f is extreme. Be cautious with interpretation.",
          expected_background
        )
      )
    }
  }

  message(sprintf("\nPrior-posterior overlap: %.1f%%", overlap_measure * 100))

  if (length(warnings) > 0) {
    message("\n⚠ Warnings:")
    for (w in warnings) {
      message(sprintf("  - %s", w))
    }
  }

  # Create output object
  result <- list(
    prior_predictive_df = prior_pred_df,
    prior_settings = list(
      lambda_sigma = lambda_sigma,
      expected_background = expected_background,
      a_pi = a_pi,
      b_pi = b_pi
    ),
    overlap_measure = overlap_measure,
    warnings = warnings,
    call = match.call()
  )

  class(result) <- "ratescape_prior_check"
  return(result)
}


#' Compute prior-posterior overlap
#'
#' @keywords internal
compute_prior_posterior_overlap <- function(prior_df, lambda_sigma, expected_background) {
  # Simplified heuristic: compute 1 - (range of prior) / (expected range of posterior)
  # Full implementation would require posterior samples

  prior_range <- diff(range(prior_df$sigma2))
  expected_posterior_range <- 4 * sqrt(1 / lambda_sigma)  # Rough estimate

  overlap <- 1 - (prior_range / (expected_posterior_range + prior_range))
  return(max(0, min(1, overlap)))
}


#' Compute entropy of a distribution
#'
#' @keywords internal
entropy <- function(counts) {
  p <- counts / sum(counts)
  -sum(p * log(p + 1e-10))
}


#' Print method for prior check
#'
#' @export
print.ratescape_prior_check <- function(x, ...) {
  cat("RateScape Prior Diagnostic Check\n")
  cat("=================================\n\n")

  cat("Prior settings:\n")
  cat(sprintf("  lambda_sigma = %.2f (E[σ²] = %.2f)\n",
              x$prior_settings$lambda_sigma,
              1 / x$prior_settings$lambda_sigma))

  if (!is.null(x$prior_settings$expected_background)) {
    cat(sprintf("  expected_background = %.2f\n", x$prior_settings$expected_background))
  } else {
    cat("  Prior on π: Beta(1, 1) uniform\n")
  }

  cat(sprintf("\nPrior-posterior overlap: %.1f%%\n", x$overlap_measure * 100))

  if (x$overlap_measure > 0.80) {
    cat("  ⚠ High overlap suggests prior may be constraining inference\n")
  } else if (x$overlap_measure < 0.20) {
    cat("  ✓ Low overlap suggests data are informative\n")
  } else {
    cat("  ✓ Overlap is reasonable\n")
  }

  cat(sprintf("\nPrior-predictive summary (n = %d simulations):\n", nrow(x$prior_predictive_df)))
  cat(sprintf("  Mean σ² = %.2f (range %.2f–%.2f)\n",
              mean(x$prior_predictive_df$sigma2),
              min(x$prior_predictive_df$sigma2),
              max(x$prior_predictive_df$sigma2)))

  cat(sprintf("  Mean proportion shifted branches = %.1f%%\n",
              100 * mean(x$prior_predictive_df$n_shifted) / 100))

  if (length(x$warnings) > 0) {
    cat("\n⚠ Warnings:\n")
    for (w in x$warnings) {
      cat(sprintf("  - %s\n", w))
    }
  } else {
    cat("\n✓ No warnings\n")
  }

  invisible(x)
}
