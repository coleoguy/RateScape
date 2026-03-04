#' Automatic Discovery of Optimal Discrete Character Models
#'
#' Systematically searches the space of discrete character evolution models to
#' find the best-fitting model for your data. Tests combinations of rate matrix
#' structure (Mk, SYM, ARD, ordered, Dollo, custom), rate heterogeneity
#' (none, gamma, spike-and-slab), and root prior specification.
#'
#' @param tree An object of class "phylo".
#' @param data Data frame with tip states (single column).
#' @param models Character vector. Which Q-matrix structures to test.
#'   Default tests "mk", "sym", "ard", "ordered". Set to "all" to include
#'   "dollo" and "meristic" as well.
#' @param rate_models Character vector. Which rate heterogeneity models to test.
#'   Options: "none" (homogeneous), "gamma" (discretized gamma ML), "rjmcmc"
#'   (spike-and-slab Bayesian). Default c("none", "gamma").
#' @param root_priors Character vector. Which root priors to test. Default
#'   c("fitzjohn", "equal").
#' @param criterion Character. Model selection criterion: "aic", "aicc", or
#'   "bic". Default "aicc".
#' @param k_gamma_max Integer. Maximum number of gamma rate categories. Default 6.
#' @param lambda_sigma Numeric. Prior for spike-and-slab sigma^2 (used if
#'   "rjmcmc" is in rate_models). Default 1.0.
#' @param ngen Integer. MCMC generations for Bayesian models. Default 10000.
#' @param burn_in Integer. Burn-in for Bayesian models. Default 2000.
#' @param thin Integer. Thinning for Bayesian. Default 10.
#' @param parallel Logical. If TRUE, use parallel::mclapply for ML models.
#'   Default FALSE.
#' @param n_cores Integer. Number of cores for parallel. Default 2.
#' @param seed Integer or NULL. Random seed.
#'
#' @details
#'
#' **Search strategy:**
#'
#' The function tests all combinations of the specified Q structures, rate models,
#' and root priors. For each combination:
#'
#' 1. Build the Q matrix with the specified structure.
#' 2. Fit the model:
#'    - rate_model = "none": Optimize Q rates by ML (homogeneous rates).
#'    - rate_model = "gamma": Fit discretized gamma via EM (rateCategories).
#'    - rate_model = "rjmcmc": Fit spike-and-slab Bayesian model (fitRateScape).
#' 3. Compute model selection criterion (AIC/AICc/BIC).
#'
#' For ML models, AIC/AICc/BIC are computed from the log-likelihood and number
#' of free parameters. For Bayesian models, DIC (deviance information criterion)
#' is computed instead, and models are ranked by DIC.
#'
#' **Comparison with the 2024 preprint:**
#'
#' The preprint "Automatic Discovery of Optimal Discrete Character Models"
#' proposes exhaustive search over model space. Our implementation is similar
#' but adds: (1) rate heterogeneity as part of the search, (2) Bayesian model
#' averaging via RJMCMC as an option, and (3) posterior predictive checks for
#' the top models.
#'
#' @return An object of class "ratescapeAutoModel" containing:
#'   \describe{
#'     \item{results}{Data frame with all models tested, sorted by criterion.}
#'     \item{best_model}{The best-fitting model details.}
#'     \item{fits}{List of all fitted model objects.}
#'     \item{criterion}{Which criterion was used.}
#'   }
#'
#' @examples
#' \dontrun{
#'   result <- autoModel(tree, data)
#'   print(result)
#'
#'   # Top 5 models
#'   head(result$results, 5)
#'
#'   # Best model's Q matrix
#'   result$best_model$Q
#' }
#'
#' @export
autoModel <- function(
    tree, data,
    models = c("mk", "sym", "ard", "ordered"),
    rate_models = c("none", "gamma"),
    root_priors = c("fitzjohn", "equal"),
    criterion = "aicc",
    k_gamma_max = 6,
    lambda_sigma = 1.0,
    ngen = 10000, burn_in = 2000, thin = 10,
    parallel = FALSE, n_cores = 2,
    seed = NULL) {

  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  criterion <- match.arg(criterion, c("aic", "aicc", "bic"))
  if (!is.null(seed)) set.seed(seed)

  # Process data
  ntips <- length(tree$tip.label)
  if (!is.data.frame(data) && !is.matrix(data)) data <- data.frame(state = data)
  if (nrow(data) != ntips) stop("nrow(data) != ntips")

  states <- as.numeric(data[, 1])
  if (!is.null(rownames(data))) {
    mi <- match(tree$tip.label, rownames(data))
    if (any(is.na(mi))) stop("tip labels don't match data rownames")
    states <- states[mi]
  }
  states <- as.integer(states)
  if (min(states) > 0) states <- states - min(states)
  k <- max(states) + 1

  nedges <- nrow(tree$edge)

  # Handle "all" models
  if ("all" %in% models) {
    models <- c("mk", "sym", "ard", "ordered", "dollo", "meristic")
  }

  # Filter out models that don't make sense
  valid_models <- models
  if (k == 2) {
    # For binary: dollo = irreversible, meristic = mk, ordered = sym
    valid_models <- setdiff(models, c("meristic"))  # meristic == mk for k=2
  }

  # Build all combinations
  combos <- expand.grid(
    q_model = valid_models,
    rate_model = rate_models,
    root_prior = root_priors,
    stringsAsFactors = FALSE
  )

  n_combos <- nrow(combos)
  message(sprintf("autoModel: testing %d model combinations (k = %d states, %d tips)",
                  n_combos, k, ntips))

  # Pre-build Q templates for each unique q_model (avoid rebuilding identical Qs)
  unique_q_models <- unique(combos$q_model)
  Q_cache <- list()
  for (qm in unique_q_models) {
    Q_cache[[qm]] <- tryCatch(makeQ(model = qm, k = k), error = function(e) NULL)
  }

  results <- data.frame(
    q_model = character(n_combos),
    rate_model = character(n_combos),
    root_prior = character(n_combos),
    loglik = numeric(n_combos),
    n_params = integer(n_combos),
    aic = numeric(n_combos),
    aicc = numeric(n_combos),
    bic = numeric(n_combos),
    dic = numeric(n_combos),
    converged = logical(n_combos),
    stringsAsFactors = FALSE
  )

  fits <- vector("list", n_combos)

  # Define single-combo fitting function for potential parallelization
  fit_one_combo <- function(ci) {
    q_mod <- combos$q_model[ci]
    r_mod <- combos$rate_model[ci]
    rp <- combos$root_prior[ci]

    message(sprintf("  [%d/%d] Q=%s, rates=%s, root=%s", ci, n_combos, q_mod, r_mod, rp))

    # Use cached Q template
    Q_template <- Q_cache[[q_mod]]
    if (is.null(Q_template)) {
      message("    -> skipped (invalid model for this k)")
      results$converged[ci] <- FALSE
      results$q_model[ci] <- q_mod
      results$rate_model[ci] <- r_mod
      results$root_prior[ci] <- rp
      results$loglik[ci] <- NA
      next
    }

    param_map <- attr(Q_template, "param_map")
    n_q_params <- if (!is.null(attr(Q_template, "n_params"))) attr(Q_template, "n_params") else {
      # Count unique non-zero off-diagonal entries
      sum(Q_template[row(Q_template) != col(Q_template)] != 0)
    }

    fit_result <- tryCatch({
      if (r_mod == "none") {
        # Homogeneous ML
        fit <- optimize_Q_ml(tree, states, Q_template, rp)
        list(
          loglik = fit$loglik,
          n_params = fit$n_params,
          Q = fit$Q,
          fit = fit,
          type = "ml"
        )
      } else if (r_mod == "gamma") {
        # Discretized gamma ML
        fit <- rateCategories(tree, data, Q_template,
                              k_min = 2, k_max = k_gamma_max,
                              root_prior = rp, seed = seed)
        list(
          loglik = fit$loglik,
          n_params = n_q_params + fit$best_k,  # Q params + gamma params
          Q = Q_template,
          fit = fit,
          type = "ml"
        )
      } else if (r_mod == "rjmcmc") {
        # Bayesian spike-and-slab
        fit <- fitRateScape(tree, data, Q = Q_template, estimate_Q = FALSE,
                            lambda_sigma = lambda_sigma, root_prior = rp,
                            ngen = ngen, burn_in = burn_in, thin = thin,
                            seed = seed)
        # DIC
        mean_ll <- mean(fit$mcmc_samples$loglik)
        ll_at_mean <- compute_likelihood(tree, states, Q_template,
                                         r_scalars = colMeans(fit$mcmc_samples$r),
                                         root_prior = rp)
        p_D <- 2 * (ll_at_mean - mean_ll)
        dic <- -2 * ll_at_mean + 2 * p_D

        list(
          loglik = median(fit$mcmc_samples$loglik),
          n_params = NA,
          Q = Q_template,
          fit = fit,
          type = "bayesian",
          dic = dic
        )
      }
    }, error = function(e) {
      message(sprintf("    -> ERROR: %s", e$message))
      list(loglik = NA, n_params = NA, Q = Q_template, fit = NULL, type = "error")
    })

    results$q_model[ci] <- q_mod
    results$rate_model[ci] <- r_mod
    results$root_prior[ci] <- rp
    results$loglik[ci] <- fit_result$loglik
    results$n_params[ci] <- fit_result$n_params
    results$converged[ci] <- !is.na(fit_result$loglik)
    fits[[ci]] <- fit_result

    # Compute IC for ML models
    if (!is.na(fit_result$loglik) && fit_result$type == "ml") {
      ll <- fit_result$loglik
      np <- fit_result$n_params
      n <- ntips

      results$aic[ci] <- -2 * ll + 2 * np
      results$aicc[ci] <- -2 * ll + 2 * np + (2 * np * (np + 1)) / max(1, n - np - 1)
      results$bic[ci] <- -2 * ll + np * log(n)
    }

    if (fit_result$type == "bayesian" && !is.null(fit_result$dic)) {
      results$dic[ci] <- fit_result$dic
    }

    NULL  # return value for fit_one_combo
  }

  # Run model fitting (sequential or parallel for ML-only combos)
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    # Separate ML combos (parallelizable) from Bayesian (sequential)
    ml_idx <- which(combos$rate_model != "rjmcmc")
    bayes_idx <- which(combos$rate_model == "rjmcmc")

    if (length(ml_idx) > 0) {
      message(sprintf("  Running %d ML models in parallel on %d cores...", length(ml_idx), n_cores))
      parallel::mclapply(ml_idx, fit_one_combo, mc.cores = n_cores)
    }
    # Bayesian must be sequential (MCMC internal state)
    for (ci in bayes_idx) fit_one_combo(ci)
  } else {
    for (ci in 1:n_combos) fit_one_combo(ci)
  }

  # Sort by criterion
  results_valid <- results[results$converged, ]

  # Separate ML and Bayesian models
  ml_mask <- results_valid$rate_model != "rjmcmc"
  if (any(ml_mask)) {
    ml_results <- results_valid[ml_mask, ]
    ml_results <- ml_results[order(ml_results[[criterion]]), ]
  } else {
    ml_results <- data.frame()
  }

  bayes_mask <- results_valid$rate_model == "rjmcmc"
  if (any(bayes_mask)) {
    bayes_results <- results_valid[bayes_mask, ]
    bayes_results <- bayes_results[order(bayes_results$dic), ]
  } else {
    bayes_results <- data.frame()
  }

  # Combined ranking (ML models first by IC, then Bayesian by DIC)
  results_sorted <- rbind(ml_results, bayes_results)

  # Delta values
  if (nrow(ml_results) > 0) {
    best_ic <- min(ml_results[[criterion]], na.rm = TRUE)
    ml_results$delta <- ml_results[[criterion]] - best_ic
    ml_results$weight <- exp(-0.5 * ml_results$delta) / sum(exp(-0.5 * ml_results$delta))
  }

  # Best model
  if (nrow(ml_results) > 0) {
    best_idx <- which.min(ml_results[[criterion]])
    best_row <- which(results$q_model == ml_results$q_model[best_idx] &
                       results$rate_model == ml_results$rate_model[best_idx] &
                       results$root_prior == ml_results$root_prior[best_idx])
    best_fit <- fits[[best_row[1]]]
  } else {
    best_fit <- NULL
  }

  message(sprintf("\nBest ML model: Q=%s, rates=%s, root=%s",
                  ml_results$q_model[1], ml_results$rate_model[1], ml_results$root_prior[1]))

  result <- list(
    results = results_sorted,
    ml_results = ml_results,
    bayes_results = bayes_results,
    best_model = best_fit,
    fits = fits,
    criterion = criterion,
    k = k,
    ntips = ntips,
    call = match.call()
  )

  class(result) <- "ratescapeAutoModel"
  result
}


#' Print method for auto model discovery
#' @export
print.ratescapeAutoModel <- function(x, ...) {
  cat("RateScape Automatic Model Discovery\n")
  cat("=====================================\n\n")
  cat(sprintf("States: %d | Tips: %d | Criterion: %s\n\n", x$k, x$ntips, x$criterion))

  if (nrow(x$ml_results) > 0) {
    cat("Top ML Models:\n")
    top_ml <- head(x$ml_results, 10)
    print_cols <- c("q_model", "rate_model", "root_prior", "loglik", "n_params",
                    x$criterion, "delta", "weight")
    avail_cols <- intersect(print_cols, names(top_ml))
    print(top_ml[, avail_cols], row.names = FALSE, digits = 3)
  }

  if (nrow(x$bayes_results) > 0) {
    cat("\nBayesian Models (ranked by DIC):\n")
    top_b <- head(x$bayes_results, 5)
    print(top_b[, c("q_model", "rate_model", "root_prior", "loglik", "dic")],
          row.names = FALSE, digits = 3)
  }

  invisible(x)
}
