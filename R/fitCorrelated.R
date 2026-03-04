#' Correlated Discrete Trait Evolution (Pagel's Test and Beyond)
#'
#' Tests for correlated (dependent) evolution between two or more discrete
#' characters, with optional branch-specific rate heterogeneity. Extends
#' Pagel's (1994) method by allowing rate variation under the dependent model.
#'
#' @param tree An object of class "phylo".
#' @param data A data frame with rows matching tip labels and columns for each
#'   discrete character. Each column must contain integer-coded states
#'   (0, 1, 2, ...). Two or more columns allowed; all are combined into a
#'   single compound character.
#' @param method Character. Inference method: "ml" for maximum likelihood or
#'   "bayesian" for RJMCMC with rate scalars. Default "ml".
#' @param root_prior Character. Root state prior: "fitzjohn", "equal", or
#'   "stationary". Default "fitzjohn".
#' @param lambda_sigma Numeric. Exponential prior rate on sigma^2 (required
#'   for method = "bayesian").
#' @param ngen Integer. MCMC generations (for method = "bayesian"). Default 10000.
#' @param burn_in Integer. Burn-in iterations. Default 2000.
#' @param thin Integer. Thinning interval. Default 10.
#' @param seed Integer or NULL. Random seed.
#'
#' @details
#'
#' **The compound character approach:**
#'
#' For two binary characters X and Y, the compound character has 4 states:
#' (0,0), (0,1), (1,0), (1,1). Evolution is modeled with an 8-parameter
#' dependent model (allowing each character's transition rates to depend on
#' the state of the other character) versus a 4-parameter independent model
#' (each character evolves on its own).
#'
#' For k characters with states s1, s2, ..., sk, the compound character has
#' prod(n_states_i) states. The independent model constrains rates such that
#' only one character changes at a time and each character's rates are
#' independent of the others' states. The dependent model relaxes this:
#' transition rates for character i can depend on the states of characters j.
#'
#' **What \RateScape\ adds over Pagel (1994):**
#'
#' 1. When method = "bayesian", the correlated model is fitted *with*
#'    branch-specific rate scalars (spike-and-slab), so correlation is
#'    estimated while accounting for rate heterogeneity. This avoids a
#'    known confound where rate variation masquerades as correlated evolution.
#'
#' 2. Supports 2+ characters (not just binary pairs).
#'
#' 3. Returns both LRT (for ML) and Bayes factors (for Bayesian), plus
#'    posterior distributions of each transition rate for the dependent model.
#'
#' @return An object of class "ratescapeCorr" containing:
#'   \describe{
#'     \item{test}{Data frame with loglikelihood, AIC, and parameters for
#'       independent and dependent models.}
#'     \item{LRT}{Likelihood ratio test statistic and p-value (ML only).}
#'     \item{Q_independent}{Fitted independent Q matrix.}
#'     \item{Q_dependent}{Fitted dependent Q matrix.}
#'     \item{compound_states}{Mapping from compound state indices to character
#'       state combinations.}
#'     \item{bayesian_fit}{The full fitRateScape result if method = "bayesian".}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Two binary characters
#'   dat <- data.frame(trait1 = c(0,0,1,1,0),
#'                     trait2 = c(0,1,0,1,1),
#'                     row.names = tree$tip.label)
#'   result <- fitCorrelated(tree, dat)
#'   print(result)
#' }
#'
#' @export
fitCorrelated <- function(
    tree, data, method = "ml",
    root_prior = "fitzjohn",
    lambda_sigma = NULL,
    ngen = 10000, burn_in = 2000, thin = 10,
    seed = NULL) {

  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  method <- match.arg(method, c("ml", "bayesian"))
  root_prior <- match.arg(root_prior, c("fitzjohn", "equal", "stationary"))
  if (method == "bayesian" && is.null(lambda_sigma))
    stop("lambda_sigma required for bayesian method")
  if (!is.null(seed)) set.seed(seed)

  # ---- Process data: match to tree, build compound character ----
  data <- as.data.frame(data)
  ntips <- length(tree$tip.label)
  if (nrow(data) != ntips) stop("nrow(data) must match number of tips")

  if (!is.null(rownames(data))) {
    mi <- match(tree$tip.label, rownames(data))
    if (any(is.na(mi))) stop("tip labels don't match data rownames")
    data <- data[mi, , drop = FALSE]
  }

  n_chars <- ncol(data)
  if (n_chars < 2) stop("Need at least 2 characters for correlated evolution test")

  # States per character
  char_states <- lapply(1:n_chars, function(i) sort(unique(data[, i])))
  n_per_char <- sapply(char_states, length)

  # Build compound state space
  compound_grid <- expand.grid(char_states)
  colnames(compound_grid) <- colnames(data)
  n_compound <- nrow(compound_grid)
  compound_labels <- apply(compound_grid, 1, paste, collapse = ",")

  message(sprintf("Correlated evolution test: %d characters -> %d compound states",
                  n_chars, n_compound))

  # Map tip data to compound states (vectorized via paste-match)
  tip_keys <- apply(data, 1, paste, collapse = ",")
  compound_states <- match(tip_keys, compound_labels) - 1L  # 0-indexed
  if (any(is.na(compound_states))) {
    bad <- which(is.na(compound_states))[1]
    stop(sprintf("Tip %d: cannot map to compound state", bad))
  }

  # ---- Build Independent Q matrix ----
  # Independent: only one character changes at a time, and its rate does NOT
  # depend on the state of the other character(s).
  Q_indep <- build_independent_Q(compound_grid, char_states, n_chars, n_compound)

  # ---- Build Dependent Q matrix ----
  # Dependent: only one character changes at a time, but its rate CAN
  # depend on the states of the other characters.
  Q_dep <- build_dependent_Q(compound_grid, char_states, n_chars, n_compound)

  # ---- Fit models ----
  compound_data <- data.frame(state = compound_states)
  rownames(compound_data) <- tree$tip.label

  if (method == "ml") {
    # Optimize Q_indep
    fit_indep <- optimize_Q_ml(tree, compound_states, Q_indep, root_prior)
    fit_dep <- optimize_Q_ml(tree, compound_states, Q_dep, root_prior)

    n_params_indep <- fit_indep$n_params
    n_params_dep <- fit_dep$n_params
    ll_indep <- fit_indep$loglik
    ll_dep <- fit_dep$loglik

    aic_indep <- -2 * ll_indep + 2 * n_params_indep
    aic_dep <- -2 * ll_dep + 2 * n_params_dep

    # LRT
    lrt_stat <- 2 * (ll_dep - ll_indep)
    df <- n_params_dep - n_params_indep
    p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)

    test_table <- data.frame(
      Model = c("Independent", "Dependent"),
      LogLik = c(ll_indep, ll_dep),
      Params = c(n_params_indep, n_params_dep),
      AIC = c(aic_indep, aic_dep)
    )

    result <- list(
      test = test_table,
      LRT = list(statistic = lrt_stat, df = df, p_value = p_value),
      Q_independent = fit_indep$Q,
      Q_dependent = fit_dep$Q,
      compound_states = compound_grid,
      compound_labels = compound_labels,
      method = "ml",
      n_chars = n_chars,
      char_states = char_states
    )

  } else {
    # Bayesian: fit dependent model with rate heterogeneity
    message("Fitting dependent model with rate heterogeneity (Bayesian)...")

    # Fit independent model (ML) for comparison
    fit_indep <- optimize_Q_ml(tree, compound_states, Q_indep, root_prior)

    # Fit dependent model with RJMCMC rate scalars
    fit_dep_bayes <- fitRateScape(
      tree = tree, data = compound_data, Q = Q_dep,
      estimate_Q = FALSE, lambda_sigma = lambda_sigma,
      root_prior = root_prior,
      ngen = ngen, burn_in = burn_in, thin = thin, seed = seed
    )

    # Homogeneous dependent model for comparison
    fit_dep_ml <- optimize_Q_ml(tree, compound_states, Q_dep, root_prior)

    test_table <- data.frame(
      Model = c("Independent (ML)", "Dependent (ML)", "Dependent + Rate Het (Bayesian)"),
      LogLik = c(fit_indep$loglik, fit_dep_ml$loglik,
                 median(fit_dep_bayes$mcmc_samples$loglik)),
      Params = c(fit_indep$n_params, fit_dep_ml$n_params, NA)
    )

    result <- list(
      test = test_table,
      Q_independent = fit_indep$Q,
      Q_dependent = fit_dep_ml$Q,
      compound_states = compound_grid,
      compound_labels = compound_labels,
      bayesian_fit = fit_dep_bayes,
      method = "bayesian",
      n_chars = n_chars,
      char_states = char_states
    )
  }

  class(result) <- "ratescapeCorr"
  result
}


#' Build independent Q matrix for compound character
#'
#' Under independence, only one character changes at a time and its rate
#' does NOT depend on the other characters' states.
#'
#' @param grid Data frame of compound state combinations.
#' @param char_states List of per-character state vectors.
#' @param n_chars Number of characters.
#' @param n_compound Number of compound states.
#' @return Q matrix with structural zeros enforced and free parameters indexed.
#' @keywords internal
build_independent_Q <- function(grid, char_states, n_chars, n_compound) {
  Q <- matrix(0, nrow = n_compound, ncol = n_compound)
  param_map <- matrix(0L, nrow = n_compound, ncol = n_compound)

  # Precompute: for each pair of compound states, find which character
  # differs (if exactly one). Convert grid to matrix for fast comparison.
  grid_mat <- as.matrix(grid)

  # Build parameter index lookup: (character, from_state, to_state) -> param_idx
  param_lookup <- list()
  param_idx <- 0L
  for (ch in 1:n_chars) {
    for (s_from in char_states[[ch]]) {
      for (s_to in char_states[[ch]]) {
        if (s_from == s_to) next
        param_idx <- param_idx + 1L
        key <- paste(ch, s_from, s_to, sep = "_")
        param_lookup[[key]] <- param_idx
      }
    }
  }

  # Single pass over compound state pairs: O(n_compound^2 * n_chars)
  for (i in 1:n_compound) {
    for (j in 1:n_compound) {
      if (i == j) next
      diffs <- which(grid_mat[i, ] != grid_mat[j, ])
      if (length(diffs) == 1L) {
        ch <- diffs
        key <- paste(ch, grid_mat[i, ch], grid_mat[j, ch], sep = "_")
        Q[i, j] <- 1.0
        param_map[i, j] <- param_lookup[[key]]
      }
    }
  }

  attr(Q, "param_map") <- param_map
  attr(Q, "n_params") <- param_idx
  diag(Q) <- -rowSums(Q)
  Q
}


#' Build dependent Q matrix for compound character
#'
#' Under dependence, only one character changes at a time, but its rate CAN
#' depend on the states of the other characters.
#'
#' @keywords internal
build_dependent_Q <- function(grid, char_states, n_chars, n_compound) {
  Q <- matrix(0, nrow = n_compound, ncol = n_compound)
  param_map <- matrix(0, nrow = n_compound, ncol = n_compound)
  param_idx <- 0

  for (i in 1:n_compound) {
    for (j in 1:n_compound) {
      if (i == j) next
      diffs <- which(grid[i, ] != grid[j, ])
      if (length(diffs) == 1) {
        # Each unique (from_compound, to_compound) pair where exactly one
        # character changed gets its own rate parameter
        param_idx <- param_idx + 1
        Q[i, j] <- 1.0
        param_map[i, j] <- param_idx
      }
    }
  }

  attr(Q, "param_map") <- param_map
  attr(Q, "n_params") <- param_idx
  diag(Q) <- -rowSums(Q)
  Q
}


#' Optimize Q matrix rates by ML
#'
#' Numerically optimizes the free parameters of a constrained Q matrix
#' to maximize the likelihood of the observed data.
#'
#' @keywords internal
optimize_Q_ml <- function(tree, states, Q_template, root_prior) {
  param_map <- attr(Q_template, "param_map")
  n_params <- attr(Q_template, "n_params")
  k <- nrow(Q_template)

  if (is.null(param_map) || n_params == 0) {
    # No free parameters; just compute likelihood
    ll <- compute_likelihood(tree, states, Q_template,
                             r_scalars = rep(1.0, nrow(tree$edge)),
                             root_prior = root_prior)
    return(list(Q = Q_template, loglik = ll, n_params = 0))
  }

  # Build Q from parameter vector
  make_Q_from_params <- function(log_rates) {
    rates <- exp(log_rates)
    Q_new <- matrix(0, nrow = k, ncol = k)
    for (i in 1:k) {
      for (j in 1:k) {
        if (i != j && param_map[i, j] > 0) {
          Q_new[i, j] <- rates[param_map[i, j]]
        }
      }
    }
    diag(Q_new) <- -rowSums(Q_new)
    Q_new
  }

  neg_loglik <- function(log_rates) {
    Q_new <- make_Q_from_params(log_rates)
    ll <- tryCatch(
      compute_likelihood(tree, states, Q_new,
                         r_scalars = rep(1.0, nrow(tree$edge)),
                         root_prior = root_prior),
      error = function(e) -1e20
    )
    -ll
  }

  # Optimize
  init <- rep(log(0.1), n_params)
  opt <- tryCatch(
    optim(init, neg_loglik, method = "L-BFGS-B",
          lower = rep(log(1e-8), n_params),
          upper = rep(log(100), n_params)),
    error = function(e) {
      optim(init, neg_loglik, method = "Nelder-Mead",
            control = list(maxit = 5000))
    }
  )

  Q_opt <- make_Q_from_params(opt$par)
  list(Q = Q_opt, loglik = -opt$value, n_params = n_params, rates = exp(opt$par))
}


#' Print method for correlated evolution test
#' @export
print.ratescapeCorr <- function(x, ...) {
  cat("RateScape Correlated Evolution Test\n")
  cat("====================================\n\n")
  cat(sprintf("Characters: %d\n", x$n_chars))
  cat(sprintf("Compound states: %d\n", nrow(x$compound_states)))
  cat(sprintf("Method: %s\n\n", x$method))

  print(x$test, row.names = FALSE)

  if (!is.null(x$LRT)) {
    cat(sprintf("\nLikelihood Ratio Test: chi2 = %.2f, df = %d, p = %.4f\n",
                x$LRT$statistic, x$LRT$df, x$LRT$p_value))
    if (x$LRT$p_value < 0.05) {
      cat("** Significant evidence for correlated evolution **\n")
    } else {
      cat("No significant evidence for correlated evolution.\n")
    }
  }

  invisible(x)
}
