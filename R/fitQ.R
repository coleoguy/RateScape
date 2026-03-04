#' Fit a Discrete Character Markov Model
#'
#' Fits a continuous-time Markov model (Mk, SYM, ARD, ordered, Dollo, meristic,
#' or custom) to discrete character data on a phylogeny using maximum likelihood
#' or Bayesian MCMC. This is the basic workhorse function for discrete trait
#' analysis—no rate heterogeneity, just a single Q matrix fit to the data.
#'
#' @param tree An object of class "phylo".
#' @param data Named vector of character states at tips (names match tip labels).
#'   Can be character, factor, or integer (0-indexed).
#' @param model Character: one of "mk", "sym", "ard", "ordered", "dollo",
#'   "irreversible", "meristic", or a Q matrix from \code{makeQ()}.
#' @param method Character: "ml" for maximum likelihood or "bayesian" for MCMC.
#' @param root_prior Character: "fitzjohn" (empirical Bayes), "equal" (flat),
#'   or "stationary" (equilibrium of Q). Default "fitzjohn".
#' @param ngen Integer. Number of MCMC generations (Bayesian only). Default 10000.
#' @param burn_in Numeric. Proportion of MCMC to discard. Default 0.25.
#' @param thin Integer. Thinning interval. Default 10.
#' @param rate_prior Character: "exponential" or "lognormal". Prior on Q rates
#'   for Bayesian. Default "exponential".
#' @param rate_prior_mean Numeric. Mean of exponential prior (or meanlog for
#'   lognormal). Default 1.0.
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param ancestral Logical. If TRUE, compute marginal ancestral state
#'   probabilities at internal nodes. Default TRUE.
#'
#' @return An object of class "ratescapeQ" with components:
#'   \item{Q}{Fitted transition rate matrix.}
#'   \item{loglik}{Log-likelihood at the MLE (ML) or posterior median (Bayesian).}
#'   \item{n_params}{Number of free parameters in Q.}
#'   \item{AIC}{AIC (ML only).}
#'   \item{AICc}{Corrected AIC (ML only).}
#'   \item{BIC}{BIC (ML only).}
#'   \item{rates}{Named vector of estimated transition rates.}
#'   \item{ancestral_states}{Matrix of marginal ancestral state probabilities
#'     at internal nodes (if \code{ancestral = TRUE}).}
#'   \item{mcmc_samples}{Data frame of posterior samples (Bayesian only).}
#'   \item{method}{Character: "ml" or "bayesian".}
#'   \item{model}{Model name or "custom".}
#'   \item{states}{Character state labels.}
#'
#' @examples
#' \dontrun{
#' # Binary character, equal-rates model
#' fit_mk <- fitQ(tree, data, model = "mk")
#' print(fit_mk)
#'
#' # All-rates-different, ML
#' fit_ard <- fitQ(tree, data, model = "ard", method = "ml")
#'
#' # Ordered model for 4-state character, Bayesian
#' fit_ord <- fitQ(tree, data, model = "ordered", method = "bayesian",
#'                 ngen = 50000)
#'
#' # Custom Q from makeQ
#' Q_custom <- makeQ(4, model = "ordered", constraints = my_constraints)
#' fit_custom <- fitQ(tree, data, model = Q_custom)
#' }
#'
#' @export
fitQ <- function(tree, data, model = "sym",
                 method = c("ml", "bayesian"),
                 root_prior = c("fitzjohn", "equal", "stationary"),
                 ngen = 10000, burn_in = 0.25, thin = 10,
                 rate_prior = c("exponential", "lognormal"),
                 rate_prior_mean = 1.0,
                 seed = NULL, ancestral = TRUE) {

  method <- match.arg(method)
  root_prior <- match.arg(root_prior)
  rate_prior <- match.arg(rate_prior)
  if (!is.null(seed)) set.seed(seed)

  # ---- Parse data ----
  if (!inherits(tree, "phylo")) stop("tree must be of class 'phylo'")
  parsed <- parse_discrete_data(tree, data)
  states_int <- parsed$states_int
  state_labels <- parsed$state_labels
  k <- parsed$k

  # ---- Build Q template ----
  if (is.matrix(model)) {
    # User passed a Q matrix directly (e.g. from makeQ)
    Q_template <- model
    model_name <- ifelse(!is.null(attr(model, "model")), attr(model, "model"), "custom")
  } else {
    model_name <- tolower(model)
    Q_template <- makeQ(k = k, model = model_name)
  }

  if (nrow(Q_template) != k) {
    stop("Q matrix dimension (", nrow(Q_template), ") does not match ",
         "number of observed states (", k, ")")
  }

  param_map <- attr(Q_template, "param_map")
  n_params <- attr(Q_template, "n_params")
  if (is.null(n_params)) n_params <- 0

  # ---- ML inference ----
  if (method == "ml") {
    fit <- optimize_Q_full(tree, states_int, Q_template, root_prior)

    # Information criteria
    n <- length(tree$tip.label)
    loglik <- fit$loglik
    aic <- -2 * loglik + 2 * n_params
    aicc <- aic + (2 * n_params * (n_params + 1)) / max(n - n_params - 1, 1)
    bic <- -2 * loglik + n_params * log(n)

    # Extract named rates
    rates_vec <- extract_rates(fit$Q, state_labels)

    # Marginal ancestral states
    anc <- NULL
    if (ancestral) {
      anc <- marginal_ancestral_states(tree, states_int, fit$Q, root_prior)
      rownames(anc) <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
      colnames(anc) <- state_labels
    }

    result <- list(
      Q = fit$Q,
      loglik = loglik,
      n_params = n_params,
      AIC = aic,
      AICc = aicc,
      BIC = bic,
      rates = rates_vec,
      ancestral_states = anc,
      method = "ml",
      model = model_name,
      states = state_labels,
      tree = tree
    )

  # ---- Bayesian inference ----
  } else {
    result <- bayesian_Q_mcmc(tree, states_int, Q_template, k, state_labels,
                              root_prior, ngen, burn_in, thin,
                              rate_prior, rate_prior_mean, ancestral,
                              model_name)
  }

  class(result) <- "ratescapeQ"
  result
}


#' Parse discrete character data into integer states
#' @keywords internal
parse_discrete_data <- function(tree, data) {
  if (!is.null(names(data))) {
    # Match data to tip order
    if (!all(tree$tip.label %in% names(data))) {
      stop("Not all tip labels found in data names")
    }
    data <- data[tree$tip.label]
  }
  if (length(data) != length(tree$tip.label)) {
    stop("data must have one entry per tip")
  }

  # Convert to integer states
  if (is.factor(data)) {
    state_labels <- levels(data)
    states_int <- as.integer(data) - 1L
  } else if (is.character(data)) {
    state_labels <- sort(unique(data))
    states_int <- match(data, state_labels) - 1L
  } else {
    states_int <- as.integer(data)
    if (min(states_int) > 0) states_int <- states_int - min(states_int)
    state_labels <- as.character(sort(unique(states_int)))
  }

  k <- length(unique(states_int))
  list(states_int = states_int, state_labels = state_labels, k = k)
}


#' Optimize Q matrix by ML (general version for fitQ)
#' @keywords internal
optimize_Q_full <- function(tree, states, Q_template, root_prior) {
  param_map <- attr(Q_template, "param_map")
  n_params <- attr(Q_template, "n_params")
  k <- nrow(Q_template)

  if (is.null(param_map) || n_params == 0) {
    ll <- compute_likelihood(tree, states, Q_template,
                             r_scalars = rep(1.0, nrow(tree$edge)),
                             root_prior = root_prior)
    return(list(Q = Q_template, loglik = ll))
  }

  # Build Q from log-rate parameter vector (VECTORIZED)
  make_Q <- function(log_rates) {
    rates <- exp(log_rates)
    Q_new <- matrix(0, nrow = k, ncol = k)
    # Vectorized: find all cells with param_map > 0 and off-diagonal
    idx <- which(param_map > 0 & row(param_map) != col(param_map))
    Q_new[idx] <- rates[param_map[idx]]
    diag(Q_new) <- -rowSums(Q_new)
    Q_new
  }

  neg_loglik <- function(log_rates) {
    Q_new <- make_Q(log_rates)
    ll <- tryCatch(
      compute_likelihood(tree, states, Q_new,
                         r_scalars = rep(1.0, nrow(tree$edge)),
                         root_prior = root_prior),
      error = function(e) -1e20
    )
    -ll
  }

  # Multiple starts to avoid local optima
  best_val <- Inf
  best_opt <- NULL

  starts <- list(rep(log(0.1), n_params))
  if (n_params <= 20) {
    starts <- c(starts, list(
      rep(log(1.0), n_params),
      rep(log(0.01), n_params),
      runif(n_params, log(0.01), log(1.0))
    ))
  }

  for (init in starts) {
    opt <- tryCatch(
      optim(init, neg_loglik, method = "L-BFGS-B",
            lower = rep(log(1e-10), n_params),
            upper = rep(log(200), n_params),
            control = list(maxit = 1000)),
      error = function(e) NULL
    )
    if (!is.null(opt) && opt$value < best_val) {
      best_val <- opt$value
      best_opt <- opt
    }
  }

  if (is.null(best_opt)) {
    warning("Optimization failed; returning template Q")
    ll <- compute_likelihood(tree, states, Q_template,
                             r_scalars = rep(1.0, nrow(tree$edge)),
                             root_prior = root_prior)
    return(list(Q = Q_template, loglik = ll))
  }

  Q_fit <- make_Q(best_opt$par)
  list(Q = Q_fit, loglik = -best_opt$value)
}


#' Extract named rate vector from fitted Q (VECTORIZED)
#' @keywords internal
extract_rates <- function(Q, state_labels) {
  k <- nrow(Q)
  mask <- (row(Q) != col(Q)) & (Q > 0)
  idx <- which(mask, arr.ind = TRUE)
  rates <- Q[mask]
  names(rates) <- paste0(state_labels[idx[, 1]], "->", state_labels[idx[, 2]])
  rates
}


#' Compute marginal ancestral state probabilities
#' @keywords internal
marginal_ancestral_states <- function(tree, states, Q, root_prior) {
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  k <- nrow(Q)
  r_scalars <- rep(1.0, nrow(tree$edge))

  # Full postorder pass
  info <- build_tree_info(tree, Q)
  L <- postorder_pass(info, states, Q, r_scalars)

  # Precompute edge lookup for O(1) edge matching instead of O(n) per edge
  edge_lookup <- paste(tree$edge[, 1], tree$edge[, 2], sep = "_")

  # Preorder pass for marginal probabilities
  # P(state at node) = P(data below | state) * P(state | data above)
  marginal <- matrix(0, nrow = tree$Nnode, ncol = k)

  # Root marginal
  root_L <- L[info$root, ]
  root_w <- get_root_weights(root_L, root_prior, Q, k)
  root_marginal <- root_L * root_w
  root_marginal <- root_marginal / sum(root_marginal)
  marginal[1, ] <- root_marginal

  # Preorder traversal
  tree_pre <- ape::reorder.phylo(tree, order = "cladewise")
  preorder_probs <- matrix(0, nrow = n_nodes, ncol = k)
  preorder_probs[info$root, ] <- root_marginal

  # Cache exponentiation info for faster matrix exponentials
  eig_Q <- info$eig

  for (i in 1:nrow(tree_pre$edge)) {
    parent <- tree_pre$edge[i, 1]
    child <- tree_pre$edge[i, 2]

    # Find original edge index using precomputed lookup (O(1) via match)
    orig_e <- match(paste(parent, child, sep = "_"), edge_lookup)
    t_branch <- tree$edge.length[orig_e]

    P <- fast_expm(eig_Q, t_branch)

    # Get parent's marginal as prior, then multiply by child's conditional
    parent_marginal <- preorder_probs[parent, ]

    if (child <= n_tips) {
      # Tip node - observed state
      child_joint <- (t(P) %*% parent_marginal) * L[child, ]
      child_joint <- as.vector(child_joint)
      if (sum(child_joint) > 0) child_joint <- child_joint / sum(child_joint)
    } else {
      # Internal node
      # Marginal = P(data below child | state) * P(state | parent marginal, branch)
      child_prior <- as.vector(t(P) %*% parent_marginal)
      child_joint <- child_prior * L[child, ]
      if (sum(child_joint) > 0) child_joint <- child_joint / sum(child_joint)
      node_idx <- child - n_tips
      marginal[node_idx, ] <- child_joint
    }
    preorder_probs[child, ] <- child_joint
  }

  marginal
}


#' Bayesian MCMC for Q matrix parameters
#' @keywords internal
bayesian_Q_mcmc <- function(tree, states, Q_template, k, state_labels,
                            root_prior, ngen, burn_in, thin,
                            rate_prior, rate_prior_mean, ancestral,
                            model_name) {

  param_map <- attr(Q_template, "param_map")
  n_params <- attr(Q_template, "n_params")
  if (is.null(n_params) || n_params == 0) {
    stop("Bayesian method requires a model with free parameters")
  }

  n_edges <- nrow(tree$edge)
  r_ones <- rep(1.0, n_edges)

  # Build Q from log-rate vector (VECTORIZED)
  make_Q <- function(log_rates) {
    rates <- exp(log_rates)
    Q_new <- matrix(0, nrow = k, ncol = k)
    idx <- which(param_map > 0 & row(param_map) != col(param_map))
    Q_new[idx] <- rates[param_map[idx]]
    diag(Q_new) <- -rowSums(Q_new)
    Q_new
  }

  # Log-prior on rates
  log_prior <- function(log_rates) {
    rates <- exp(log_rates)
    if (rate_prior == "exponential") {
      sum(dexp(rates, rate = 1 / rate_prior_mean, log = TRUE) + log_rates)
    } else {
      sum(dlnorm(rates, meanlog = log(rate_prior_mean), sdlog = 1.0, log = TRUE))
    }
  }

  # Initialize at ML estimate
  ml_fit <- optimize_Q_full(tree, states, Q_template, root_prior)
  current_log_rates <- log(pmax(extract_param_values(ml_fit$Q, param_map, n_params), 1e-8))
  current_Q <- ml_fit$Q
  current_ll <- ml_fit$loglik
  current_lp <- log_prior(current_log_rates)

  # Pre-build tree info ONCE to avoid repeated O(n) passes in MCMC
  tree_info <- build_tree_info(tree, current_Q)

  # MCMC storage
  n_samples <- floor((ngen * (1 - burn_in)) / thin)
  samples_loglik <- numeric(n_samples)
  samples_rates <- matrix(0, nrow = n_samples, ncol = n_params)

  # Adaptive proposal widths
  prop_sd <- rep(0.1, n_params)
  accept_count <- rep(0, n_params)
  attempt_count <- rep(0, n_params)

  burn_n <- floor(ngen * burn_in)
  sample_idx <- 0

  message("Running MCMC for Q-matrix estimation (", ngen, " generations)...")

  for (gen in 1:ngen) {
    # Update each parameter one at a time (component-wise MH)
    for (p in 1:n_params) {
      attempt_count[p] <- attempt_count[p] + 1

      proposed_log_rates <- current_log_rates
      proposed_log_rates[p] <- current_log_rates[p] + rnorm(1, 0, prop_sd[p])

      # Bounds check
      if (proposed_log_rates[p] < log(1e-10) || proposed_log_rates[p] > log(200)) next

      Q_prop <- make_Q(proposed_log_rates)

      # Use pre-built tree info for likelihood computation
      # Rebuild eigendecomp for new Q (necessary for new rates)
      tree_info_prop <- build_tree_info(tree, Q_prop)
      ll_prop <- tryCatch(
        {
          L <- postorder_pass(tree_info_prop, states, Q_prop, r_ones)
          root_L <- L[tree_info_prop$root, ]
          root_w <- get_root_weights(root_L, root_prior, Q_prop, k)
          sum(log(pmax(root_L * root_w, 1e-300)))
        },
        error = function(e) -1e20
      )
      lp_prop <- log_prior(proposed_log_rates)

      log_alpha <- (ll_prop + lp_prop) - (current_ll + current_lp)

      if (log(runif(1)) < log_alpha) {
        current_log_rates <- proposed_log_rates
        current_Q <- Q_prop
        current_ll <- ll_prop
        current_lp <- lp_prop
        tree_info <- tree_info_prop
        accept_count[p] <- accept_count[p] + 1
      }
    }

    # Adaptive tuning every 200 generations during burn-in
    if (gen <= burn_n && gen %% 200 == 0) {
      for (p in 1:n_params) {
        if (attempt_count[p] > 0) {
          acc_rate <- accept_count[p] / attempt_count[p]
          if (acc_rate < 0.2) prop_sd[p] <- prop_sd[p] * 0.8
          if (acc_rate > 0.5) prop_sd[p] <- prop_sd[p] * 1.2
        }
      }
      accept_count[] <- 0
      attempt_count[] <- 0
    }

    # Store samples
    if (gen > burn_n && (gen - burn_n) %% thin == 0) {
      sample_idx <- sample_idx + 1
      if (sample_idx <= n_samples) {
        samples_loglik[sample_idx] <- current_ll
        samples_rates[sample_idx, ] <- exp(current_log_rates)
      }
    }

    if (gen %% 2000 == 0) {
      message("  Generation ", gen, "/", ngen, " | logLik = ",
              round(current_ll, 2))
    }
  }

  # Trim if needed
  if (sample_idx < n_samples) {
    samples_loglik <- samples_loglik[1:sample_idx]
    samples_rates <- samples_rates[1:sample_idx, , drop = FALSE]
  }

  # Summary Q from posterior medians
  median_rates <- apply(samples_rates, 2, median)
  Q_median <- make_Q(log(median_rates))

  ll_median <- compute_likelihood(tree, states, Q_median, r_ones, root_prior)

  # Rate names (VECTORIZED)
  rate_names <- make_rate_names(param_map, state_labels, k)
  colnames(samples_rates) <- rate_names

  # Marginal ancestral states at posterior median Q
  anc <- NULL
  if (ancestral) {
    anc <- marginal_ancestral_states(tree, states, Q_median, root_prior)
    rownames(anc) <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
    colnames(anc) <- state_labels
  }

  # Posterior summary
  rate_summary <- data.frame(
    parameter = rate_names,
    mean = colMeans(samples_rates),
    median = apply(samples_rates, 2, median),
    lower_95 = apply(samples_rates, 2, quantile, 0.025),
    upper_95 = apply(samples_rates, 2, quantile, 0.975)
  )

  mcmc_df <- as.data.frame(samples_rates)
  mcmc_df$loglik <- samples_loglik

  list(
    Q = Q_median,
    loglik = ll_median,
    n_params = n_params,
    rates = setNames(median_rates, rate_names),
    rate_summary = rate_summary,
    ancestral_states = anc,
    mcmc_samples = mcmc_df,
    method = "bayesian",
    model = model_name,
    states = state_labels,
    tree = tree
  )
}


#' Extract parameter values from Q using param_map (VECTORIZED)
#' @keywords internal
extract_param_values <- function(Q, param_map, n_params) {
  k <- nrow(Q)
  vals <- numeric(n_params)
  # Vectorized: find all off-diagonal cells with param_map > 0
  idx <- which(param_map > 0 & row(param_map) != col(param_map))
  vals[param_map[idx]] <- Q[idx]
  vals
}


#' Make rate names from param_map (VECTORIZED)
#' @keywords internal
make_rate_names <- function(param_map, state_labels, k) {
  n_params <- max(param_map)
  names_vec <- character(n_params)
  # Vectorized: find all off-diagonal cells with param_map > 0
  idx <- which(param_map > 0 & row(param_map) != col(param_map))
  idx_mat <- arrayInd(idx, dim(param_map))
  for (r in 1:nrow(idx_mat)) {
    i <- idx_mat[r, 1]
    j <- idx_mat[r, 2]
    p <- param_map[i, j]
    candidate <- paste0(state_labels[i], "->", state_labels[j])
    if (names_vec[p] == "") {
      names_vec[p] <- candidate
    } else if (!grepl(candidate, names_vec[p], fixed = TRUE)) {
      # Multiple transitions sharing same rate class
      names_vec[p] <- paste0(names_vec[p], "|", candidate)
    }
  }
  names_vec
}


#' Print method for ratescapeQ objects
#' @export
print.ratescapeQ <- function(x, ...) {
  cat("=== RateScape: Discrete Character Model ===\n")
  cat("Model:", x$model, " | Method:", toupper(x$method), "\n")
  cat("States:", paste(x$states, collapse = ", "),
      "(k =", length(x$states), ")\n")
  cat("Free parameters:", x$n_params, "\n\n")

  cat("Transition rate matrix (Q):\n")
  Q_print <- x$Q
  rownames(Q_print) <- colnames(Q_print) <- x$states
  print(round(Q_print, 5))
  cat("\n")

  if (x$method == "ml") {
    cat("Log-likelihood:", round(x$loglik, 3), "\n")
    cat("AIC:", round(x$AIC, 2), " | AICc:", round(x$AICc, 2),
        " | BIC:", round(x$BIC, 2), "\n")
  } else {
    cat("Posterior median log-likelihood:", round(x$loglik, 3), "\n")
    cat("\nRate posterior summaries:\n")
    print(x$rate_summary, digits = 4, row.names = FALSE)
  }

  if (!is.null(x$ancestral_states)) {
    cat("\nAncestral state probabilities available for",
        nrow(x$ancestral_states), "internal nodes.\n")
  }

  invisible(x)
}
