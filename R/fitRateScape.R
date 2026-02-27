#' Fit Branch-Specific Rate Variation Model
#'
#' Fits a model of discrete character evolution where each branch of the
#' phylogeny has a rate scalar relative to a background Q matrix. Uses a
#' spike-and-slab prior to distinguish background-rate branches from those
#' experiencing genuine rate shifts.
#'
#' The core likelihood computation and Gibbs sweep are implemented in C++
#' (via Rcpp/RcppArmadillo) for performance. A pure-R fallback is available
#' if the compiled code is not loaded.
#'
#' @param tree A phylo object (ultrametric or non-ultrametric)
#' @param data Named vector of character states (integer, 1-indexed).
#'   Names must match tree tip labels.
#' @param Q Either a k x k rate matrix (fixed Q, two-step approach) or
#'   a function that takes a parameter vector and returns a Q matrix
#'   (co-estimated Q, joint approach). See \code{\link{makeQ}} for helpers.
#' @param Q_params If Q is a function, a named list with elements:
#'   \code{init} (initial parameter vector), \code{lower} (lower bounds),
#'   \code{upper} (upper bounds), \code{proposal_sd} (proposal std devs).
#' @param ngen Number of MCMC generations. Default 10000.
#' @param sample_freq Sampling frequency. Default 10.
#' @param burnin Fraction of samples to discard as burn-in. Default 0.25.
#' @param prior_pi_a Shape1 parameter for Beta prior on pi (spike probability).
#'   Default 1 (uniform).
#' @param prior_pi_b Shape2 parameter for Beta prior on pi. Default 1.
#' @param prior_sigma2_rate Rate parameter for Exponential prior on sigma^2
#'   (slab variance). Default 1.
#' @param root_prior Root state prior: "equal", "stationary", or numeric vector.
#' @param use_cpp Logical. Use C++ backend? Default TRUE (falls back to R
#'   automatically if compiled code not available).
#' @param verbose Logical, print progress? Default TRUE.
#' @return An object of class "ratescape" containing posterior samples,
#'   rate estimates, and model fit information.
#' @export
#' @examples
#' \dontrun{
#' library(ape)
#' # Simple example with a 3-state equal-rates model
#' tree <- rtree(30)
#' Q <- matrix(0.1, 3, 3)
#' diag(Q) <- -0.2
#' data <- simRateScape(tree, Q, rates = rep(1, Nedge(tree)))
#' fit <- fitRateScape(tree, data, Q, ngen = 5000)
#' plotRateTree(fit)
#' }
fitRateScape <- function(tree, data, Q,
                         Q_params = NULL,
                         ngen = 10000,
                         sample_freq = 10,
                         burnin = 0.25,
                         prior_pi_a = 1,
                         prior_pi_b = 1,
                         prior_sigma2_rate = 1,
                         root_prior = "equal",
                         use_cpp = TRUE,
                         verbose = TRUE) {

  # --- Input validation ---
  if (!inherits(tree, "phylo")) stop("tree must be a phylo object")
  if (is.null(names(data))) stop("data must be a named vector")
  if (!all(tree$tip.label %in% names(data))) {
    stop("All tip labels must have corresponding entries in data")
  }

  co_estimate_Q <- is.function(Q)

  if (!co_estimate_Q) {
    if (!is.matrix(Q)) stop("Q must be a matrix or a function")
    nstates <- nrow(Q)
    Q_fixed <- Q
  } else {
    if (is.null(Q_params)) stop("Q_params required when Q is a function")
    Q_func <- Q
    q_current <- Q_params$init
    nstates <- nrow(Q_func(q_current))
    Q_fixed <- NULL
  }

  nedge <- ape::Nedge(tree)
  ntip <- ape::Ntip(tree)
  nnode <- ape::Nnode(tree)
  nsamples <- floor(ngen / sample_freq)

  # Ensure data is integer 1-indexed, ordered by tip label
  tip_states <- as.integer(data[tree$tip.label])
  names(tip_states) <- tree$tip.label

  if (any(tip_states < 1) || any(tip_states > nstates)) {
    stop("All states must be between 1 and ", nstates)
  }

  # --- Detect C++ availability ---
  cpp_available <- use_cpp && .check_cpp_available()
  if (use_cpp && !cpp_available) {
    if (verbose) message("C++ backend not available, using pure R (slower)")
  } else if (cpp_available && verbose) {
    message("Using C++ backend (Rcpp/RcppArmadillo)")
  }

  # --- Prepare tree for likelihood computation ---
  tree_po <- ape::reorder.phylo(tree, "postorder")

  # Root prior vector
  pi_root <- .resolve_root_prior(root_prior, nstates,
                                  if (co_estimate_Q) Q_func(q_current) else Q_fixed)

  # --- Initialize MCMC ---
  rates_current <- rep(1.0, nedge)
  z_current <- rep(0L, nedge)
  pi_current <- 0.8
  sigma2_current <- 0.5

  if (co_estimate_Q) {
    Q_current <- Q_func(q_current)
  } else {
    Q_current <- Q_fixed
  }

  # Compute initial log-likelihood
  llik_current <- .compute_llik(tree_po, tip_states, Q_current, rates_current,
                                 pi_root, ntip, nnode, cpp_available)

  # --- Storage ---
  rate_samples <- matrix(NA, nrow = nsamples, ncol = nedge)
  z_samples <- matrix(NA_integer_, nrow = nsamples, ncol = nedge)
  pi_samples <- numeric(nsamples)
  sigma2_samples <- numeric(nsamples)
  llik_samples <- numeric(nsamples)

  if (co_estimate_Q) {
    q_samples <- matrix(NA, nrow = nsamples, ncol = length(q_current))
  }

  rate_proposal_sd <- 0.3
  sigma2_proposal_sd <- 0.2
  n_rate_proposed <- 0
  n_rate_accepted <- 0
  n_sigma2_proposed <- 0
  n_sigma2_accepted <- 0
  sample_idx <- 0

  # --- MCMC loop ---
  if (verbose) cat("Running RateScape MCMC for", ngen, "generations\n")

  for (gen in 1:ngen) {

    # =====================================================
    # Steps 1 & 2: Gibbs sweep (z indicators) + MH (rates)
    # =====================================================
    if (cpp_available && !co_estimate_Q) {
      # === C++ FAST PATH ===
      # Entire Gibbs + MH sweep in one C++ call
      sweep_result <- gibbs_sweep_cpp(
        parent_vec = tree_po$edge[, 1],
        child_vec = tree_po$edge[, 2],
        edge_lengths = tree_po$edge.length,
        tip_states = tip_states,
        Q_mat = Q_current,
        rates = rates_current,
        z_indicators = z_current,
        pi_val = pi_current,
        sigma2 = sigma2_current,
        root_prior = pi_root,
        ntip = ntip,
        nnode = nnode,
        rate_proposal_sd = rate_proposal_sd
      )

      rates_current <- sweep_result$rates
      z_current <- sweep_result$z
      llik_current <- sweep_result$loglik
      n_rate_accepted <- n_rate_accepted + sweep_result$n_accepted
      n_rate_proposed <- n_rate_proposed + sweep_result$n_proposed

    } else {
      # === PURE R PATH ===
      # Step 1: Update z_i (Gibbs)
      for (e in 1:nedge) {
        saved_rate <- rates_current[e]

        # Spike likelihood
        rates_current[e] <- 1.0
        llik_spike <- .compute_llik(tree_po, tip_states, Q_current,
                                     rates_current, pi_root, ntip, nnode,
                                     cpp_available)

        # Slab likelihood
        if (z_current[e] == 1) {
          rates_current[e] <- saved_rate
          llik_slab <- .compute_llik(tree_po, tip_states, Q_current,
                                      rates_current, pi_root, ntip, nnode,
                                      cpp_available)
          slab_rate <- saved_rate
        } else {
          slab_rate <- exp(stats::rnorm(1, 0, sqrt(sigma2_current)))
          rates_current[e] <- slab_rate
          llik_slab <- .compute_llik(tree_po, tip_states, Q_current,
                                      rates_current, pi_root, ntip, nnode,
                                      cpp_available)
        }

        log_prior_slab <- stats::dlnorm(slab_rate, 0, sqrt(sigma2_current),
                                         log = TRUE)

        log_prob_spike <- log(pi_current) + llik_spike
        log_prob_slab <- log(1 - pi_current) + llik_slab + log_prior_slab

        max_log <- max(log_prob_spike, log_prob_slab)
        prob_spike <- exp(log_prob_spike - max_log) /
          (exp(log_prob_spike - max_log) + exp(log_prob_slab - max_log))

        if (stats::runif(1) < prob_spike) {
          z_current[e] <- 0L
          rates_current[e] <- 1.0
          llik_current <- llik_spike
        } else {
          z_current[e] <- 1L
          rates_current[e] <- slab_rate
          llik_current <- llik_slab
        }
      }

      # Step 2: Update r_i for slab branches (MH)
      slab_edges <- which(z_current == 1)
      if (length(slab_edges) > 0) {
        for (e in slab_edges) {
          n_rate_proposed <- n_rate_proposed + 1
          log_r_proposed <- stats::rnorm(1, log(rates_current[e]),
                                          rate_proposal_sd)
          r_proposed <- exp(log_r_proposed)

          rates_proposed <- rates_current
          rates_proposed[e] <- r_proposed
          llik_proposed <- .compute_llik(tree_po, tip_states, Q_current,
                                          rates_proposed, pi_root, ntip, nnode,
                                          cpp_available)

          log_prior_curr <- stats::dlnorm(rates_current[e], 0,
                                           sqrt(sigma2_current), log = TRUE)
          log_prior_prop <- stats::dlnorm(r_proposed, 0,
                                           sqrt(sigma2_current), log = TRUE)

          log_alpha <- (llik_proposed - llik_current) +
            (log_prior_prop - log_prior_curr)

          if (log(stats::runif(1)) < log_alpha) {
            rates_current[e] <- r_proposed
            llik_current <- llik_proposed
            n_rate_accepted <- n_rate_accepted + 1
          }
        }
      }
    }

    # =====================================================
    # Step 3: Update pi (Gibbs, Beta-Binomial conjugacy)
    # =====================================================
    n_spike <- sum(z_current == 0)
    n_slab <- sum(z_current == 1)
    pi_current <- stats::rbeta(1, prior_pi_a + n_spike, prior_pi_b + n_slab)

    # =====================================================
    # Step 4: Update sigma2 (MH)
    # =====================================================
    n_sigma2_proposed <- n_sigma2_proposed + 1
    log_sigma2_proposed <- stats::rnorm(1, log(sigma2_current), sigma2_proposal_sd)
    sigma2_proposed <- exp(log_sigma2_proposed)

    log_prior_s2_current <- -prior_sigma2_rate * sigma2_current
    log_prior_s2_proposed <- -prior_sigma2_rate * sigma2_proposed

    slab_edges <- which(z_current == 1)
    if (length(slab_edges) > 0) {
      log_lk_rates_current <- sum(stats::dlnorm(rates_current[slab_edges], 0,
                                                 sqrt(sigma2_current), log = TRUE))
      log_lk_rates_proposed <- sum(stats::dlnorm(rates_current[slab_edges], 0,
                                                  sqrt(sigma2_proposed), log = TRUE))
    } else {
      log_lk_rates_current <- 0
      log_lk_rates_proposed <- 0
    }

    log_alpha_s2 <- (log_lk_rates_proposed - log_lk_rates_current) +
      (log_prior_s2_proposed - log_prior_s2_current) +
      (log_sigma2_proposed - log(sigma2_current))

    if (log(stats::runif(1)) < log_alpha_s2) {
      sigma2_current <- sigma2_proposed
      n_sigma2_accepted <- n_sigma2_accepted + 1
    }

    # =====================================================
    # Step 5: Update Q parameters if co-estimating
    # =====================================================
    if (co_estimate_Q) {
      for (p in seq_along(q_current)) {
        q_proposed <- q_current
        q_proposed[p] <- stats::rnorm(1, q_current[p], Q_params$proposal_sd[p])

        if (q_proposed[p] >= Q_params$lower[p] &&
            q_proposed[p] <= Q_params$upper[p]) {

          Q_proposed <- Q_func(q_proposed)
          llik_proposed <- .compute_llik(tree_po, tip_states, Q_proposed,
                                          rates_current, pi_root, ntip, nnode,
                                          cpp_available)

          if (log(stats::runif(1)) < (llik_proposed - llik_current)) {
            q_current[p] <- q_proposed[p]
            Q_current <- Q_proposed
            llik_current <- llik_proposed
          }
        }
      }
    }

    # === Sample storage ===
    if (gen %% sample_freq == 0) {
      sample_idx <- sample_idx + 1
      rate_samples[sample_idx, ] <- rates_current
      z_samples[sample_idx, ] <- z_current
      pi_samples[sample_idx] <- pi_current
      sigma2_samples[sample_idx] <- sigma2_current
      llik_samples[sample_idx] <- llik_current
      if (co_estimate_Q) q_samples[sample_idx, ] <- q_current
    }

    # Progress
    if (verbose && gen %% max(1, ngen %/% 10) == 0) {
      pct <- round(100 * gen / ngen)
      cat(sprintf("  %d%% | logLik: %.2f | shifted: %d/%d | pi: %.3f | sigma2: %.3f\n",
                  pct, llik_current, n_slab, nedge, pi_current, sigma2_current))
    }
  }

  # --- Post-processing ---
  burnin_samples <- floor(nsamples * burnin)
  post_idx <- (burnin_samples + 1):nsamples

  post_rates <- rate_samples[post_idx, , drop = FALSE]
  post_z <- z_samples[post_idx, , drop = FALSE]
  post_pi <- pi_samples[post_idx]
  post_sigma2 <- sigma2_samples[post_idx]
  post_llik <- llik_samples[post_idx]

  rate_means <- colMeans(post_rates)
  rate_medians <- apply(post_rates, 2, stats::median)
  shift_probs <- colMeans(post_z)

  prior_odds <- mean(1 - post_pi) / mean(post_pi)
  posterior_odds <- shift_probs / pmax(1 - shift_probs, 1e-10)
  bayes_factors <- posterior_odds / prior_odds

  result <- list(
    tree = tree,
    data = data,
    Q = if (co_estimate_Q) Q_func(colMeans(q_samples[post_idx, , drop = FALSE])) else Q_fixed,
    nstates = nstates,
    rate_means = rate_means,
    rate_medians = rate_medians,
    shift_probs = shift_probs,
    bayes_factors = bayes_factors,
    pi_mean = mean(post_pi),
    sigma2_mean = mean(post_sigma2),
    rate_samples = post_rates,
    z_samples = post_z,
    pi_samples = post_pi,
    sigma2_samples = post_sigma2,
    llik_samples = post_llik,
    acceptance_rate = n_rate_accepted / max(n_rate_proposed, 1),
    sigma2_acceptance = n_sigma2_accepted / max(n_sigma2_proposed, 1),
    ngen = ngen,
    sample_freq = sample_freq,
    burnin = burnin,
    co_estimated_Q = co_estimate_Q,
    q_samples = if (co_estimate_Q) q_samples[post_idx, , drop = FALSE] else NULL,
    used_cpp = cpp_available
  )

  class(result) <- "ratescape"
  return(result)
}


# ============================================================
# Internal helper functions
# ============================================================

#' Check if C++ compiled code is available
#' @keywords internal
.check_cpp_available <- function() {
  tryCatch({
    is.loaded("_RateScape_pruning_llik_cpp")
  }, error = function(e) FALSE)
}

#' Resolve root prior specification to a numeric vector
#' @keywords internal
.resolve_root_prior <- function(root_prior, nstates, Q) {
  if (is.numeric(root_prior)) {
    if (length(root_prior) != nstates) {
      stop("root_prior vector must have length equal to number of states")
    }
    return(root_prior / sum(root_prior))
  }
  if (root_prior == "equal") {
    return(rep(1 / nstates, nstates))
  }
  if (root_prior == "stationary") {
    return(stationary_dist(Q))
  }
  stop("root_prior must be 'equal', 'stationary', or a numeric vector")
}

#' Compute log-likelihood dispatching to C++ or R
#' @keywords internal
.compute_llik <- function(tree, tip_states, Q, rates, root_prior_vec,
                           ntip, nnode, use_cpp) {
  if (use_cpp) {
    pruning_llik_cpp(
      parent_vec = tree$edge[, 1],
      child_vec = tree$edge[, 2],
      edge_lengths = tree$edge.length,
      tip_states = tip_states,
      Q_mat = Q,
      rates = rates,
      root_prior = root_prior_vec,
      ntip = ntip,
      nnode = nnode
    )
  } else {
    pruning_likelihood_scaled(tree, tip_states, Q, rates, root_prior_vec)
  }
}


#' Print method for ratescape objects
#' @param x A ratescape object
#' @param ... Additional arguments (ignored)
#' @export
print.ratescape <- function(x, ...) {
  cat("RateScape model fit\n")
  cat("-------------------\n")
  cat("Tree:", ape::Ntip(x$tree), "tips,", ape::Nedge(x$tree), "edges\n")
  cat("States:", x$nstates, "\n")
  cat("Q matrix:", ifelse(x$co_estimated_Q, "co-estimated", "fixed"), "\n")
  cat("Backend:", ifelse(isTRUE(x$used_cpp), "C++ (Rcpp)", "pure R"), "\n")
  cat("MCMC:", x$ngen, "generations, sampled every", x$sample_freq, "\n")
  cat("\nHyperparameters (posterior means):\n")
  cat("  pi (background prob):", round(x$pi_mean, 3), "\n")
  cat("  sigma2 (slab variance):", round(x$sigma2_mean, 3), "\n")
  n_shifted <- sum(x$shift_probs > 0.5)
  n_strong <- sum(x$bayes_factors > 10)
  cat("\nRate shifts:\n")
  cat("  Branches with P(shift) > 0.5:", n_shifted, "/",
      ape::Nedge(x$tree), "\n")
  cat("  Branches with BF > 10:", n_strong, "/",
      ape::Nedge(x$tree), "\n")
  cat("  Rate acceptance:", round(x$acceptance_rate, 3), "\n")
  cat("  Mean log-likelihood:", round(mean(x$llik_samples), 2), "\n")
  invisible(x)
}
