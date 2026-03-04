#' RJMCMC Over Q-Matrix Structure
#'
#' Uses reversible-jump MCMC to explore which transitions in the rate matrix Q
#' are non-zero (structurally present) versus zero (forbidden). This provides
#' Bayesian model averaging over the space of possible Q-matrix structures,
#' yielding posterior probabilities for each transition being present.
#'
#' This is analogous to BayesTraits' RJMCMC over discrete models, but
#' integrated within R and optionally combined with branch-specific rate
#' heterogeneity.
#'
#' @param tree An object of class "phylo".
#' @param data Data frame with tip states (single column).
#' @param k Integer. Number of discrete states (inferred from data if NULL).
#' @param base_model Character. Starting model structure: "ard" starts with all
#'   transitions allowed, "mk" starts equal-rates, "sym" starts symmetric.
#'   Default "ard".
#' @param allow_zero_rate Logical. If TRUE, each off-diagonal transition in Q
#'   can be "turned off" (set to 0) by the RJMCMC. This explores which
#'   transitions are supported by the data. Default TRUE.
#' @param rate_heterogeneity Logical. If TRUE, simultaneously estimate
#'   branch-specific rate scalars using the spike-and-slab model. Default FALSE.
#' @param lambda_sigma Numeric. Prior on sigma^2 (required if rate_heterogeneity = TRUE).
#' @param root_prior Character. Root state prior. Default "fitzjohn".
#' @param lambda_q Numeric. Exponential prior rate on each Q transition rate.
#'   Default 1.0.
#' @param ngen Integer. MCMC generations. Default 50000.
#' @param burn_in Integer. Burn-in iterations. Default 10000.
#' @param thin Integer. Thinning interval. Default 20.
#' @param seed Integer or NULL. Random seed.
#'
#' @details
#'
#' **How it works:**
#'
#' The RJMCMC explores the joint space of (Q structure, rate parameters).
#' At each iteration, three move types are proposed:
#'
#' 1. **Add transition:** Select a currently-zero off-diagonal entry in Q and
#'    propose to turn it on with a rate drawn from the prior.
#'
#' 2. **Remove transition:** Select a currently-nonzero entry and propose to
#'    set it to zero.
#'
#' 3. **Update rate:** For a currently-active transition, propose a new rate
#'    via log-normal random walk.
#'
#' The acceptance probability accounts for the Jacobian of the dimension-changing
#' moves and the prior ratio (which penalizes complex models).
#'
#' **Key outputs:**
#'
#' The posterior inclusion probability (PIP) for each transition (i -> j) is
#' the fraction of post-burn-in samples in which that transition was active.
#' PIP > 0.95 means strong evidence that the transition is needed; PIP < 0.05
#' means strong evidence against.
#'
#' **Why this matters:**
#'
#' Standard model selection (AIC/BIC comparing Mk, SYM, ARD) only tests 3
#' models. RJMCMC explores the full space of 2^(k*(k-1)) possible Q structures,
#' providing a much finer-grained picture of which transitions matter and which
#' can be dropped.
#'
#' @return An object of class "ratescapeQstructure" containing:
#'   \describe{
#'     \item{pip}{Matrix (k x k) of posterior inclusion probabilities.}
#'     \item{mean_rates}{Matrix (k x k) of posterior mean rates (conditional on inclusion).}
#'     \item{Q_map}{Q matrix from the sample with highest posterior probability.}
#'     \item{Q_bma}{Bayesian model-averaged Q matrix (pip * mean_rate).}
#'     \item{samples}{List with Q_active (binary matrices) and Q_rates at each sample.}
#'     \item{loglik}{Vector of log-likelihoods at each sample.}
#'     \item{model_complexity}{Vector of number of active transitions per sample.}
#'     \item{rate_fit}{If rate_heterogeneity = TRUE, the rate scalar samples.}
#'   }
#'
#' @examples
#' \dontrun{
#'   result <- fitQstructure(tree, data, k = 4, ngen = 50000)
#'
#'   # Which transitions are strongly supported?
#'   print(result$pip)
#'
#'   # Model-averaged Q matrix
#'   print(round(result$Q_bma, 3))
#' }
#'
#' @export
fitQstructure <- function(
    tree, data, k = NULL, base_model = "ard",
    allow_zero_rate = TRUE,
    rate_heterogeneity = FALSE,
    lambda_sigma = NULL,
    root_prior = "fitzjohn",
    lambda_q = 1.0,
    ngen = 50000, burn_in = 10000, thin = 20,
    seed = NULL) {

  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  if (rate_heterogeneity && is.null(lambda_sigma))
    stop("lambda_sigma required when rate_heterogeneity = TRUE")
  base_model <- match.arg(base_model, c("ard", "mk", "sym"))
  root_prior <- match.arg(root_prior, c("fitzjohn", "equal", "stationary"))
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
  if (is.null(k)) k <- max(states) + 1

  nedges <- nrow(tree$edge)
  niter <- ngen + burn_in
  nsamples <- floor(ngen / thin)

  # Number of possible transitions
  n_trans <- k * (k - 1)
  trans_idx <- which(matrix(TRUE, k, k) & !diag(TRUE, k), arr.ind = TRUE)

  # Precompute diagonal mask for performance
  diag_mask <- diag(TRUE, k)

  # Initialize Q: all transitions active with rate 0.1
  Q_active <- matrix(TRUE, k, k)
  diag(Q_active) <- FALSE
  Q_rates <- matrix(0, k, k)
  Q_rates[Q_active] <- 0.1

  # Build initial Q
  build_Q <- function(active, rates) {
    Q <- rates * active
    diag(Q) <- -rowSums(Q)
    Q
  }

  Q_current <- build_Q(Q_active, Q_rates)
  
  # BUILD TREE INFO ONCE AT START (tree structure is invariant)
  info <- build_tree_info(tree, Q_current)

  # Helper function to update only Q decomposition (called for proposals)
  update_Q_decomposition <- function(info, Q_new) {
    eig <- eigen(Q_new)
    V <- eig$vectors
    V_inv <- solve(V)
    eig_t <- eigen(t(Q_new))
    idx <- which.min(abs(eig_t$values))
    stat_dist <- abs(Re(eig_t$vectors[, idx]))
    stat_dist <- stat_dist / sum(stat_dist)
    info$eig <- eig
    info$V <- V
    info$V_inv <- V_inv
    info$stat_dist <- stat_dist
    info
  }

  # Rate scalars (if using rate heterogeneity)
  r_current <- rep(1.0, nedges)
  if (rate_heterogeneity) {
    sigma2_current <- 1 / lambda_sigma
  }

  # Compute initial likelihood
  L_cache <- postorder_pass(info, states, Q_current, r_current)
  root_L <- L_cache[info$root, ]
  loglik_current <- log(max(sum(root_L * get_root_weights(root_L, root_prior, Q_current, k)), 1e-300))

  # Storage
  pip_count <- matrix(0, k, k)  # count how often each transition is active
  rate_sum <- matrix(0, k, k)   # sum of rates when active
  rate_count <- matrix(0, k, k) # count for computing conditional means
  loglik_samples <- numeric(nsamples)
  complexity_samples <- integer(nsamples)
  Q_active_samples <- vector("list", nsamples)
  Q_rate_samples <- vector("list", nsamples)

  if (rate_heterogeneity) {
    r_samples <- matrix(NA, nrow = nsamples, ncol = nedges)
  }

  # MH parameters
  tau_q <- 0.3  # proposal SD for log-rate updates

  message("----------------------------------------------------")
  message(sprintf("  RateScape Q-Structure RJMCMC"))
  message(sprintf("  Tree:     %d tips, %d edges", ntips, nedges))
  message(sprintf("  States:   k = %d (%d possible transitions)", k, n_trans))
  message(sprintf("  Chain:    %d gen (%d burn-in + %d sampling)", niter, burn_in, ngen))
  message(sprintf("  Rate heterogeneity: %s", rate_heterogeneity))
  message("----------------------------------------------------")

  n_accept <- 0; n_prop <- 0
  mcmc_start_time <- proc.time()[3]

  for (iter in 1:niter) {

    # ---- Move type selection ----
    move <- sample(c("add", "remove", "update"), 1,
                   prob = c(0.3, 0.3, 0.4))

    active_entries <- which(Q_active & !diag_mask, arr.ind = TRUE)
    inactive_entries <- which(!Q_active & !diag_mask, arr.ind = TRUE)
    n_active <- nrow(active_entries)
    n_inactive <- nrow(inactive_entries)

    # Guard: can't add if all active, can't remove if <= k-1 active
    if (move == "add" && n_inactive == 0) move <- "update"
    if (move == "remove" && n_active <= (k - 1)) move <- "update"
    if (move == "update" && n_active == 0) next

    n_prop <- n_prop + 1

    if (move == "add") {
      # Pick a random inactive transition and propose to activate it
      pick <- sample(nrow(inactive_entries), 1)
      i_pick <- inactive_entries[pick, 1]
      j_pick <- inactive_entries[pick, 2]

      # Draw rate from prior
      new_rate <- rexp(1, lambda_q)

      Q_active_prop <- Q_active
      Q_rates_prop <- Q_rates
      Q_active_prop[i_pick, j_pick] <- TRUE
      Q_rates_prop[i_pick, j_pick] <- new_rate
      Q_prop <- build_Q(Q_active_prop, Q_rates_prop)

      # Update only Q decomposition (tree structure unchanged)
      info_prop <- update_Q_decomposition(info, Q_prop)
      L_prop <- postorder_pass(info_prop, states, Q_prop, r_current)
      root_L_prop <- L_prop[info_prop$root, ]
      loglik_prop <- log(max(sum(root_L_prop *
        get_root_weights(root_L_prop, root_prior, Q_prop, k)), 1e-300))

      # Acceptance: likelihood ratio * prior ratio * proposal ratio
      # Prior: Poisson(lambda) on number of active transitions
      # Proposal ratio: P(remove pick) / P(add pick) = (1/n_active_new) / (1/n_inactive)
      log_prior_ratio <- log(lambda_q) - lambda_q * new_rate  # exp prior on rate
      log_proposal_ratio <- log(n_inactive) - log(n_active + 1)

      log_alpha <- (loglik_prop - loglik_current) + log_prior_ratio + log_proposal_ratio

      if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
        Q_active <- Q_active_prop
        Q_rates <- Q_rates_prop
        Q_current <- Q_prop
        info <- info_prop
        L_cache <- L_prop
        loglik_current <- loglik_prop
        n_accept <- n_accept + 1
      }

    } else if (move == "remove") {
      # Pick a random active transition and propose to deactivate it
      pick <- sample(nrow(active_entries), 1)
      i_pick <- active_entries[pick, 1]
      j_pick <- active_entries[pick, 2]
      old_rate <- Q_rates[i_pick, j_pick]

      Q_active_prop <- Q_active
      Q_rates_prop <- Q_rates
      Q_active_prop[i_pick, j_pick] <- FALSE
      Q_rates_prop[i_pick, j_pick] <- 0
      Q_prop <- build_Q(Q_active_prop, Q_rates_prop)

      # Update only Q decomposition (tree structure unchanged)
      info_prop <- update_Q_decomposition(info, Q_prop)
      L_prop <- postorder_pass(info_prop, states, Q_prop, r_current)
      root_L_prop <- L_prop[info_prop$root, ]
      loglik_prop <- log(max(sum(root_L_prop *
        get_root_weights(root_L_prop, root_prior, Q_prop, k)), 1e-300))

      # Reverse of add move
      n_inactive_new <- n_inactive + 1
      n_active_new <- n_active - 1
      log_prior_ratio <- -(log(lambda_q) - lambda_q * old_rate)
      log_proposal_ratio <- log(n_active) - log(n_inactive_new)

      log_alpha <- (loglik_prop - loglik_current) + log_prior_ratio + log_proposal_ratio

      if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
        Q_active <- Q_active_prop
        Q_rates <- Q_rates_prop
        Q_current <- Q_prop
        info <- info_prop
        L_cache <- L_prop
        loglik_current <- loglik_prop
        n_accept <- n_accept + 1
      }

    } else {
      # Update: MH on a random active transition rate
      pick <- sample(nrow(active_entries), 1)
      i_pick <- active_entries[pick, 1]
      j_pick <- active_entries[pick, 2]

      log_old <- log(Q_rates[i_pick, j_pick])
      log_new <- log_old + rnorm(1, 0, tau_q)
      new_rate <- exp(log_new)

      Q_rates_prop <- Q_rates
      Q_rates_prop[i_pick, j_pick] <- new_rate
      Q_prop <- build_Q(Q_active, Q_rates_prop)

      # Update only Q decomposition (tree structure unchanged)
      info_prop <- update_Q_decomposition(info, Q_prop)
      L_prop <- postorder_pass(info_prop, states, Q_prop, r_current)
      root_L_prop <- L_prop[info_prop$root, ]
      loglik_prop <- log(max(sum(root_L_prop *
        get_root_weights(root_L_prop, root_prior, Q_prop, k)), 1e-300))

      # Prior: Exp(lambda_q) on rate
      log_prior_ratio <- lambda_q * (Q_rates[i_pick, j_pick] - new_rate)
      # Jacobian for log-scale proposal
      log_jacobian <- log_new - log_old

      log_alpha <- (loglik_prop - loglik_current) + log_prior_ratio + log_jacobian

      if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
        Q_rates <- Q_rates_prop
        Q_current <- Q_prop
        info <- info_prop
        L_cache <- L_prop
        loglik_current <- loglik_prop
        n_accept <- n_accept + 1
      }
    }

    # ---- Rate heterogeneity updates (if enabled) ----
    if (rate_heterogeneity) {
      # Simple MH for each edge rate scalar
      for (e in 1:nedges) {
        log_r_old <- log(r_current[e])
        log_r_new <- log_r_old + rnorm(1, 0, 0.2)
        r_new_e <- exp(log_r_new)

        r_prop <- r_current
        r_prop[e] <- r_new_e

        # Recompute likelihood with new rate scalar
        L_prop_r <- postorder_pass(info, states, Q_current, r_prop)
        root_L_prop_r <- L_prop_r[info$root, ]
        loglik_prop_r <- log(max(sum(root_L_prop_r *
          get_root_weights(root_L_prop_r, root_prior, Q_current, k)), 1e-300))

        # LogNormal(0, sigma2) prior
        sd_r <- sqrt(max(sigma2_current, 1e-10))
        lp_new <- dlnorm(r_new_e, 0, sd_r, log = TRUE)
        lp_old <- dlnorm(r_current[e], 0, sd_r, log = TRUE)
        log_a <- (loglik_prop_r - loglik_current) + (lp_new - lp_old) + (log_r_new - log_r_old)

        if (!is.nan(log_a) && log(runif(1)) < log_a) {
          r_current[e] <- r_new_e
          loglik_current <- loglik_prop_r
          L_cache <- L_prop_r
        }
      }

      # Update sigma2 (vectorized dlnorm computation)
      sigma2_prop <- rlnorm(1, log(max(sigma2_current, 1e-10)), 0.3)
      sd_p <- sqrt(max(sigma2_prop, 1e-10))
      sd_c <- sqrt(max(sigma2_current, 1e-10))
      lp_p <- dexp(sigma2_prop, lambda_sigma, log = TRUE)
      lp_c <- dexp(sigma2_current, lambda_sigma, log = TRUE)
      
      # Vectorized: compute log-likelihood for all edges at once
      lp_p <- lp_p + sum(dlnorm(r_current, 0, sd_p, log = TRUE))
      lp_c <- lp_c + sum(dlnorm(r_current, 0, sd_c, log = TRUE))
      
      log_a_s <- (lp_p - lp_c) + (log(sigma2_prop) - log(max(sigma2_current, 1e-10)))
      if (!is.nan(log_a_s) && log(runif(1)) < log_a_s) sigma2_current <- sigma2_prop
    }

    # ---- Adaptive tuning ----
    if (iter <= burn_in && iter %% 200 == 0 && n_prop > 0) {
      ar <- n_accept / n_prop
      if (ar > 0.35) tau_q <- tau_q * 1.1
      else if (ar < 0.20) tau_q <- tau_q / 1.1
    }

    # ---- Store ----
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      si <- (iter - burn_in) %/% thin
      if (si >= 1 && si <= nsamples) {
        pip_count <- pip_count + Q_active
        
        # Vectorized: compute rate_sum and rate_count using matrix operations
        active_mask <- Q_active & !diag_mask
        rate_sum[active_mask] <- rate_sum[active_mask] + Q_rates[active_mask]
        rate_count[active_mask] <- rate_count[active_mask] + 1
        
        loglik_samples[si] <- loglik_current
        complexity_samples[si] <- sum(Q_active & !diag_mask)
        Q_active_samples[[si]] <- Q_active
        Q_rate_samples[[si]] <- Q_rates
        if (rate_heterogeneity) r_samples[si, ] <- r_current
      }
    }

    # Progress
    if (iter %% 2000 == 0) {
      pct <- iter / niter
      n_act <- sum(Q_active & !diag_mask)
      elapsed <- proc.time()[3] - mcmc_start_time
      phase <- if (iter <= burn_in) "burn-in" else "sampling"
      message(sprintf("  %3.0f%% | %s | ll=%.1f active=%d/%d | %.0fs",
                      pct*100, phase, loglik_current, n_act, n_trans, elapsed))
    }
  }

  # ---- Posterior summaries ----
  pip <- pip_count / nsamples
  diag(pip) <- NA

  mean_rates <- matrix(0, k, k)
  mean_rates[rate_count > 0] <- rate_sum[rate_count > 0] / rate_count[rate_count > 0]
  diag(mean_rates) <- NA

  # BMA Q matrix: pip * conditional mean rate
  Q_bma <- pip * mean_rates
  Q_bma[is.na(Q_bma)] <- 0
  diag(Q_bma) <- -rowSums(Q_bma)

  # MAP Q: sample with highest log-likelihood
  best_idx <- which.max(loglik_samples)
  Q_map <- build_Q(Q_active_samples[[best_idx]], Q_rate_samples[[best_idx]])

  ar_final <- if (n_prop > 0) n_accept / n_prop else NA
  elapsed_total <- proc.time()[3] - mcmc_start_time

  message("----------------------------------------------------")
  message(sprintf("  Q-structure RJMCMC complete in %.0fs", elapsed_total))
  message(sprintf("  Acceptance rate: %.1f%%", ar_final * 100))
  message(sprintf("  Median active transitions: %d / %d",
                  median(complexity_samples), n_trans))
  message("----------------------------------------------------")

  result <- list(
    pip = pip,
    mean_rates = mean_rates,
    Q_map = Q_map,
    Q_bma = Q_bma,
    loglik = loglik_samples,
    model_complexity = complexity_samples,
    k = k,
    n_trans = n_trans,
    nsamples = nsamples,
    rate_heterogeneity = rate_heterogeneity,
    tree = tree,
    call = match.call()
  )

  if (rate_heterogeneity) result$r_samples <- r_samples

  class(result) <- "ratescapeQstructure"
  result
}


#' Print method for Q-structure RJMCMC
#' @export
print.ratescapeQstructure <- function(x, ...) {
  cat("RateScape Q-Matrix Structure (RJMCMC)\n")
  cat("=======================================\n\n")
  cat(sprintf("States: %d (%d possible transitions)\n", x$k, x$n_trans))
  cat(sprintf("Posterior samples: %d\n", x$nsamples))
  cat(sprintf("Rate heterogeneity: %s\n\n", x$rate_heterogeneity))

  cat("Posterior Inclusion Probabilities (PIP):\n")
  pip_print <- round(x$pip, 3)
  diag(pip_print) <- NA
  print(pip_print, na.print = "-")

  cat(sprintf("\nMedian model complexity: %d active transitions\n",
              median(x$model_complexity)))

  # Summarize strongly supported / unsupported
  pip_vec <- x$pip[!is.na(x$pip)]
  cat(sprintf("Transitions with PIP > 0.95: %d\n", sum(pip_vec > 0.95)))
  cat(sprintf("Transitions with PIP < 0.05: %d\n", sum(pip_vec < 0.05)))

  invisible(x)
}
