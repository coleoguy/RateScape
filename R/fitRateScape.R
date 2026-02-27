#' Bayesian Spike-and-Slab Fitting for Rate Heterogeneity
#'
#' Joint (z_i, r_i) reversible-jump MCMC with incremental cache updates.
#'
#' @param tree phylo object.
#' @param data Data frame with tip states.
#' @param Q Rate matrix.
#' @param estimate_Q Co-estimate Q? Default FALSE.
#' @param lambda_sigma Exp prior rate on sigma^2. NO DEFAULT.
#' @param expected_background Expected background proportion.
#' @param root_prior "fitzjohn", "equal", or "stationary".
#' @param q_prior "diffuse", "moderate", "informative", or numeric.
#' @param ngen MCMC generations.
#' @param burn_in Burn-in iterations.
#' @param thin Thinning interval.
#' @param tau_init Initial MH proposal SD.
#' @param target_acceptance Target MH acceptance rate.
#' @param seed Random seed.
#'
#' @return ratescapeFit object.
#' @export
fitRateScape <- function(
    tree, data, Q = NULL, estimate_Q = FALSE,
    lambda_sigma,
    expected_background = NULL,
    root_prior = "fitzjohn",
    q_prior = "diffuse",
    ngen = 10000, burn_in = 2000, thin = 10,
    tau_init = 0.3, target_acceptance = 0.30,
    seed = NULL) {

  if (missing(lambda_sigma)) stop("lambda_sigma required (no default).")
  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  if (!is.data.frame(data) && !is.matrix(data)) stop("data must be data.frame/matrix")

  ntips <- length(tree$tip.label)
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

  if (!is.null(seed)) set.seed(seed)
  if (is.null(Q)) {
    if (!estimate_Q) stop("Q is NULL but estimate_Q = FALSE")
    Q <- makeQ(model = "mk", k = k)
  }
  root_prior <- match.arg(root_prior, c("fitzjohn", "equal", "stationary"))

  if (is.numeric(q_prior)) { lambda_Q <- q_prior
  } else {
    lambda_Q <- switch(match.arg(q_prior, c("diffuse","moderate","informative")),
                       diffuse=0.1, moderate=1.0, informative=10.0)
  }

  if (!is.null(expected_background)) {
    kappa <- 10; a_pi <- expected_background*kappa; b_pi <- (1-expected_background)*kappa
  } else { a_pi <- 1; b_pi <- 1 }

  nedges <- nrow(tree$edge)
  niter <- ngen + burn_in
  nsamples <- floor(ngen / thin)

  info <- build_tree_info(tree, Q)

  # Initialize
  r_current <- rep(1.0, nedges)
  z_current <- rep(1L, nedges)
  pi_current <- 0.5
  sigma2_current <- 1 / lambda_sigma

  r_samples <- matrix(NA, nrow = nsamples, ncol = nedges)
  z_samples <- matrix(NA, nrow = nsamples, ncol = nedges)
  pi_samples <- numeric(nsamples)
  sigma2_samples <- numeric(nsamples)
  loglik_samples <- numeric(nsamples)

  n_accept_r <- 0; n_prop_r <- 0
  n_accept_rj <- 0; n_prop_rj <- 0
  tau_current <- tau_init

  message("──────────────────────────────────────────────────────")
  message(sprintf("  RateScape MCMC"))
  message(sprintf("  Tree:     %d tips, %d edges", ntips, nedges))
  message(sprintf("  States:   k = %d", nrow(Q)))
  message(sprintf("  Chain:    %d generations (%d burn-in + %d sampling)",
                  niter, burn_in, ngen))
  message(sprintf("  Thinning: every %d → %d posterior samples", thin, nsamples))
  message("──────────────────────────────────────────────────────")

  # Initial likelihood + cache
  L_cache <- postorder_pass(info, states, Q, r_current)
  root_L <- L_cache[info$root, ]
  loglik_current <- log(max(sum(root_L * get_root_weights(root_L, root_prior, Q, k)), 1e-300))

  message(sprintf("  Initial log-likelihood: %.2f", loglik_current))
  message("  Starting burn-in...")
  mcmc_start_time <- proc.time()[3]
  last_report_time <- mcmc_start_time

  for (iter in 1:niter) {
    sd_slab <- sqrt(max(sigma2_current, 1e-10))

    # === JOINT (z_i, r_i) reversible-jump ===
    for (edge_idx in 1:nedges) {
      n_prop_rj <- n_prop_rj + 1

      if (z_current[edge_idx] == 1) {
        # Spike -> Slab: draw r from prior
        r_new <- rlnorm(1, 0, sd_slab)
        res <- compute_likelihood_one_edge(
          info, states, Q, r_current, L_cache, edge_idx, r_new, root_prior, return_L = TRUE)
        log_alpha <- (res$loglik - loglik_current) +
          log(1 - pi_current + 1e-300) - log(pi_current + 1e-300)

        if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
          z_current[edge_idx] <- 0L
          r_current[edge_idx] <- r_new
          loglik_current <- res$loglik
          L_cache <- res$L  # Incremental update!
          n_accept_rj <- n_accept_rj + 1

          # Inner MH refinement: let r explore before next RJ can flip back
          for (.mh in 1:2) {
            log_r_old <- log(r_current[edge_idx])
            log_r_prop <- log_r_old + rnorm(1, 0, tau_current)
            r_prop <- exp(log_r_prop)
            res_mh <- compute_likelihood_one_edge(
              info, states, Q, r_current, L_cache, edge_idx, r_prop, root_prior, return_L = TRUE)
            lp_new <- dlnorm(r_prop, 0, sd_slab, log = TRUE)
            lp_old <- dlnorm(r_current[edge_idx], 0, sd_slab, log = TRUE)
            la <- (res_mh$loglik - loglik_current) + (lp_new - lp_old) + (log_r_prop - log_r_old)
            if (!is.nan(la) && log(runif(1)) < la) {
              r_current[edge_idx] <- r_prop
              loglik_current <- res_mh$loglik
              L_cache <- res_mh$L
            }
          }
        }
      } else {
        # Slab -> Spike: set r=1
        res <- compute_likelihood_one_edge(
          info, states, Q, r_current, L_cache, edge_idx, 1.0, root_prior, return_L = TRUE)
        log_alpha <- (res$loglik - loglik_current) +
          log(pi_current + 1e-300) - log(1 - pi_current + 1e-300)

        if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
          z_current[edge_idx] <- 1L
          r_current[edge_idx] <- 1.0
          loglik_current <- res$loglik
          L_cache <- res$L
          n_accept_rj <- n_accept_rj + 1
        }
      }
    }

    # === MH for r_i in slab ===
    slab_edges <- which(z_current == 0)
    for (edge_idx in slab_edges) {
      n_prop_r <- n_prop_r + 1

      log_r_old <- log(r_current[edge_idx])
      log_r_new <- log_r_old + rnorm(1, 0, tau_current)
      r_new <- exp(log_r_new)

      res <- compute_likelihood_one_edge(
        info, states, Q, r_current, L_cache, edge_idx, r_new, root_prior, return_L = TRUE)

      lp_new <- dlnorm(r_new, 0, sd_slab, log = TRUE)
      lp_old <- dlnorm(r_current[edge_idx], 0, sd_slab, log = TRUE)
      log_alpha <- (res$loglik - loglik_current) + (lp_new - lp_old) + (log_r_new - log_r_old)

      if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
        r_current[edge_idx] <- r_new
        loglik_current <- res$loglik
        L_cache <- res$L
        n_accept_r <- n_accept_r + 1
      }
    }

    # === Gibbs for pi ===
    n_bg <- sum(z_current == 1)
    pi_current <- rbeta(1, a_pi + n_bg, b_pi + nedges - n_bg)

    # === MH for sigma^2 ===
    sigma2_prop <- rlnorm(1, log(max(sigma2_current, 1e-10)), 0.5)
    lp_p <- dexp(sigma2_prop, lambda_sigma, log = TRUE)
    lp_c <- dexp(sigma2_current, lambda_sigma, log = TRUE)
    sd_p <- sqrt(max(sigma2_prop, 1e-10))
    for (e in slab_edges) {
      lp_p <- lp_p + dlnorm(r_current[e], 0, sd_p, log = TRUE)
      lp_c <- lp_c + dlnorm(r_current[e], 0, sd_slab, log = TRUE)
    }
    log_a_s <- (lp_p - lp_c) + (log(sigma2_prop) - log(max(sigma2_current, 1e-10)))
    if (!is.nan(log_a_s) && log(runif(1)) < log_a_s) sigma2_current <- sigma2_prop

    # === Adaptive tuning ===
    if (iter <= burn_in && iter %% 100 == 0 && n_prop_r > 0) {
      ar <- n_accept_r / n_prop_r
      if (ar > target_acceptance + 0.05) tau_current <- tau_current * 1.1
      else if (ar < target_acceptance - 0.05) tau_current <- tau_current / 1.1
      n_accept_r <- 0; n_prop_r <- 0
    }

    # === Store ===
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      si <- (iter - burn_in) %/% thin
      if (si >= 1 && si <= nsamples) {
        r_samples[si, ] <- r_current
        z_samples[si, ] <- z_current
        pi_samples[si] <- pi_current
        sigma2_samples[si] <- sigma2_current
        loglik_samples[si] <- loglik_current
      }
    }

    # === Progress reporting ===
    if (iter %% 500 == 0) {
      now <- proc.time()[3]
      elapsed <- now - mcmc_start_time
      pct <- iter / niter
      eta <- elapsed / pct - elapsed  # estimated time remaining
      phase <- if (iter <= burn_in) "burn-in" else "sampling"

      # Progress bar: 30 characters wide
      bar_width <- 30
      filled <- round(pct * bar_width)
      bar <- paste0("[",
                    paste(rep("=", filled), collapse = ""),
                    if (filled < bar_width) ">",
                    paste(rep(" ", max(0, bar_width - filled - 1)), collapse = ""),
                    "]")

      # Format time nicely
      fmt_time <- function(secs) {
        if (secs < 60) return(sprintf("%.0fs", secs))
        if (secs < 3600) return(sprintf("%.0fm %02.0fs", secs %/% 60, secs %% 60))
        return(sprintf("%.0fh %02.0fm", secs %/% 3600, (secs %% 3600) %/% 60))
      }

      n_slab <- sum(z_current == 0)
      message(sprintf("\r  %s %3.0f%% | %s | ll=%.1f pi=%.3f s2=%.2f slab=%d/%d | %s elapsed, ~%s left",
                      bar, pct * 100, phase,
                      loglik_current, pi_current, sigma2_current,
                      n_slab, nedges,
                      fmt_time(elapsed), fmt_time(eta)),
              appendLF = TRUE)

      # Announce phase transition
      if (iter == burn_in) {
        ar_burn <- if (n_prop_r > 0) n_accept_r / n_prop_r else NA
        ar_rj_burn <- if (n_prop_rj > 0) n_accept_rj / n_prop_rj else NA
        message(sprintf("  ── Burn-in complete (tau=%.3f, MH accept=%.1f%%, RJ accept=%.1f%%) ──",
                        tau_current,
                        if (!is.na(ar_burn)) ar_burn * 100 else 0,
                        if (!is.na(ar_rj_burn)) ar_rj_burn * 100 else 0))
        message("  Starting posterior sampling...")
      }
    }
  }

  # ── Completion summary ──
  mcmc_elapsed <- proc.time()[3] - mcmc_start_time
  fmt_time <- function(secs) {
    if (secs < 60) return(sprintf("%.1fs", secs))
    if (secs < 3600) return(sprintf("%.0fm %02.0fs", secs %/% 60, secs %% 60))
    return(sprintf("%.0fh %02.0fm", secs %/% 3600, (secs %% 3600) %/% 60))
  }

  ar_r <- if (n_prop_r > 0) n_accept_r/n_prop_r else NA
  ar_rj <- if (n_prop_rj > 0) n_accept_rj/n_prop_rj else NA

  message("──────────────────────────────────────────────────────")
  message(sprintf("  MCMC complete in %s", fmt_time(mcmc_elapsed)))
  message(sprintf("  Final log-likelihood: %.2f", loglik_current))
  message(sprintf("  MH acceptance rate:   %s",
                  if (!is.na(ar_r)) sprintf("%.1f%%", ar_r * 100) else "N/A"))
  message(sprintf("  RJ acceptance rate:   %s",
                  if (!is.na(ar_rj)) sprintf("%.1f%%", ar_rj * 100) else "N/A"))
  message(sprintf("  Final tau:            %.4f", tau_current))
  message(sprintf("  Slab edges:           %d / %d (%.0f%%)",
                  sum(z_current == 0), nedges, 100 * sum(z_current == 0) / nedges))
  message(sprintf("  Posterior samples:    %d", nsamples))
  message("──────────────────────────────────────────────────────")

  # Post-hoc: mean(r) = 1
  for (s in 1:nsamples) {
    rm <- mean(r_samples[s, ]); if (rm > 0) r_samples[s, ] <- r_samples[s, ] / rm
  }

  pi_prior <- a_pi/(a_pi+b_pi); pi_post <- mean(pi_samples, na.rm=TRUE)
  bf <- ((1-pi_post)/(pi_post+1e-10)) / ((1-pi_prior)/(pi_prior+1e-10))

  result <- list(
    mcmc_samples = list(r=r_samples, z=z_samples, pi=pi_samples,
                        sigma2=sigma2_samples, loglik=loglik_samples),
    tree=tree, data=states, Q_fixed=Q, call=match.call(),
    acceptance_rate=ar_r, rj_acceptance_rate=ar_rj, tau_final=tau_current,
    prior_settings=list(lambda_sigma=lambda_sigma, expected_background=expected_background,
                        a_pi=a_pi, b_pi=b_pi, lambda_Q=lambda_Q, root_prior=root_prior),
    bayes_factor=bf)
  class(result) <- "ratescapeFit"
  result
}
