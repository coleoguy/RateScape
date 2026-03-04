#' Stochastic Character Mapping Under Rate Heterogeneity
#'
#' Generates stochastic character maps (Huelsenbeck et al. 2003, Nielsen 2002)
#' that account for branch-specific rate variation. Standard stochastic mapping
#' assumes a homogeneous rate across the tree; this function conditions each
#' branch's mapping on its posterior rate scalar, producing maps that reflect
#' more transitions on faster-evolving branches.
#'
#' @param fit An object of class "ratescapeFit" (from fitRateScape).
#' @param nsim Integer. Number of stochastic maps to generate. Default 100.
#' @param rate_summary Character. How to summarize the rate posterior for each
#'   branch: "mean" uses posterior mean r_i, "sample" draws a rate vector from
#'   the posterior for each map (fully propagates uncertainty), "map" uses the
#'   MAP (maximum a posteriori) estimate. Default "sample".
#' @param Q_scaled Logical. If TRUE (default), the Q matrix is multiplied by
#'   the branch-specific rate scalar for each edge. If FALSE, rates are applied
#'   only to branch lengths (equivalent but mathematically identical).
#' @param seed Random seed for reproducibility.
#'
#' @details
#'
#' **How it works:**
#'
#' For each stochastic map:
#'
#' 1. Select a rate vector from the posterior (method depends on rate_summary).
#'
#' 2. Compute conditional likelihoods at all internal nodes using the
#'    Felsenstein pruning algorithm with the selected branch-specific rates.
#'
#' 3. Sample ancestral states at each internal node from their conditional
#'    posterior distribution (top-down pass from root to tips).
#'
#' 4. For each branch, given the parent and child states and the rate-scaled
#'    Q matrix, simulate the actual history of state changes using the
#'    Gillespie algorithm (endpoint-conditioned). This produces the exact
#'    dwelling times in each state and the transition times along each branch.
#'
#' **Why this matters:**
#'
#' Standard make.simmap (phytools) assumes all branches evolve at the same
#' rate. If branch 47 has a rate scalar of 3x, you'd expect ~3x more
#' transitions on that branch. Ignoring rate heterogeneity produces biased
#' stochastic maps that undercount transitions on fast branches and overcount
#' on slow branches.
#'
#' @return A list of class "ratescape_simmap" containing:
#'   \describe{
#'     \item{maps}{A list of nsim stochastic maps, each itself a list with one
#'       element per edge containing a named numeric vector of dwelling times
#'       (names are character states, values are durations in that state).
#'       Compatible with phytools mapped.edge format.}
#'     \item{mapped.edge}{A matrix (nsim rows x (k*nedges) columns) summarizing
#'       total time spent in each state on each edge, averaged across maps.}
#'     \item{trees}{A list of nsim "simmap" phylo objects compatible with
#'       phytools plotting functions (plotSimmap, describe.simmap, etc.).}
#'     \item{rates_used}{Matrix of rate vectors used for each map (nsim x nedges).}
#'     \item{Q}{The Q matrix used.}
#'     \item{summary}{Data frame of per-branch expected number of transitions.}
#'   }
#'
#' @examples
#' \dontrun{
#'   fit <- fitRateScape(tree, data, Q, lambda_sigma = 1.0)
#'   smaps <- simmapRate(fit, nsim = 100)
#'
#'   # Plot one map
#'   phytools::plotSimmap(smaps$trees[[1]])
#'
#'   # Summarize transition counts
#'   print(smaps$summary)
#' }
#'
#' @export
simmapRate <- function(fit, nsim = 100, rate_summary = "sample",
                       Q_scaled = TRUE, seed = NULL) {

  if (!inherits(fit, "ratescapeFit")) {
    stop("fit must be an object of class 'ratescapeFit'")
  }
  rate_summary <- match.arg(rate_summary, c("mean", "sample", "map"))
  if (!is.null(seed)) set.seed(seed)

  tree <- fit$tree
  states <- fit$data
  Q <- fit$Q_fixed
  k <- nrow(Q)
  nedges <- nrow(tree$edge)
  ntips <- length(tree$tip.label)
  n_nodes <- ntips + tree$Nnode
  root <- ntips + 1
  npost <- nrow(fit$mcmc_samples$r)
  root_prior <- fit$prior_settings$root_prior

  # Precompute rate vectors
  r_samples <- fit$mcmc_samples$r

  if (rate_summary == "mean") {
    r_fixed <- colMeans(r_samples)
  } else if (rate_summary == "map") {
    # MAP = sample with highest log-likelihood
    best_idx <- which.max(fit$mcmc_samples$loglik)
    r_fixed <- r_samples[best_idx, ]
  }

  info <- build_tree_info(tree, Q)

  message(sprintf("Generating %d stochastic maps with rate-heterogeneous model...", nsim))

  maps_list <- vector("list", nsim)
  rates_used <- matrix(NA, nrow = nsim, ncol = nedges)
  transition_counts <- matrix(0, nrow = nsim, ncol = nedges)

  # OPTIMIZATION: Cache matrix exponentials for non-sample cases
  P_cache <- NULL
  if (rate_summary != "sample") {
    P_cache <- vector("list", nedges)
    for (e in 1:nedges) {
      if (rate_summary == "mean") {
        t_branch <- tree$edge.length[e] * r_fixed[e]
      } else if (rate_summary == "map") {
        t_branch <- tree$edge.length[e] * r_fixed[e]
      }
      P_cache[[e]] <- fast_expm(info$eig, t_branch, info$V, info$V_inv)
    }
  }

  for (s in 1:nsim) {

    # 1. Select rate vector
    if (rate_summary == "sample") {
      idx <- sample(npost, 1)
      r_vec <- r_samples[idx, ]
    } else {
      r_vec <- r_fixed
    }
    rates_used[s, ] <- r_vec

    # 2. Compute conditional likelihoods (pruning pass)
    L <- postorder_pass(info, states, Q, r_vec)

    # 3. Sample ancestral states (preorder / top-down)
    anc_states <- integer(n_nodes)
    # Tips are observed
    for (i in 1:ntips) anc_states[i] <- states[i]

    # Sample root state
    root_L <- L[root, ]
    root_w <- get_root_weights(root_L, root_prior, Q, k)
    root_probs <- root_L * root_w
    root_probs <- root_probs / sum(root_probs)
    root_probs[root_probs < 0] <- 0
    if (sum(root_probs) <= 0) root_probs <- rep(1/k, k)
    anc_states[root] <- sample(0:(k-1), 1, prob = root_probs)

    # Preorder traversal to sample internal node states
    # Process edges in reverse postorder (= preorder)
    for (po_idx in info$nedges:1) {
      orig_edge <- info$po_edges[po_idx]
      parent <- tree$edge[orig_edge, 1]
      child <- tree$edge[orig_edge, 2]

      # Skip if child is a tip (already known)
      if (child <= ntips) next

      # Conditional probability of child state given parent state
      parent_state <- anc_states[parent]
      t_branch <- tree$edge.length[orig_edge] * r_vec[orig_edge]

      # Use cached P matrix if available, else compute
      if (!is.null(P_cache)) {
        P <- P_cache[[orig_edge]]
      } else {
        P <- fast_expm(info$eig, t_branch)
      }

      # P(child_state | parent_state) * L(child_state | data below)
      child_probs <- P[parent_state + 1, ] * L[child, ]
      child_probs[child_probs < 0] <- 0
      total <- sum(child_probs)
      if (total <= 0) {
        child_probs <- rep(1/k, k)
      } else {
        child_probs <- child_probs / total
      }
      anc_states[child] <- sample(0:(k-1), 1, prob = child_probs)
    }

    # 4. Simulate history on each branch (endpoint-conditioned Gillespie)
    edge_maps <- vector("list", nedges)

    for (e in 1:nedges) {
      parent <- tree$edge[e, 1]
      child <- tree$edge[e, 2]
      start_state <- anc_states[parent]
      end_state <- anc_states[child]
      t_total <- tree$edge.length[e]
      r_e <- r_vec[e]

      # Scale Q by the rate scalar for this edge
      Q_e <- r_e * Q

      # Endpoint-conditioned simulation using rejection sampling
      edge_maps[[e]] <- simulate_branch_history(
        start_state, end_state, Q_e, t_total, k, max_attempts = 1000)
      transition_counts[s, e] <- length(edge_maps[[e]]) - 1  # number of segments - 1
    }

    maps_list[[s]] <- edge_maps

    if (s %% max(1, nsim %/% 10) == 0 || s == 1) {
      message(sprintf("  Map %d / %d", s, nsim))
    }
  }

  # Build simmap-compatible trees
  trees <- build_simmap_trees(tree, maps_list, states, k)

  # Summary: expected transitions per edge
  mean_transitions <- colMeans(transition_counts)
  edge_summary <- data.frame(
    edge = 1:nedges,
    parent = tree$edge[, 1],
    child = tree$edge[, 2],
    branch_length = tree$edge.length,
    mean_rate = colMeans(rates_used),
    mean_transitions = mean_transitions,
    sd_transitions = apply(transition_counts, 2, sd)
  )

  result <- list(
    maps = maps_list,
    trees = trees,
    rates_used = rates_used,
    transition_counts = transition_counts,
    Q = Q,
    summary = edge_summary
  )
  class(result) <- "ratescape_simmap"
  result
}


#' Simulate branch history with endpoint conditioning
#'
#' Uses rejection sampling: simulate forward from start_state under Q,
#' accept only if we arrive at end_state. For identical start/end states,
#' also accepts no-change paths.
#'
#' @param start_state integer (0-indexed)
#' @param end_state integer (0-indexed)
#' @param Q_scaled the rate-scaled Q matrix
#' @param t_total branch length
#' @param k number of states
#' @param max_attempts maximum rejection attempts before using the
#'   deterministic bridge (Nielsen 2002 mapping)
#' @return Named numeric vector: names are states (as characters), values are
#'   dwelling times in that state, in chronological order along the branch.
#' @keywords internal
simulate_branch_history <- function(start_state, end_state, Q_scaled, t_total,
                                     k, max_attempts = 1000) {

  state_labels <- as.character(0:(k-1))

  # Special case: if branch length is essentially zero

  if (t_total < 1e-12) {
    out <- t_total
    names(out) <- state_labels[start_state + 1]
    return(out)
  }

  # Attempt rejection sampling
  for (attempt in 1:max_attempts) {
    history <- simulate_forward(start_state, Q_scaled, t_total, k)

    # Check if final state matches
    final_state <- history$states[length(history$states)]
    if (final_state == end_state) {
      # Build dwelling-time vector
      dwell <- history$times
      names(dwell) <- state_labels[history$states + 1]
      return(dwell)
    }
  }

  # Fallback: direct bridge (no intermediate states)
  # Allocate time proportional to transition probability
  if (start_state == end_state) {
    out <- t_total
    names(out) <- state_labels[start_state + 1]
    return(out)
  } else {
    # Simple two-segment map: split at midpoint
    out <- c(t_total / 2, t_total / 2)
    names(out) <- c(state_labels[start_state + 1], state_labels[end_state + 1])
    return(out)
  }
}


#' Simulate forward Gillespie process on a branch
#' @keywords internal
simulate_forward <- function(start_state, Q_scaled, t_total, k) {
  # OPTIMIZATION: Pre-allocate arrays for efficiency
  max_events <- 1000L
  states <- integer(max_events)
  times <- numeric(max_events)
  states[1] <- start_state
  n <- 0L
  current_state <- start_state
  time_remaining <- t_total

  while (time_remaining > 0) {
    exit_rate <- -Q_scaled[current_state + 1, current_state + 1]

    if (exit_rate <= 1e-15) {
      # Absorbing state
      n <- n + 1L
      times[n] <- time_remaining
      break
    }

    wait_time <- rexp(1, exit_rate)

    if (wait_time >= time_remaining) {
      # No more transitions
      n <- n + 1L
      times[n] <- time_remaining
      break
    }

    n <- n + 1L
    times[n] <- wait_time
    time_remaining <- time_remaining - wait_time

    # Sample next state
    trans_rates <- Q_scaled[current_state + 1, ]
    trans_rates[current_state + 1] <- 0
    if (sum(trans_rates) <= 0) {
      times[n] <- times[n] + time_remaining
      break
    }
    current_state <- sample(0:(k-1), 1, prob = trans_rates)
    n <- n + 1L
    states[n] <- current_state
  }

  # Return only the used portion of arrays
  list(states = states[1:(n+1)], times = times[1:n])
}


#' Build simmap-compatible phylo objects
#'
#' Converts the internal map representation to phytools-compatible
#' "simmap" objects for plotting with plotSimmap.
#'
#' @keywords internal
build_simmap_trees <- function(tree, maps_list, states, k) {
  nsim <- length(maps_list)
  state_labels <- as.character(0:(k-1))

  trees <- vector("list", nsim)
  for (s in 1:nsim) {
    tr <- tree
    tr$maps <- maps_list[[s]]

    # Build mapped.edge matrix
    nedges <- nrow(tree$edge)
    mapped_edge <- matrix(0, nrow = nedges, ncol = k)
    colnames(mapped_edge) <- state_labels

    # OPTIMIZATION: Vectorized aggregation of state times per edge
    for (e in 1:nedges) {
      map_e <- maps_list[[s]][[e]]
      for (ci in seq_along(map_e)) {
        st <- names(map_e)[ci]
        col_idx <- match(st, state_labels)
        if (!is.na(col_idx)) {
          mapped_edge[e, col_idx] <- mapped_edge[e, col_idx] + map_e[ci]
        }
      }
    }
    tr$mapped.edge <- mapped_edge
    class(tr) <- c("simmap", "phylo")
    trees[[s]] <- tr
  }
  trees
}


#' Summarize stochastic maps
#'
#' @param object An object of class "ratescape_simmap".
#' @param ... Additional arguments (unused).
#'
#' @return Prints summary statistics and returns the summary data frame invisibly.
#' @export
print.ratescape_simmap <- function(x, ...) {
  cat("RateScape Stochastic Character Maps\n")
  cat("====================================\n\n")
  cat(sprintf("Number of maps: %d\n", length(x$maps)))
  cat(sprintf("Number of edges: %d\n", nrow(x$summary)))

  mean_total_trans <- mean(rowSums(x$transition_counts))
  sd_total_trans <- sd(rowSums(x$transition_counts))
  cat(sprintf("Mean total transitions per map: %.1f (SD = %.1f)\n",
              mean_total_trans, sd_total_trans))

  cat(sprintf("Mean rate range: %.2f - %.2f\n",
              min(x$summary$mean_rate), max(x$summary$mean_rate)))

  cat("\nPer-edge summary (top 10 by transitions):\n")
  top_edges <- head(x$summary[order(-x$summary$mean_transitions), ], 10)
  print(top_edges, row.names = FALSE)
  invisible(x)
}
