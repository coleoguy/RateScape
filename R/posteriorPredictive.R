#' Posterior Predictive Model Adequacy Checks
#'
#' Simulates datasets from the posterior distribution and compares summary
#' statistics to the observed data. This provides a formal assessment of
#' whether the fitted model adequately captures the patterns in the data.
#'
#' @param fit An object of class "ratescapeFit" (from fitRateScape).
#' @param observed_data Data frame with tip states matching the original data.
#' @param nsim Integer. Number of posterior predictive simulations. Default 500.
#' @param test_statistics Character vector. Which test statistics to compute.
#'   Options include:
#'   \describe{
#'     \item{"transitions"}{Total expected number of state changes on the tree
#'       (estimated by parsimony). Tests whether the model generates datasets
#'       with similar amounts of character change.}
#'     \item{"entropy"}{Shannon entropy of the tip-state frequency distribution.
#'       Tests whether the model produces realistic state frequencies.}
#'     \item{"consistency"}{Retention index (RI). Tests whether the model
#'       captures the phylogenetic structure of the data.}
#'     \item{"rate_dispersion"}{Variance of reconstructed rates across branches.
#'       Tests whether the model captures the amount of rate heterogeneity.}
#'     \item{"max_clade_freq"}{Maximum frequency of any single state within
#'       the largest clade. Tests for spatial clustering of states.}
#'   }
#'   Default is all of the above.
#' @param seed Random seed for reproducibility.
#'
#' @details
#'
#' **Posterior predictive checking workflow:**
#'
#' 1. For each simulation s = 1, ..., nsim:
#'    a. Draw a parameter vector (r, z, pi, sigma2) from the posterior samples.
#'    b. Simulate a new dataset on the observed tree under these parameters.
#'    c. Compute each test statistic T_s on the simulated dataset.
#'
#' 2. Compute each test statistic T_obs on the observed dataset.
#'
#' 3. For each statistic, compute the posterior predictive p-value:
#'    p = Pr(T_sim >= T_obs), the proportion of simulated values that exceed
#'    the observed value. Extreme p-values (< 0.05 or > 0.95) indicate
#'    model inadequacy for that aspect of the data.
#'
#' **Interpretation:**
#'
#' - p near 0.5: model generates data consistent with observations for this
#'   statistic.
#' - p < 0.05: observed statistic is LARGER than most simulated values,
#'   indicating the model underestimates this feature (e.g., too few transitions).
#' - p > 0.95: observed statistic is SMALLER than most simulated values,
#'   indicating the model overestimates this feature.
#'
#' @return An object of class "ratescape_ppc" containing:
#'   \describe{
#'     \item{observed}{Named vector of observed test statistics.}
#'     \item{simulated}{Matrix (nsim x n_statistics) of simulated test statistics.}
#'     \item{p_values}{Named vector of posterior predictive p-values.}
#'     \item{effect_sizes}{Named vector of effect sizes (obs - mean(sim)) / sd(sim).}
#'     \item{adequate}{Logical vector: TRUE if 0.05 < p < 0.95 for each statistic.}
#'     \item{nsim}{Number of simulations performed.}
#'   }
#'
#' @examples
#' \dontrun{
#'   fit <- fitRateScape(tree, data, Q, lambda_sigma = 1.0)
#'   ppc <- posteriorPredictive(fit, data, nsim = 200)
#'   print(ppc)
#'   plot(ppc)
#' }
#'
#' @export
posteriorPredictive <- function(
    fit,
    observed_data,
    nsim = 500,
    test_statistics = c("transitions", "entropy", "consistency",
                        "rate_dispersion", "max_clade_freq"),
    seed = NULL) {

  if (!inherits(fit, "ratescapeFit")) {
    stop("fit must be an object of class 'ratescapeFit'")
  }
  if (!is.null(seed)) set.seed(seed)

  tree <- fit$tree
  Q <- fit$Q_fixed
  k <- nrow(Q)
  nedges <- nrow(tree$edge)
  ntips <- length(tree$tip.label)
  npost <- nrow(fit$mcmc_samples$r)

  # Process observed data
  obs_states <- as.integer(as.numeric(observed_data[, 1]))
  if (!is.null(rownames(observed_data))) {
    mi <- match(tree$tip.label, rownames(observed_data))
    if (!any(is.na(mi))) obs_states <- obs_states[mi]
  }
  if (min(obs_states) > 0) obs_states <- obs_states - min(obs_states)

  # Compute observed test statistics
  obs_stats <- compute_test_statistics(obs_states, tree, Q, k, test_statistics)

  message(sprintf("Posterior predictive check: %d simulations x %d statistics",
                  nsim, length(test_statistics)))

  # Simulate from posterior
  sim_stats <- matrix(NA, nrow = nsim, ncol = length(test_statistics))
  colnames(sim_stats) <- test_statistics

  for (s in 1:nsim) {
    # Draw parameter vector from posterior
    idx <- sample(npost, 1)
    r_vec <- fit$mcmc_samples$r[idx, ]
    z_vec <- fit$mcmc_samples$z[idx, ]
    pi_s <- fit$mcmc_samples$pi[idx]
    sigma2_s <- fit$mcmc_samples$sigma2[idx]

    # Simulate data under this parameter vector
    sim_result <- simRateScape(
      tree = tree,
      Q = Q,
      r_scalars = r_vec,
      nrep = 1,
      return_rates = FALSE
    )
    sim_states <- as.integer(sim_result$data[1, ])

    # Compute statistics
    sim_stats[s, ] <- compute_test_statistics(sim_states, tree, Q, k, test_statistics)

    if (s %% max(1, nsim %/% 10) == 0 || s == 1) {
      message(sprintf("  Simulation %d / %d", s, nsim))
    }
  }

  # Compute p-values and effect sizes
  p_values <- numeric(length(test_statistics))
  effect_sizes <- numeric(length(test_statistics))
  names(p_values) <- names(effect_sizes) <- test_statistics

  for (j in seq_along(test_statistics)) {
    sim_vals <- sim_stats[, j]
    obs_val <- obs_stats[j]
    p_values[j] <- mean(sim_vals >= obs_val, na.rm = TRUE)
    sim_sd <- sd(sim_vals, na.rm = TRUE)
    effect_sizes[j] <- if (sim_sd > 0) (obs_val - mean(sim_vals, na.rm = TRUE)) / sim_sd else 0
  }

  adequate <- p_values > 0.05 & p_values < 0.95

  result <- list(
    observed = obs_stats,
    simulated = sim_stats,
    p_values = p_values,
    effect_sizes = effect_sizes,
    adequate = adequate,
    nsim = nsim,
    test_statistics = test_statistics
  )
  class(result) <- "ratescape_ppc"
  result
}


#' Compute test statistics for a character vector on a tree
#'
#' @param states Integer vector of tip states (0-indexed).
#' @param tree phylo object.
#' @param Q Rate matrix.
#' @param k Number of states.
#' @param which_stats Character vector of statistics to compute.
#' @return Named numeric vector of test statistics.
#' @keywords internal
compute_test_statistics <- function(states, tree, Q, k, which_stats) {
  out <- numeric(length(which_stats))
  names(out) <- which_stats

  for (stat in which_stats) {
    out[stat] <- switch(stat,
      "transitions" = count_parsimony_changes(states, tree, k),
      "entropy" = compute_state_entropy(states, k),
      "consistency" = compute_retention_index(states, tree, k),
      "rate_dispersion" = compute_rate_dispersion(states, tree, Q, k),
      "max_clade_freq" = compute_max_clade_freq(states, tree, k),
      NA
    )
  }
  out
}


#' Count parsimony changes (Fitch algorithm)
#' @keywords internal
count_parsimony_changes <- function(states, tree, k) {
  ntips <- length(tree$tip.label)
  n_nodes <- ntips + tree$Nnode
  nedges <- nrow(tree$edge)

  # Fitch sets
  fitch_sets <- vector("list", n_nodes)
  for (i in 1:ntips) {
    fitch_sets[[i]] <- states[i]
  }

  tree_po <- ape::reorder.phylo(tree, "postorder")
  n_changes <- 0

  for (i in 1:nedges) {
    parent <- tree_po$edge[i, 1]
    child <- tree_po$edge[i, 2]

    if (is.null(fitch_sets[[parent]])) {
      fitch_sets[[parent]] <- fitch_sets[[child]]
    } else {
      isect <- intersect(fitch_sets[[parent]], fitch_sets[[child]])
      if (length(isect) > 0) {
        fitch_sets[[parent]] <- isect
      } else {
        fitch_sets[[parent]] <- union(fitch_sets[[parent]], fitch_sets[[child]])
        n_changes <- n_changes + 1
      }
    }
  }
  n_changes
}


#' Compute Shannon entropy of state distribution
#' @keywords internal
compute_state_entropy <- function(states, k) {
  counts <- tabulate(states + 1, nbins = k)
  p <- counts / sum(counts)
  -sum(p[p > 0] * log(p[p > 0]))
}


#' Compute retention index
#' @keywords internal
compute_retention_index <- function(states, tree, k) {
  observed_changes <- count_parsimony_changes(states, tree, k)
  n_states_present <- length(unique(states))
  min_changes <- max(0, n_states_present - 1)

  # Max changes: number of edges is an upper bound
  ntips <- length(tree$tip.label)
  max_changes <- ntips - 1  # worst case

  if (max_changes <= min_changes) return(1.0)

  ri <- (max_changes - observed_changes) / (max_changes - min_changes)
  max(0, min(1, ri))
}


#' Estimate rate dispersion from data
#'
#' Uses a quick local estimate: for each internal edge, count state changes
#' in the subtending clade and normalize by clade size.
#' @keywords internal
compute_rate_dispersion <- function(states, tree, Q, k) {
  # Simple proxy: variance of per-edge "change indicator" (does parent != child
  # under parsimony reconstruction)
  ntips <- length(tree$tip.label)
  n_nodes <- ntips + tree$Nnode
  nedges <- nrow(tree$edge)

  # Quick parsimony reconstruction
  node_states <- integer(n_nodes)
  node_states[1:ntips] <- states

  tree_po <- ape::reorder.phylo(tree, "postorder")
  for (i in 1:nedges) {
    parent <- tree_po$edge[i, 1]
    child <- tree_po$edge[i, 2]
    if (node_states[parent] == 0 && parent > ntips) {
      node_states[parent] <- node_states[child]
    }
  }

  # OPTIMIZATION: Vectorized edge-level change indicators
  change_indicator <- as.integer(node_states[tree$edge[, 1]] != node_states[tree$edge[, 2]])

  var(change_indicator)
}


#' Compute maximum state frequency in largest clade
#' @keywords internal
compute_max_clade_freq <- function(states, tree, k) {
  ntips <- length(tree$tip.label)
  nedges <- nrow(tree$edge)

  # OPTIMIZATION: Precompute children list for each node to avoid O(n²) lookups
  children_of <- vector("list", ntips + tree$Nnode)
  for (e in 1:nedges) {
    p <- tree$edge[e, 1]
    children_of[[p]] <- c(children_of[[p]], tree$edge[e, 2])
  }

  # Find all internal nodes and their descendant tips
  max_freq <- 0
  for (e in 1:nedges) {
    child <- tree$edge[e, 2]
    if (child > ntips) {
      # Get tips descending from this node using cached children list
      desc_tips <- get_descendants_fast(child, children_of, ntips)
      if (length(desc_tips) >= 4) {  # only check clades of reasonable size
        desc_states <- states[desc_tips]
        freq_table <- tabulate(desc_states + 1, nbins = k)
        max_freq <- max(max_freq, max(freq_table) / length(desc_tips))
      }
    }
  }
  max_freq
}


#' Get all tip descendants of an internal node (fast version with cached children)
#' @keywords internal
get_descendants_fast <- function(node, children_of, ntips) {
  tips <- integer(0)
  stack <- node
  while (length(stack) > 0) {
    current <- stack[length(stack)]
    stack <- stack[-length(stack)]
    children_vec <- children_of[[current]]
    if (!is.null(children_vec)) {
      for (child in children_vec) {
        if (child <= ntips) {
          tips <- c(tips, child)
        } else {
          stack <- c(stack, child)
        }
      }
    }
  }
  tips
}


#' Print method for posterior predictive check
#' @export
print.ratescape_ppc <- function(x, ...) {
  cat("RateScape Posterior Predictive Check\n")
  cat("=====================================\n\n")
  cat(sprintf("Simulations: %d\n\n", x$nsim))

  cat("Statistic            Observed    Sim Mean    Sim SD    p-value    Adequate?\n")
  cat("--------------------------------------------------------------------------\n")

  for (stat in x$test_statistics) {
    sim_mean <- mean(x$simulated[, stat], na.rm = TRUE)
    sim_sd <- sd(x$simulated[, stat], na.rm = TRUE)
    adequate_str <- if (x$adequate[stat]) "  OK" else "  ***"
    cat(sprintf("%-20s %8.3f    %8.3f    %6.3f    %7.3f    %s\n",
                stat, x$observed[stat], sim_mean, sim_sd,
                x$p_values[stat], adequate_str))
  }

  n_inadequate <- sum(!x$adequate)
  if (n_inadequate == 0) {
    cat("\nAll statistics adequate (0.05 < p < 0.95).\n")
  } else {
    cat(sprintf("\n%d statistic(s) flagged (p < 0.05 or p > 0.95).\n", n_inadequate))
    cat("Consider model misspecification for flagged statistics.\n")
  }

  invisible(x)
}


#' Plot posterior predictive distributions
#'
#' Creates a multi-panel plot showing the distribution of each test statistic
#' from the posterior predictive simulations, with the observed value marked.
#'
#' @param x An object of class "ratescape_ppc".
#' @param y Unused.
#' @param ... Additional arguments passed to hist().
#'
#' @export
plot.ratescape_ppc <- function(x, y = NULL, ...) {
  n_stats <- length(x$test_statistics)
  n_col <- min(3, n_stats)
  n_row <- ceiling(n_stats / n_col)

  old_par <- par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))
  on.exit(par(old_par))

  for (stat in x$test_statistics) {
    sim_vals <- x$simulated[, stat]
    obs_val <- x$observed[stat]

    hist(sim_vals, main = stat, xlab = "Value",
         col = "lightblue", border = "white", ...)
    abline(v = obs_val, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = sprintf("p = %.3f", x$p_values[stat]),
           bty = "n", cex = 0.9)
  }
}
