#' Compute log-likelihood under branch-specific rate model
#'
#' Uses Felsenstein's pruning algorithm with branch-specific rate scalars
#' applied to a single Q matrix. Each branch i has effective rate matrix
#' r_i * Q, which is equivalent to using Q with branch length r_i * t_i.
#'
#' @param tree A phylo object (must be in pruningwise order)
#' @param tip_states Integer vector of tip states (1-indexed), named by tip labels
#' @param Q Instantaneous rate matrix (k x k)
#' @param rates Numeric vector of rate scalars, one per edge (length = Nedge(tree))
#' @param root_prior Character or numeric vector. "equal" for equal frequencies,
#'   "stationary" for stationary distribution of Q, or a numeric vector of length k.
#' @return Log-likelihood (numeric scalar)
#' @keywords internal
pruning_likelihood <- function(tree, tip_states, Q, rates, root_prior = "equal") {
  ntip <- ape::Ntip(tree)
  nnode <- ape::Nnode(tree)
  nstates <- nrow(Q)
  nedge <- ape::Nedge(tree)

  # Reorder tree for pruning
  tree <- ape::reorder.phylo(tree, "postorder")

  # Initialize partial likelihood arrays
  # Rows = states, Cols = nodes (tips + internal)
  partial_lk <- matrix(0, nrow = nstates, ncol = ntip + nnode)

  # Set tip likelihoods
  for (i in 1:ntip) {
    state <- tip_states[tree$tip.label[i]]
    partial_lk[state, i] <- 1
  }

  # Precompute transition probability matrices for each edge
  P_list <- vector("list", nedge)
  for (e in 1:nedge) {
    eff_t <- rates[e] * tree$edge.length[e]
    P_list[[e]] <- expm::expm(Q * eff_t)
  }

  # Traverse in postorder
  for (e in seq(1, nedge, by = 2)) {
    # In postorder, edges come in pairs sharing a parent
    parent1 <- tree$edge[e, 1]
    child1 <- tree$edge[e, 2]
    parent2 <- tree$edge[e + 1, 1]
    child2 <- tree$edge[e + 1, 2]

    # Both edges should share the same parent in postorder traversal
    # Compute conditional likelihoods
    lk_child1 <- P_list[[e]] %*% partial_lk[, child1]
    lk_child2 <- P_list[[e + 1]] %*% partial_lk[, child2]

    partial_lk[, parent1] <- partial_lk[, parent1] + log(pmax(lk_child1, .Machine$double.xmin))
    # Wait -- we need to be careful. The standard pruning algorithm multiplies.
    # Let me redo this properly.
  }

  # Actually, let me rewrite this more carefully using the standard approach
  # Reset
  partial_lk <- matrix(0, nrow = nstates, ncol = ntip + nnode)

  # Set tip likelihoods
  for (i in 1:ntip) {
    state <- tip_states[tree$tip.label[i]]
    partial_lk[state, i] <- 1
  }

  # Track which children have been processed for each internal node
  # Use a clean postorder traversal
  node_children_lk <- list()

  for (e in 1:nedge) {
    parent <- tree$edge[e, 1]
    child <- tree$edge[e, 2]

    # Conditional likelihood of data below child given parent state
    cond_lk <- as.vector(P_list[[e]] %*% partial_lk[, child])

    # Multiply into parent's partial likelihood
    if (is.null(node_children_lk[[as.character(parent)]])) {
      node_children_lk[[as.character(parent)]] <- cond_lk
    } else {
      node_children_lk[[as.character(parent)]] <-
        node_children_lk[[as.character(parent)]] * cond_lk
    }

    # Check if parent now has all children processed
    # by checking if this is the last edge with this parent
    remaining <- sum(tree$edge[, 1] == parent & (1:nedge) > e)
    if (remaining == 0) {
      partial_lk[, parent] <- node_children_lk[[as.character(parent)]]
    }
  }

  # Root likelihood
  root_node <- ntip + 1
  root_lk <- partial_lk[, root_node]

  # Apply root prior
  if (is.character(root_prior) && root_prior == "equal") {
    pi_root <- rep(1 / nstates, nstates)
  } else if (is.character(root_prior) && root_prior == "stationary") {
    pi_root <- stationary_dist(Q)
  } else {
    pi_root <- root_prior
  }

  log_lik <- log(sum(root_lk * pi_root))
  return(log_lik)
}


#' Compute stationary distribution of a rate matrix
#'
#' @param Q A valid rate matrix
#' @return Numeric vector of stationary frequencies
#' @keywords internal
stationary_dist <- function(Q) {
  nstates <- nrow(Q)
  # Solve pi * Q = 0 subject to sum(pi) = 1
  A <- t(Q)
  A[nstates, ] <- 1
  b <- rep(0, nstates)
  b[nstates] <- 1
  pi <- solve(A, b)
  pi <- pmax(pi, 0)
  pi <- pi / sum(pi)
  return(pi)
}


#' Optimized pruning likelihood using scaling to avoid underflow
#'
#' For large state spaces (e.g., chromosome number), the standard pruning
#' algorithm can underflow. This version uses per-node log-scaling factors.
#'
#' @inheritParams pruning_likelihood
#' @return Log-likelihood (numeric scalar)
#' @keywords internal
pruning_likelihood_scaled <- function(tree, tip_states, Q, rates,
                                      root_prior = "equal") {
  ntip <- ape::Ntip(tree)
  nnode <- ape::Nnode(tree)
  nstates <- nrow(Q)
  nedge <- ape::Nedge(tree)

  tree <- ape::reorder.phylo(tree, "postorder")

  # Partial likelihoods and scaling factors
  partial_lk <- matrix(0, nrow = nstates, ncol = ntip + nnode)
  log_scale <- numeric(ntip + nnode)

  # Tip likelihoods
  for (i in 1:ntip) {
    state <- tip_states[tree$tip.label[i]]
    partial_lk[state, i] <- 1
  }

  # Precompute P matrices
  P_list <- vector("list", nedge)
  for (e in 1:nedge) {
    eff_t <- rates[e] * tree$edge.length[e]
    P_list[[e]] <- expm::expm(Q * eff_t)
  }

  # Accumulate children contributions
  node_children_lk <- list()
  node_children_logscale <- numeric(ntip + nnode)

  for (e in 1:nedge) {
    parent <- tree$edge[e, 1]
    child <- tree$edge[e, 2]

    cond_lk <- as.vector(P_list[[e]] %*% partial_lk[, child])

    if (is.null(node_children_lk[[as.character(parent)]])) {
      node_children_lk[[as.character(parent)]] <- cond_lk
      node_children_logscale[parent] <- log_scale[child]
    } else {
      node_children_lk[[as.character(parent)]] <-
        node_children_lk[[as.character(parent)]] * cond_lk
      node_children_logscale[parent] <-
        node_children_logscale[parent] + log_scale[child]
    }

    remaining <- sum(tree$edge[, 1] == parent & (1:nedge) > e)
    if (remaining == 0) {
      raw_lk <- node_children_lk[[as.character(parent)]]
      max_lk <- max(raw_lk)
      if (max_lk > 0) {
        partial_lk[, parent] <- raw_lk / max_lk
        log_scale[parent] <- log(max_lk) + node_children_logscale[parent]
      } else {
        partial_lk[, parent] <- raw_lk
        log_scale[parent] <- node_children_logscale[parent]
      }
    }
  }

  # Root
  root_node <- ntip + 1
  root_lk <- partial_lk[, root_node]

  if (is.character(root_prior) && root_prior == "equal") {
    pi_root <- rep(1 / nstates, nstates)
  } else if (is.character(root_prior) && root_prior == "stationary") {
    pi_root <- stationary_dist(Q)
  } else {
    pi_root <- root_prior
  }

  log_lik <- log(sum(root_lk * pi_root)) + log_scale[root_node]
  return(log_lik)
}
