#' Compute Likelihood of Discrete Character Data Under Rate Scalars
#'
#' Implements the Felsenstein pruning algorithm with efficient partial-likelihood
#' caching for branch-specific rate scalar estimation.
#'
#' @param tree An object of class "phylo".
#' @param data Numeric vector of character states at the tips (0-indexed).
#' @param Q Transition rate matrix (k x k).
#' @param r_scalars Numeric vector of rate multipliers (one per edge).
#' @param root_prior Character: "fitzjohn", "equal", or "stationary".
#' @param log_scale Logical. If TRUE (default), return log-likelihood.
#'
#' @return Log-likelihood scalar.
#' @export
compute_likelihood <- function(tree, data, Q, r_scalars,
                                root_prior = "fitzjohn", log_scale = TRUE) {
  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
  if (!is.matrix(Q) || nrow(Q) != ncol(Q)) stop("Q must be square")

  k <- nrow(Q)
  n_tips <- length(tree$tip.label)
  if (length(data) != n_tips) stop("data length must match tips")

  data <- as.integer(data)
  if (length(r_scalars) == 1) r_scalars <- rep(r_scalars, nrow(tree$edge))
  if (length(r_scalars) != nrow(tree$edge)) stop("r_scalars length mismatch")
  if (any(r_scalars <= 0)) stop("r_scalars must be positive")

  root_prior <- match.arg(root_prior, c("fitzjohn", "equal", "stationary"))

  info <- build_tree_info(tree, Q)
  L <- postorder_pass(info, data, Q, r_scalars)

  root_L <- L[info$root, ]
  root_weights <- get_root_weights(root_L, root_prior, Q, k)

  total_lik <- sum(root_L * root_weights)
  if (total_lik <= 0) return(-1e20)
  return(log(total_lik))
}


#' Build tree info structure for efficient traversal
#'
#' All O(n) or O(n^2) preprocessing happens here once. Subsequent
#' likelihood computations reuse this cache.
#'
#' @keywords internal
build_tree_info <- function(tree, Q) {
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  nedges <- nrow(tree$edge)
  root <- n_tips + 1

  # Postorder traversal
 tree_po <- ape::reorder.phylo(tree, order = "postorder")

  # Vectorized edge matching: O(n) via paste-match instead of O(n^2) loop
  po_keys <- paste(tree_po$edge[, 1], tree_po$edge[, 2], sep = "_")
  orig_keys <- paste(tree$edge[, 1], tree$edge[, 2], sep = "_")
  po_to_orig <- match(po_keys, orig_keys)

  # Children list: which original edge indices lead FROM each node
  children_of <- vector("list", n_nodes)
  for (i in 1:n_nodes) children_of[[i]] <- integer(0)
  # Vectorized: group edges by parent
  parents <- tree$edge[, 1]
  for (e in seq_len(nedges)) {
    p <- parents[e]
    children_of[[p]] <- c(children_of[[p]], e)
  }

  # Parent edge for each node (vectorized assignment)
  parent_edge <- rep(NA_integer_, n_nodes)
  parent_edge[tree$edge[, 2]] <- seq_len(nedges)

  # Ancestor paths: pre-allocate with known max depth
  # Use integer vectors, avoid growing with c()
  ancestors <- vector("list", nedges)
  for (e in seq_len(nedges)) {
    node <- parents[e]
    path <- integer(0)
    while (!is.na(parent_edge[node])) {
      path[length(path) + 1L] <- parent_edge[node]
      node <- tree$edge[parent_edge[node], 1]
    }
    ancestors[[e]] <- path
  }

  # Eigendecomposition of Q + cache the inverse (avoids solve() on every expm call)
  eig <- eigen(Q)
  V <- eig$vectors
  V_inv <- solve(V)

  # Precompute stationary distribution for root_prior = "stationary"
  eig_t <- eigen(t(Q))
  idx <- which.min(abs(eig_t$values))
  stat_dist <- abs(Re(eig_t$vectors[, idx]))
  stat_dist <- stat_dist / sum(stat_dist)

  list(
    tree = tree,
    n_tips = n_tips,
    n_nodes = n_nodes,
    nedges = nedges,
    root = root,
    po_edges = po_to_orig,
    children_of = children_of,
    parent_edge = parent_edge,
    ancestors = ancestors,
    eig = eig,
    V = V,
    V_inv = V_inv,
    stat_dist = stat_dist
  )
}


#' Compute matrix exponential using cached eigendecomposition
#'
#' Uses pre-cached V and V_inv from build_tree_info to avoid
#' calling solve() on every branch.
#'
#' @keywords internal
fast_expm <- function(eig, t, V = NULL, V_inv = NULL) {
  k <- length(eig$values)
  if (t <= 0) return(diag(k))

  if (is.null(V)) V <- eig$vectors
  if (is.null(V_inv)) V_inv <- solve(V)

  # exp(D * t) as column-scaling: V %*% diag(exp_D) %*% V_inv
  # Efficient: scale columns of V by exp_D, then multiply by V_inv
  exp_D <- exp(eig$values * t)
  P <- Re((V * rep(exp_D, each = k)) %*% V_inv)
  P[P < 0] <- 0
  # Row-normalize to ensure valid transition probabilities
  rs <- rowSums(P)
  P <- P / rs
  P
}


#' Full postorder pass returning conditional likelihoods at all nodes
#' @keywords internal
postorder_pass <- function(info, data, Q, r_scalars) {
  k <- nrow(Q)
  n_nodes <- info$n_nodes
  n_tips <- info$n_tips
  nedges <- info$nedges

  # Pre-allocate L matrix
  L <- matrix(1.0, nrow = n_nodes, ncol = k)

  # Tips: set observed state columns (vectorized)
  tip_idx <- cbind(seq_len(n_tips), data + 1L)
  L[seq_len(n_tips), ] <- 0.0
  L[tip_idx] <- 1.0

  # Cache V and V_inv for reuse across all branches
  V <- info$V
  V_inv <- info$V_inv
  eig <- info$eig
  edge_mat <- info$tree$edge
  edge_lengths <- info$tree$edge.length

  # Postorder traversal
  po_edges <- info$po_edges
  for (po_idx in seq_len(nedges)) {
    orig_edge <- po_edges[po_idx]
    parent <- edge_mat[orig_edge, 1]
    child <- edge_mat[orig_edge, 2]
    t_branch <- edge_lengths[orig_edge] * r_scalars[orig_edge]

    P <- fast_expm(eig, t_branch, V, V_inv)
    # P %*% L[child, ] as matrix-vector product
    L[parent, ] <- L[parent, ] * (P %*% L[child, ])
  }

  L
}


#' Get root prior weights
#' @keywords internal
get_root_weights <- function(root_L, root_prior, Q, k, stat_dist = NULL) {
  if (root_prior == "fitzjohn") {
    s <- sum(root_L)
    if (s <= 0) return(rep(1 / k, k))
    return(root_L / s)
  } else if (root_prior == "equal") {
    return(rep(1 / k, k))
  } else {
    # Stationary: use precomputed if available
    if (!is.null(stat_dist)) return(stat_dist)
    eig_t <- eigen(t(Q))
    idx <- which.min(abs(eig_t$values))
    s <- abs(Re(eig_t$vectors[, idx]))
    return(s / sum(s))
  }
}


#' Compute likelihood with one edge changed (efficient partial update)
#'
#' Only recomputes the affected parent node and its ancestors to the root.
#' Returns both the log-likelihood and the updated L matrix.
#'
#' @keywords internal
compute_likelihood_one_edge <- function(info, data, Q, r_scalars, L_cache,
                                         edge_idx, new_r, root_prior,
                                         return_L = FALSE) {
  k <- nrow(Q)
  L <- L_cache  # R copies on modify (copy-on-write)
  V <- info$V
  V_inv <- info$V_inv
  eig <- info$eig
  edge_mat <- info$tree$edge
  edge_lengths <- info$tree$edge.length

  parent <- edge_mat[edge_idx, 1]

  # Recompute parent from its children
  L[parent, ] <- 1.0
  for (ce in info$children_of[[parent]]) {
    cn <- edge_mat[ce, 2]
    r_val <- if (ce == edge_idx) new_r else r_scalars[ce]
    tb <- edge_lengths[ce] * r_val
    P <- fast_expm(eig, tb, V, V_inv)
    L[parent, ] <- L[parent, ] * (P %*% L[cn, ])
  }

  # Propagate up to root through ancestor edges
  for (ae in info$ancestors[[edge_idx]]) {
    ap <- edge_mat[ae, 1]
    L[ap, ] <- 1.0
    for (ce in info$children_of[[ap]]) {
      cn <- edge_mat[ce, 2]
      tb <- edge_lengths[ce] * r_scalars[ce]
      P <- fast_expm(eig, tb, V, V_inv)
      L[ap, ] <- L[ap, ] * (P %*% L[cn, ])
    }
  }

  root_L <- L[info$root, ]
  rw <- get_root_weights(root_L, root_prior, Q, k, info$stat_dist)
  tl <- sum(root_L * rw)
  ll <- if (tl <= 0) -1e20 else log(tl)

  if (return_L) return(list(loglik = ll, L = L))
  ll
}


#' Fallback matrix exponential (standalone, no cached decomposition)
#' @keywords internal
matrix_expm <- function(Q, t) {
  if (t <= 0) return(diag(nrow(Q)))
  eig <- eigen(Q)
  V <- eig$vectors
  V_inv <- solve(V)
  k <- nrow(Q)
  exp_D <- exp(eig$values * t)
  P <- Re((V * rep(exp_D, each = k)) %*% V_inv)
  P[P < 0] <- 0
  P <- P / rowSums(P)
  P
}
