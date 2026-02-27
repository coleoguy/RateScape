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

  # Build tree info cache if not present
  info <- build_tree_info(tree, Q)

  # Full postorder pass
  L <- postorder_pass(info, data, Q, r_scalars)

  # Root likelihood
  root_L <- L[info$root, ]
  root_weights <- get_root_weights(root_L, root_prior, Q, k)

  total_lik <- sum(root_L * root_weights)
  if (total_lik <= 0) return(-1e20)
  return(log(total_lik))
}


#' Build tree info structure for efficient traversal
#' @keywords internal
build_tree_info <- function(tree, Q) {
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  nedges <- nrow(tree$edge)
  root <- n_tips + 1

  # Get postorder edge sequence
  tree_po <- ape::reorder.phylo(tree, order = "postorder")

  # Map postorder edges back to original edge indices
  po_to_orig <- integer(nedges)
  for (i in 1:nedges) {
    p <- tree_po$edge[i, 1]
    c <- tree_po$edge[i, 2]
    orig <- which(tree$edge[, 1] == p & tree$edge[, 2] == c)
    po_to_orig[i] <- orig
  }

  # Children list (in original edge indexing)
  children_of <- vector("list", n_nodes)
  for (i in 1:n_nodes) children_of[[i]] <- integer(0)
  for (e in 1:nedges) {
    p <- tree$edge[e, 1]
    children_of[[p]] <- c(children_of[[p]], e)
  }

  # Parent edge for each node
  parent_edge <- integer(n_nodes)
  parent_edge[] <- NA
  for (e in 1:nedges) {
    child <- tree$edge[e, 2]
    parent_edge[child] <- e
  }

  # Path from each edge to root (for partial updates)
  # For each edge e (parent -> child), the ancestors are the edges
  # from parent up to root
  ancestors <- vector("list", nedges)
  for (e in 1:nedges) {
    node <- tree$edge[e, 1]  # parent node of this edge
    path <- integer(0)
    while (!is.na(parent_edge[node])) {
      path <- c(path, parent_edge[node])
      node <- tree$edge[parent_edge[node], 1]
    }
    ancestors[[e]] <- path
  }

  # Eigendecomposition of Q for fast matrix expm
  eig <- eigen(Q)

  list(
    tree = tree,
    n_tips = n_tips,
    n_nodes = n_nodes,
    nedges = nedges,
    root = root,
    po_edges = po_to_orig,  # postorder edge indices (in original numbering)
    po_parents = tree_po$edge[, 1],
    po_children = tree_po$edge[, 2],
    children_of = children_of,
    parent_edge = parent_edge,
    ancestors = ancestors,
    eig = eig
  )
}


#' Compute matrix exponential using cached eigendecomposition
#' @keywords internal
fast_expm <- function(eig, t) {
  if (t <= 0) return(diag(length(eig$values)))
  V <- eig$vectors
  exp_D <- exp(eig$values * t)
  P <- Re(V %*% (exp_D * solve(V)))  # equivalent to V %*% diag(exp_D) %*% V^{-1} but faster
  P[P < 0] <- 0
  P <- P / rowSums(P)
  P
}


#' Full postorder pass returning conditional likelihoods at all nodes
#' @keywords internal
postorder_pass <- function(info, data, Q, r_scalars) {
  k <- nrow(Q)
  L <- matrix(0, nrow = info$n_nodes, ncol = k)

  # Tips
  for (i in 1:info$n_tips) {
    L[i, data[i] + 1] <- 1.0
  }

  # Internal nodes: init to 1 (multiplicative identity)
  for (i in (info$n_tips + 1):info$n_nodes) {
    L[i, ] <- 1.0
  }

  # Postorder traversal
  for (po_idx in 1:info$nedges) {
    orig_edge <- info$po_edges[po_idx]
    parent <- info$tree$edge[orig_edge, 1]
    child <- info$tree$edge[orig_edge, 2]
    t_branch <- info$tree$edge.length[orig_edge] * r_scalars[orig_edge]

    P <- fast_expm(info$eig, t_branch)
    child_contrib <- as.vector(P %*% L[child, ])
    L[parent, ] <- L[parent, ] * child_contrib
  }

  L
}


#' Get root prior weights
#' @keywords internal
get_root_weights <- function(root_L, root_prior, Q, k) {
  if (root_prior == "fitzjohn") {
    w <- root_L / sum(root_L)
    if (any(is.nan(w))) w <- rep(1/k, k)
    return(w)
  } else if (root_prior == "equal") {
    return(rep(1/k, k))
  } else {
    eig <- eigen(t(Q))
    idx <- which.min(abs(eig$values))
    stat <- abs(Re(eig$vectors[, idx]))
    return(stat / sum(stat))
  }
}


#' Compute likelihood with one edge changed (efficient partial update)
#'
#' Returns both the log-likelihood and the updated L matrix.
#'
#' @keywords internal
compute_likelihood_one_edge <- function(info, data, Q, r_scalars, L_cache,
                                         edge_idx, new_r, root_prior,
                                         return_L = FALSE) {
  k <- nrow(Q)
  L <- L_cache  # R copies on modify

  parent <- info$tree$edge[edge_idx, 1]

  # Recompute parent node from scratch
  L[parent, ] <- 1.0
  for (ce in info$children_of[[parent]]) {
    cn <- info$tree$edge[ce, 2]
    tb <- info$tree$edge.length[ce] * (if (ce == edge_idx) new_r else r_scalars[ce])
    P <- fast_expm(info$eig, tb)
    L[parent, ] <- L[parent, ] * as.vector(P %*% L[cn, ])
  }

  # Recompute ancestors
  for (ae in info$ancestors[[edge_idx]]) {
    ap <- info$tree$edge[ae, 1]
    L[ap, ] <- 1.0
    for (ce in info$children_of[[ap]]) {
      cn <- info$tree$edge[ce, 2]
      tb <- info$tree$edge.length[ce] * r_scalars[ce]
      P <- fast_expm(info$eig, tb)
      L[ap, ] <- L[ap, ] * as.vector(P %*% L[cn, ])
    }
  }

  root_L <- L[info$root, ]
  rw <- get_root_weights(root_L, root_prior, Q, k)
  tl <- sum(root_L * rw)
  ll <- if (tl <= 0) -1e20 else log(tl)

  if (return_L) return(list(loglik = ll, L = L))
  ll
}


#' Fallback matrix_expm using eigendecomposition
#' @keywords internal
matrix_expm <- function(Q, t) {
  if (t <= 0) return(diag(nrow(Q)))
  eig <- eigen(Q)
  V <- eig$vectors
  exp_D <- exp(eig$values * t)
  P <- Re(V %*% diag(exp_D) %*% solve(V))
  P[P < 0] <- 0
  P <- P / rowSums(P)
  P
}
