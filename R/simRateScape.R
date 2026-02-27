#' Simulate Discrete Character Data Under Branch-Specific Rates
#'
#' Simulates discrete character evolution along a phylogeny where each branch
#' has its own rate scalar applied to a background Q matrix.
#'
#' @param tree A phylo object
#' @param Q A k x k instantaneous rate matrix
#' @param rates Numeric vector of rate scalars, one per edge. Default: all 1.
#' @param root_state Integer (1-indexed) specifying root state. Default: random
#'   draw from stationary distribution of Q.
#' @param nsim Number of independent simulations. Default 1.
#' @return If nsim = 1, a named integer vector of tip states. If nsim > 1,
#'   a matrix with nsim columns.
#' @export
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- rtree(50)
#' Q <- matrix(0.5, 3, 3)
#' diag(Q) <- -1.0
#'
#' # Simulate under homogeneous rate
#' data_homo <- simRateScape(tree, Q)
#'
#' # Simulate with some branches 3x faster
#' rates <- rep(1, Nedge(tree))
#' rates[1:10] <- 3.0
#' data_het <- simRateScape(tree, Q, rates = rates)
#' }
simRateScape <- function(tree, Q, rates = NULL, root_state = NULL, nsim = 1) {

  if (!inherits(tree, "phylo")) stop("tree must be a phylo object")
  if (!is.matrix(Q)) stop("Q must be a matrix")

  nstates <- nrow(Q)
  ntip <- ape::Ntip(tree)
  nedge <- ape::Nedge(tree)

  if (is.null(rates)) rates <- rep(1.0, nedge)
  if (length(rates) != nedge) stop("rates must have length equal to Nedge(tree)")

  # Reorder for preorder traversal
  tree <- ape::reorder.phylo(tree, "cladewise")

  results <- matrix(NA_integer_, nrow = ntip, ncol = nsim)
  rownames(results) <- tree$tip.label

  for (s in 1:nsim) {
    # Node states
    node_states <- integer(ntip + ape::Nnode(tree))

    # Root state
    if (is.null(root_state)) {
      pi <- stationary_dist(Q)
      node_states[ntip + 1] <- sample(1:nstates, 1, prob = pi)
    } else {
      node_states[ntip + 1] <- root_state
    }

    # Traverse in preorder (cladewise)
    for (e in 1:nedge) {
      parent <- tree$edge[e, 1]
      child <- tree$edge[e, 2]
      t_eff <- rates[e] * tree$edge.length[e]

      # Transition probability matrix
      P <- expm::expm(Q * t_eff)

      # Draw child state from parent
      parent_state <- node_states[parent]
      node_states[child] <- sample(1:nstates, 1, prob = P[parent_state, ])
    }

    results[, s] <- node_states[1:ntip]
  }

  if (nsim == 1) {
    out <- results[, 1]
    names(out) <- tree$tip.label
    return(out)
  } else {
    return(results)
  }
}


#' Simulate Rate Scalars from Spike-and-Slab Prior
#'
#' Generates branch-specific rate scalars from the spike-and-slab prior
#' used in the RateScape model. Useful for simulation studies.
#'
#' @param nedge Number of edges
#' @param pi Probability of spike (background rate). Default 0.7.
#' @param sigma2 Variance of log-normal slab. Default 0.5.
#' @return Numeric vector of rate scalars
#' @export
simRates <- function(nedge, pi = 0.7, sigma2 = 0.5) {
  z <- stats::rbinom(nedge, 1, prob = 1 - pi)  # 1 = slab (shifted)
  rates <- rep(1.0, nedge)
  n_slab <- sum(z)
  if (n_slab > 0) {
    rates[z == 1] <- stats::rlnorm(n_slab, meanlog = 0, sdlog = sqrt(sigma2))
  }
  return(rates)
}
