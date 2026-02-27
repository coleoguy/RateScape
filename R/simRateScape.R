#' Simulate Discrete Character Data with Branch Rate Heterogeneity
#'
#' Generates discrete character datasets under the RateScape spike-and-slab
#' or discretized gamma model. Useful for prior-predictive simulation, validation,
#' and simulation studies.
#'
#' @param tree An object of class "phylo".
#' @param Q Transition rate matrix (k × k).
#' @param lambda_sigma Numeric. Rate parameter for exponential prior on σ².
#'   Ignored if rates are provided directly via r_scalars.
#' @param pi Numeric in [0, 1]. Mixing proportion: probability that a branch
#'   is at background rate (r = 1). Default is 0.5.
#' @param nrep Integer. Number of independent replicates to simulate. Default is 1.
#' @param r_scalars Numeric vector. Pre-specified branch rate scalars. If provided,
#'   overrides the spike-and-slab prior. Default is NULL (sample from prior).
#' @param root_state Integer or NULL. Root character state (0 to k-1).
#'   If NULL, sampled from the root prior.
#' @param return_rates Logical. If TRUE, also return the branch-specific rates.
#'   Default is FALSE.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @details
#'
#' **Simulation procedure:**
#'
#' 1. If r_scalars is not provided, sample from the spike-and-slab prior:
#'    - z_i ~ Bernoulli(π) (spike indicator for each branch)
#'    - r_i = 1 if z_i = 1 (background rate)
#'    - r_i ~ LogNormal(0, σ²) if z_i = 0 (heterogeneous rate), where σ² ~ Exp(λ_σ)
#'
#' 2. If a root_state is not provided, sample from the root prior.
#'
#' 3. Evolve the character down the tree using the Gillespie algorithm:
#'    - At each branch with state s, compute the exit rate from s under
#'      the scaled rate matrix r_i * Q.
#'    - Time to transition is exponential with this rate.
#'    - Destination state is sampled from the scaled off-diagonals of Q.
#'
#' 4. Repeat for nrep independent simulations.
#'
#' @return A list containing:
#'   - `data`: A matrix of simulated character data (nrep rows × n_tips columns),
#'     with rows corresponding to replicates and columns to tips.
#'   - `rates`: (if return_rates = TRUE) List of branch rate scalars for each replicate.
#'   - `root_states`: Root states used for each replicate.
#'   - `tree`: The input tree.
#'   - `call`: The function call.
#'
#' @examples
#' \dontrun{
#'   tree <- ape::rtree(50, rooted = TRUE)
#'   Q <- makeQ(model = "mk", k = 2)
#'
#'   # Simulate with default spike-and-slab prior
#'   sim1 <- simRateScape(
#'     tree = tree,
#'     Q = Q,
#'     lambda_sigma = 1.0,
#'     pi = 0.7,
#'     nrep = 10
#'   )
#'   head(sim1$data)
#'
#'   # Simulate with fixed rates
#'   r_fixed <- rep(1, nrow(tree$edge))
#'   r_fixed[1:5] <- 2.0  # Five branches evolving 2× faster
#'   sim2 <- simRateScape(
#'     tree = tree,
#'     Q = Q,
#'     r_scalars = r_fixed,
#'     nrep = 10,
#'     return_rates = TRUE
#'   )
#' }
#'
#' @export
simRateScape <- function(
    tree,
    Q,
    lambda_sigma = NULL,
    pi = 0.5,
    nrep = 1,
    r_scalars = NULL,
    root_state = NULL,
    return_rates = FALSE,
    seed = NULL) {

  if (!inherits(tree, "phylo")) {
    stop("tree must be an object of class 'phylo'")
  }

  if (!is.matrix(Q) || nrow(Q) != ncol(Q)) {
    stop("Q must be a square matrix")
  }

  k <- nrow(Q)
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  nedges <- nrow(tree$edge)

  if (!is.null(r_scalars)) {
    if (length(r_scalars) != nedges) {
      stop("r_scalars must have length equal to number of edges")
    }
    if (any(r_scalars <= 0)) {
      stop("All r_scalars must be positive")
    }
  } else {
    if (is.null(lambda_sigma)) {
      stop("Either r_scalars or lambda_sigma must be provided")
    }
  }

  if (pi < 0 || pi > 1) {
    stop("pi must be in [0, 1]")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Initialize output
  data_matrix <- matrix(NA, nrow = nrep, ncol = n_tips)
  rates_list <- if (return_rates) vector("list", nrep) else NULL
  root_states <- numeric(nrep)

  message(sprintf(
    "Simulating %d replicate(s) on tree with %d tips and %d states",
    nrep, n_tips, k
  ))

  # Main simulation loop
  for (rep in 1:nrep) {

    # Sample or use provided rates
    if (is.null(r_scalars)) {
      z <- rbinom(nedges, 1, pi)  # Spike indicators
      sigma2 <- rexp(1, lambda_sigma)
      r <- numeric(nedges)
      for (i in 1:nedges) {
        r[i] <- if (z[i] == 1) 1.0 else rlnorm(1, 0, sqrt(sigma2))
      }
    } else {
      r <- r_scalars
    }

    if (return_rates) {
      rates_list[[rep]] <- r
    }

    # Sample root state
    if (is.null(root_state)) {
      # FitzJohn prior: uniform for simplicity (full version uses likelihood weighting)
      root_states[rep] <- sample(0:(k - 1), 1)
    } else {
      root_states[rep] <- root_state
    }

    # Evolve character down the tree
    states_at_nodes <- integer(n_nodes)
    states_at_nodes[n_tips + 1] <- root_states[rep]  # Root node state

    # Post-order traversal: evolve states down the tree
    for (edge_idx in 1:nedges) {
      parent_node <- tree$edge[edge_idx, 1]
      child_node <- tree$edge[edge_idx, 2]
      edge_len <- tree$edge.length[edge_idx]

      parent_state <- states_at_nodes[parent_node]

      # Evolve under scaled Q
      scaled_Q <- r[edge_idx] * Q
      child_state <- evolve_state(parent_state, scaled_Q, edge_len, k)

      states_at_nodes[child_node] <- child_state
    }

    # Extract tip states
    data_matrix[rep, ] <- states_at_nodes[1:n_tips]

    if (rep %% max(1, nrep / 10) == 0 || rep == 1) {
      message(sprintf("  Replicate %d / %d", rep, nrep))
    }
  }

  # Return result
  result <- list(
    data = data_matrix,
    root_states = root_states,
    tree = tree,
    call = match.call()
  )

  if (return_rates) {
    result$rates <- rates_list
  }

  class(result) <- c("ratescape_sim", "list")
  return(result)
}


#' Evolve character state down a single branch
#'
#' Implements Gillespie algorithm for character evolution under a scaled rate matrix.
#'
#' @keywords internal
evolve_state <- function(current_state, scaled_Q, branch_length, k) {

  time_remaining <- branch_length
  state <- current_state

  while (time_remaining > 0) {

    # Exit rate from current state
    exit_rate <- -scaled_Q[state + 1, state + 1]

    if (exit_rate <= 0) {
      # Absorbing state; remain in this state
      break
    }

    # Time to next transition (exponential)
    time_to_transition <- rexp(1, exit_rate)

    if (time_to_transition > time_remaining) {
      # No transition before branch end
      break
    }

    # Time has passed; now transition
    time_remaining <- time_remaining - time_to_transition

    # Sample next state from transition probabilities
    transition_rates <- scaled_Q[state + 1, -( state + 1)]
    if (sum(transition_rates) <= 0) {
      break
    }

    state <- sample(setdiff(0:(k - 1), state), 1, prob = transition_rates)
  }

  return(state)
}
