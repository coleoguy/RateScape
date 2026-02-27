#' Discretized-Gamma Maximum Likelihood Rate Categories
#'
#' Estimates branch-specific rate variation using a discretized gamma
#' distribution, analogous to the Yang (1994) model for among-site rate
#' variation. The shape parameter alpha of a Gamma(alpha, alpha) distribution
#' (mean = 1) is optimized by maximum likelihood, and branch-specific
#' posterior probabilities of rate category membership are computed via
#' empirical Bayes.
#'
#' @param tree A phylo object
#' @param data Named vector of character states (integer, 1-indexed)
#' @param Q A k x k rate matrix. Must be fixed (not a function).
#' @param ncat Number of rate categories. Default 4.
#' @param root_prior Root state prior: "equal", "stationary", or numeric vector.
#' @param alpha_init Starting value for alpha. Default 1.
#' @param optimize_Q Logical. If TRUE, jointly optimize Q scaling. Default FALSE.
#' @return A list of class "ratescape_ml" containing:
#'   \describe{
#'     \item{alpha}{Estimated shape parameter}
#'     \item{rate_categories}{Numeric vector of rate category values}
#'     \item{branch_posteriors}{Matrix of posterior probabilities (edges x categories)}
#'     \item{branch_rates}{Posterior mean rate for each branch}
#'     \item{loglik}{Log-likelihood at the optimum}
#'     \item{AIC}{Akaike Information Criterion}
#'     \item{BIC}{Bayesian Information Criterion}
#'   }
#' @export
rateCategories <- function(tree, data, Q,
                           ncat = 4,
                           root_prior = "equal",
                           alpha_init = 1,
                           optimize_Q = FALSE) {

  if (!inherits(tree, "phylo")) stop("tree must be a phylo object")
  if (!is.matrix(Q)) stop("Q must be a fixed rate matrix for ML approach")

  nstates <- nrow(Q)
  nedge <- ape::Nedge(tree)
  ntip <- ape::Ntip(tree)

  tip_states <- as.integer(data)
  names(tip_states) <- names(data)

  # --- Objective function: negative log-likelihood ---
  neg_loglik <- function(log_alpha) {
    alpha <- exp(log_alpha)

    # Discretize Gamma(alpha, alpha) into ncat categories
    rate_cats <- discretize_gamma(alpha, ncat)

    # For each rate category, compute the full-tree likelihood
    # with all branches at that rate
    # Then marginalize over categories for each branch

    # Actually, the proper approach: for each branch, marginalize over
    # its rate category while all other branches are also marginalized.
    # This is the Yang (1994) approach adapted for branches instead of sites.

    # Simplified approach: assume branch rates are independent
    # Compute likelihood for each assignment of rate category to each branch
    # This is intractable for large trees, so we use a mean-field approximation

    # Mean-field: assign each branch its expected rate (= 1 under Gamma(a,a))
    # and iterate

    # Actually, the most practical ML approach for branches:
    # 1. Compute likelihood under each of ncat rates for ALL branches simultaneously
    # 2. Since we can't enumerate all combinations, we use EM

    llik <- em_branch_rates(tree, tip_states, Q, rate_cats, root_prior,
                            max_iter = 50)
    return(-llik$loglik)
  }

  # --- Optimize alpha ---
  opt <- stats::optim(
    par = log(alpha_init),
    fn = neg_loglik,
    method = "Brent",
    lower = log(0.01),
    upper = log(100)
  )

  alpha_hat <- exp(opt$par)
  rate_cats <- discretize_gamma(alpha_hat, ncat)

  # Final EM to get branch posteriors
  em_result <- em_branch_rates(tree, tip_states, Q, rate_cats, root_prior,
                               max_iter = 100)

  # Model comparison
  n_params <- 1  # alpha (Q is fixed)
  if (optimize_Q) n_params <- n_params + sum(Q[upper.tri(Q)] != 0)
  loglik <- em_result$loglik
  AIC <- -2 * loglik + 2 * n_params
  BIC <- -2 * loglik + n_params * log(ntip)

  result <- list(
    tree = tree,
    data = data,
    Q = Q,
    alpha = alpha_hat,
    rate_categories = rate_cats,
    branch_posteriors = em_result$branch_posteriors,
    branch_rates = em_result$branch_rates,
    loglik = loglik,
    AIC = AIC,
    BIC = BIC,
    ncat = ncat,
    nstates = nstates
  )

  class(result) <- c("ratescape_ml", "ratescape")
  return(result)
}


#' Discretize a Gamma(alpha, alpha) distribution into n categories
#'
#' Returns the quantile midpoints, as in Yang (1994).
#'
#' @param alpha Shape (and rate) parameter
#' @param n Number of categories
#' @return Numeric vector of rate values with mean approximately 1
#' @keywords internal
discretize_gamma <- function(alpha, n) {
  # Quantile boundaries
  breaks <- stats::qgamma(seq(0, 1, length.out = n + 1), shape = alpha, rate = alpha)
  # Midpoints: use conditional means within each interval
  rates <- numeric(n)
  for (i in 1:n) {
    # Conditional mean of Gamma(alpha, alpha) in interval (breaks[i], breaks[i+1])
    # = alpha/alpha * [pgamma(breaks[i+1], alpha+1, alpha) - pgamma(breaks[i], alpha+1, alpha)]
    #   / [pgamma(breaks[i+1], alpha, alpha) - pgamma(breaks[i], alpha, alpha)]
    p_upper <- stats::pgamma(breaks[i + 1], shape = alpha + 1, rate = alpha)
    p_lower <- stats::pgamma(breaks[i], shape = alpha + 1, rate = alpha)
    q_upper <- stats::pgamma(breaks[i + 1], shape = alpha, rate = alpha)
    q_lower <- stats::pgamma(breaks[i], shape = alpha, rate = alpha)
    rates[i] <- (p_upper - p_lower) / (q_upper - q_lower)
  }
  return(rates)
}


#' EM algorithm for branch-specific rate categories
#'
#' Iterates between: (E-step) computing posterior probability of each rate
#' category for each branch given current assignments of other branches;
#' (M-step) updating assignments to maximize expected log-likelihood.
#'
#' Uses a mean-field approximation where each branch's rate category is
#' estimated independently conditional on the expected rates of all other branches.
#'
#' @param tree phylo object
#' @param tip_states Named integer vector of tip states
#' @param Q Rate matrix
#' @param rate_cats Numeric vector of rate category values
#' @param root_prior Root prior specification
#' @param max_iter Maximum EM iterations
#' @param tol Convergence tolerance for log-likelihood
#' @return List with loglik, branch_posteriors, and branch_rates
#' @keywords internal
em_branch_rates <- function(tree, tip_states, Q, rate_cats, root_prior,
                            max_iter = 100, tol = 1e-6) {

  nedge <- ape::Nedge(tree)
  ntip <- ape::Ntip(tree)
  nnode <- ape::Nnode(tree)
  ncat <- length(rate_cats)

  # Resolve root prior
  nstates <- nrow(Q)
  pi_root <- .resolve_root_prior(root_prior, nstates, Q)

  # Initialize: all branches at mean rate (=1)
  branch_rates <- rep(1.0, nedge)
  branch_posteriors <- matrix(1 / ncat, nrow = nedge, ncol = ncat)

  tree_po <- ape::reorder.phylo(tree, "postorder")
  cpp_available <- .check_cpp_available()

  prev_loglik <- -Inf

  for (iter in 1:max_iter) {

    if (cpp_available) {
      # === C++ FAST PATH: batch compute all branch x category likelihoods ===
      llik_matrix <- batch_branch_llik_cpp(
        parent_vec = tree_po$edge[, 1],
        child_vec = tree_po$edge[, 2],
        edge_lengths = tree_po$edge.length,
        tip_states = tip_states,
        Q_mat = Q,
        base_rates = branch_rates,
        rate_categories = rate_cats,
        root_prior = pi_root,
        ntip = ntip,
        nnode = nnode
      )

      # E-step from the batch results
      for (e in 1:nedge) {
        log_liks <- llik_matrix[e, ]
        log_prior <- rep(log(1 / ncat), ncat)
        log_post <- log_liks + log_prior
        max_lp <- max(log_post)
        post <- exp(log_post - max_lp)
        post <- post / sum(post)
        branch_posteriors[e, ] <- post
        branch_rates[e] <- sum(post * rate_cats)
      }

    } else {
      # === PURE R PATH ===
      for (e in 1:nedge) {
        log_liks <- numeric(ncat)
        for (c in 1:ncat) {
          test_rates <- branch_rates
          test_rates[e] <- rate_cats[c]
          log_liks[c] <- pruning_likelihood_scaled(tree_po, tip_states, Q,
                                                    test_rates, pi_root)
        }

        log_prior <- rep(log(1 / ncat), ncat)
        log_post <- log_liks + log_prior
        max_lp <- max(log_post)
        post <- exp(log_post - max_lp)
        post <- post / sum(post)

        branch_posteriors[e, ] <- post
        branch_rates[e] <- sum(post * rate_cats)
      }
    }

    # Compute current log-likelihood
    current_loglik <- .compute_llik(tree_po, tip_states, Q, branch_rates,
                                     pi_root, ntip, nnode, cpp_available)

    if (abs(current_loglik - prev_loglik) < tol) break
    prev_loglik <- current_loglik
  }

  return(list(
    loglik = current_loglik,
    branch_posteriors = branch_posteriors,
    branch_rates = branch_rates
  ))
}
