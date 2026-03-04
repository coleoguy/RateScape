#' Threshold Model for Discrete Characters
#'
#' Fits Felsenstein's threshold model, where discrete character states arise from
#' an underlying continuous "liability" that evolves by Brownian motion on the
#' phylogeny. The observed discrete states are determined by where the liability
#' falls relative to one or more thresholds. This provides an alternative to
#' Markov models when discrete states are believed to reflect an underlying
#' quantitative genetic trait.
#'
#' @param tree An object of class "phylo".
#' @param data Data frame with tip states (single column). States must be
#'   integers 0, 1, 2, ..., (k-1) where k is the number of states.
#' @param k Integer. Number of discrete states (inferred from data if NULL).
#' @param ngen Integer. MCMC generations. Default 50000.
#' @param burn_in Integer. Burn-in iterations. Default 10000.
#' @param thin Integer. Thinning interval. Default 20.
#' @param sigma2_prior Character or numeric. Prior on the BM rate sigma^2.
#'   "diffuse" (default) uses Exp(0.1), "moderate" uses Exp(1), or provide
#'   a numeric value for a custom Exp rate parameter.
#' @param root_prior Character. Prior on root liability: "diffuse" uses
#'   N(0, 100) (default), "moderate" uses N(0, 10).
#' @param seed Integer or NULL. Random seed.
#'
#' @details
#'
#' **The threshold model:**
#'
#' 1. An unobserved continuous liability L(t) evolves by Brownian motion with
#'    rate sigma^2 on the phylogeny. At the root, L ~ N(mu_root, sigma^2_root).
#'
#' 2. For k discrete states, there are k-1 thresholds:
#'    t_1 < t_2 < ... < t_{k-1}
#'
#' 3. The observed state at tip i is:
#'    - State 0 if L_i < t_1
#'    - State j if t_j <= L_i < t_{j+1}  (for j = 1, ..., k-2)
#'    - State k-1 if L_i >= t_{k-1}
#'
#' For identifiability, we fix t_1 = 0 and sigma^2 is estimated. The remaining
#' k-2 thresholds (for k > 2) and the root liability are free parameters.
#'
#' **MCMC inference:**
#'
#' The liabilities at all internal nodes and tips are latent variables. We use
#' data augmentation MCMC:
#'
#' 1. Sample liabilities at tips (constrained to match their observed state).
#' 2. Sample liabilities at internal nodes from their BM conditional.
#' 3. Update sigma^2 from its full conditional.
#' 4. Update thresholds (for k > 2) via MH proposals.
#' 5. Update root liability via MH.
#'
#' **Why this is useful:**
#'
#' The Mk model treats discrete states as fundamentally categorical. The
#' threshold model instead posits a latent continuous process. This is
#' appropriate when:
#'
#' - The discrete states represent categories along a continuum (e.g., body
#'   size classes, petal number, habitat categories)
#' - You want to estimate the phylogenetic signal of the underlying liability
#' - You need ancestral reconstructions in liability space
#' - Heritability / variance component estimation is the goal
#'
#' **Comparison with Mk:**
#'
#' The threshold model and Mk model make different assumptions. The threshold
#' model naturally produces ordered transitions (because the liability must
#' cross thresholds in sequence), while Mk allows any transition. Use
#' \code{compareModels} to formally compare the two approaches.
#'
#' @return An object of class "ratescapeThreshold" containing:
#'   \describe{
#'     \item{sigma2}{Posterior samples of BM rate.}
#'     \item{thresholds}{Matrix of posterior threshold samples (nsamples x (k-2)).}
#'     \item{root_liability}{Posterior samples of root liability.}
#'     \item{tip_liabilities}{Matrix (nsamples x ntips) of tip liability samples.}
#'     \item{node_liabilities}{Matrix (nsamples x n_internal) of internal node
#'       liability samples.}
#'     \item{ancestral_states}{Matrix (nsamples x n_internal) of predicted
#'       ancestral discrete states (from liability thresholding).}
#'     \item{loglik}{Posterior log-likelihood samples.}
#'     \item{tree}{The input phylo object.}
#'     \item{k}{Number of discrete states.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Binary threshold model
#'   fit <- fitThreshold(tree, data)
#'   summary(fit)
#'
#'   # Multi-state (ordered) threshold model
#'   fit3 <- fitThreshold(tree, data, k = 3)
#' }
#'
#' @export
fitThreshold <- function(
    tree, data, k = NULL,
    ngen = 50000, burn_in = 10000, thin = 20,
    sigma2_prior = "diffuse",
    root_prior = "diffuse",
    seed = NULL) {

  if (!inherits(tree, "phylo")) stop("tree must be 'phylo'")
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

  if (k < 2) stop("k must be at least 2")

  # Prior setup
  if (is.character(sigma2_prior)) {
    sigma2_prior <- match.arg(sigma2_prior, c("diffuse", "moderate"))
    lambda_s2 <- switch(sigma2_prior, diffuse = 0.1, moderate = 1.0)
  } else {
    lambda_s2 <- as.numeric(sigma2_prior)
  }
  root_prior <- match.arg(root_prior, c("diffuse", "moderate"))
  root_sd <- switch(root_prior, diffuse = 10, moderate = sqrt(10))

  # Tree structure
  n_nodes <- ntips + tree$Nnode
  root <- ntips + 1
  nedges <- nrow(tree$edge)

  # Parent-child structure
  parent_of <- integer(n_nodes)
  parent_of[root] <- 0
  for (e in 1:nedges) {
    parent_of[tree$edge[e, 2]] <- tree$edge[e, 1]
  }

  # Branch lengths indexed by child node
  bl <- numeric(n_nodes)
  for (e in 1:nedges) {
    bl[tree$edge[e, 2]] <- tree$edge.length[e]
  }

  # Children of each node
  children <- vector("list", n_nodes)
  for (i in 1:n_nodes) children[[i]] <- integer(0)
  for (e in 1:nedges) {
    children[[tree$edge[e, 1]]] <- c(children[[tree$edge[e, 1]]], tree$edge[e, 2])
  }

  # Thresholds: fix t_1 = 0 for identifiability
  # Free thresholds: t_2, t_3, ..., t_{k-1} (for k > 2)
  n_free_thresh <- max(0, k - 2)
  if (n_free_thresh > 0) {
    # Initialize thresholds at equal spacing
    thresh_free <- seq(1, n_free_thresh, length.out = n_free_thresh)
  } else {
    thresh_free <- numeric(0)
  }

  get_thresholds <- function(free_t) {
    c(0, free_t)  # t_1 = 0, then free thresholds
  }

  # Map liability to discrete state
  liability_to_state <- function(l, thresholds) {
    # State j if thresholds[j] <= l < thresholds[j+1]
    # State 0 if l < thresholds[1]; State k-1 if l >= thresholds[k-1]
    s <- sum(l >= thresholds)
    min(s, k - 1)
  }

  # Initialize liabilities
  sigma2_current <- 1.0
  liabilities <- numeric(n_nodes)

  # Root liability
  liabilities[root] <- 0.5

  # Initialize internal nodes by averaging children
  tree_po <- ape::reorder.phylo(tree, "postorder")
  for (i in 1:nedges) {
    parent <- tree_po$edge[i, 1]
    child <- tree_po$edge[i, 2]
    if (child <= ntips) {
      # Tips: initialize to center of their state's interval
      thresholds <- get_thresholds(thresh_free)
      lower <- if (states[child] == 0) -3 else thresholds[states[child]]
      upper <- if (states[child] == k-1) 3 + thresholds[length(thresholds)] else thresholds[states[child] + 1]
      liabilities[child] <- (lower + upper) / 2
    }
  }

  # Initialize internal nodes from tips upward
  for (i in 1:nedges) {
    parent <- tree_po$edge[i, 1]
    child <- tree_po$edge[i, 2]
    if (child > ntips) {
      child_vals <- liabilities[children[[child]]]
      liabilities[child] <- mean(child_vals)
    }
  }

  # MCMC setup
  niter <- ngen + burn_in
  nsamples <- floor(ngen / thin)

  sigma2_samples <- numeric(nsamples)
  thresh_samples <- if (n_free_thresh > 0) matrix(NA, nsamples, n_free_thresh) else NULL
  root_samples <- numeric(nsamples)
  tip_liab_samples <- matrix(NA, nsamples, ntips)
  node_liab_samples <- matrix(NA, nsamples, tree$Nnode)
  anc_state_samples <- matrix(NA, nsamples, tree$Nnode)

  # MH proposal scales
  tau_thresh <- 0.3
  tau_root <- 0.5

  message("----------------------------------------------------")
  message(sprintf("  RateScape Threshold Model MCMC"))
  message(sprintf("  Tree:     %d tips", ntips))
  message(sprintf("  States:   k = %d (%d free thresholds)", k, n_free_thresh))
  message(sprintf("  Chain:    %d gen (%d burn-in + %d sampling)", niter, burn_in, ngen))
  message("----------------------------------------------------")

  mcmc_start_time <- proc.time()[3]

  for (iter in 1:niter) {
    thresholds <- get_thresholds(thresh_free)

    # ---- 1. Sample tip liabilities (truncated normal) ----
    for (tip in 1:ntips) {
      parent <- parent_of[tip]
      parent_liab <- liabilities[parent]
      t_branch <- bl[tip]
      cond_var <- sigma2_current * t_branch

      # Liability must be in [lower, upper) for observed state
      s <- states[tip]
      lower <- if (s == 0) -Inf else thresholds[s]
      upper <- if (s == k - 1) Inf else thresholds[s + 1]

      # Draw from truncated normal: N(parent_liab, cond_var) in [lower, upper]
      liabilities[tip] <- rtruncnorm(parent_liab, sqrt(cond_var), lower, upper)
    }

    # ---- 2. Sample internal node liabilities ----
    # For each internal node (not root): conditional on parent and children
    # Preorder traversal (root -> tips)
    for (po_idx in nedges:1) {
      orig_e <- tree_po$edge[po_idx, 1]  # parent in postorder
      child <- tree_po$edge[po_idx, 2]

      if (child <= ntips) next  # skip tips
      if (child == root) next   # root handled separately

      parent <- parent_of[child]
      t_from_parent <- bl[child]

      # BM conditional: given parent value and children values
      child_nodes <- children[[child]]
      child_bls <- bl[child_nodes]

      # Precision from parent
      prec_parent <- 1 / (sigma2_current * t_from_parent)
      # Precision from each child
      prec_children <- 1 / (sigma2_current * child_bls)

      total_prec <- prec_parent + sum(prec_children)
      cond_mean <- (prec_parent * liabilities[parent] +
                     sum(prec_children * liabilities[child_nodes])) / total_prec
      cond_var <- 1 / total_prec

      liabilities[child] <- rnorm(1, cond_mean, sqrt(cond_var))
    }

    # ---- 3. Sample root liability ----
    # Conditional on children and prior
    child_nodes <- children[[root]]
    child_bls <- bl[child_nodes]
    prec_children <- 1 / (sigma2_current * child_bls)
    prec_prior <- 1 / root_sd^2

    total_prec <- prec_prior + sum(prec_children)
    cond_mean <- (prec_prior * 0 +
                   sum(prec_children * liabilities[child_nodes])) / total_prec
    cond_var <- 1 / total_prec
    liabilities[root] <- rnorm(1, cond_mean, sqrt(cond_var))

    # ---- 4. Gibbs update for sigma^2 ----
    # sigma^2 | liabilities ~ InverseGamma
    # Vectorized sum of squared contrasts
    parent_liabs <- liabilities[tree$edge[, 1]]
    child_liabs <- liabilities[tree$edge[, 2]]
    bls <- tree$edge.length
    valid <- bls > 0
    ss <- sum((child_liabs[valid] - parent_liabs[valid])^2 / bls[valid])

    # Posterior: InverseGamma(a_post, b_post)
    a_post <- nedges / 2 + 1  # shape (adding 1 for prior)
    b_post <- ss / 2 + lambda_s2  # rate (adding prior contribution)
    sigma2_current <- 1 / rgamma(1, shape = a_post, rate = b_post)
    sigma2_current <- max(sigma2_current, 1e-10)

    # ---- 5. MH for thresholds (k > 2) ----
    if (n_free_thresh > 0) {
      for (ti in 1:n_free_thresh) {
        thresh_prop <- thresh_free
        thresh_prop[ti] <- thresh_free[ti] + rnorm(1, 0, tau_thresh)

        # Maintain ordering: thresh_free must be increasing
        all_thresh <- c(0, thresh_prop)
        if (is.unsorted(all_thresh)) next  # reject if ordering violated

        # Vectorized consistency check: all tips must remain in correct state
        tip_liabs <- liabilities[1:ntips]
        all_bounds <- c(-Inf, all_thresh, Inf)
        s_pred <- findInterval(tip_liabs, all_bounds) - 1L
        if (all(s_pred == states)) {
          thresh_free[ti] <- thresh_prop[ti]
        }
      }
    }

    # ---- Store ----
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      si <- (iter - burn_in) %/% thin
      if (si >= 1 && si <= nsamples) {
        sigma2_samples[si] <- sigma2_current
        if (n_free_thresh > 0) thresh_samples[si, ] <- thresh_free
        root_samples[si] <- liabilities[root]
        tip_liab_samples[si, ] <- liabilities[1:ntips]
        node_liab_samples[si, ] <- liabilities[(ntips+1):n_nodes]

        # Ancestral states from liabilities
        curr_thresh <- get_thresholds(thresh_free)
        for (n_idx in 1:tree$Nnode) {
          node <- ntips + n_idx
          anc_state_samples[si, n_idx] <- liability_to_state(liabilities[node], curr_thresh)
        }
      }
    }

    # Progress
    if (iter %% 5000 == 0) {
      pct <- iter / niter
      phase <- if (iter <= burn_in) "burn-in" else "sampling"
      elapsed <- proc.time()[3] - mcmc_start_time
      message(sprintf("  %3.0f%% | %s | sigma2=%.3f root=%.2f | %.0fs",
                      pct*100, phase, sigma2_current, liabilities[root], elapsed))
    }
  }

  elapsed_total <- proc.time()[3] - mcmc_start_time
  message(sprintf("  Threshold MCMC complete in %.0fs", elapsed_total))

  result <- list(
    sigma2 = sigma2_samples,
    thresholds = thresh_samples,
    root_liability = root_samples,
    tip_liabilities = tip_liab_samples,
    node_liabilities = node_liab_samples,
    ancestral_states = anc_state_samples,
    tree = tree,
    k = k,
    states = states,
    nsamples = nsamples,
    call = match.call()
  )

  class(result) <- "ratescapeThreshold"
  result
}


#' Truncated normal sampler
#'
#' Samples from N(mu, sigma) truncated to [a, b] using inverse-CDF method.
#'
#' @keywords internal
rtruncnorm <- function(mu, sigma, a, b) {
  # Handle edge cases
  if (sigma <= 0) return(mu)
  if (a == -Inf && b == Inf) return(rnorm(1, mu, sigma))

  Fa <- pnorm(a, mu, sigma)
  Fb <- pnorm(b, mu, sigma)

  # Guard against both being 0 or 1
  if (Fb - Fa < 1e-15) {
    return((a + b) / 2)  # fallback: midpoint
  }

  u <- runif(1, Fa, Fb)
  qnorm(u, mu, sigma)
}


#' Print method for threshold model
#' @export
print.ratescapeThreshold <- function(x, ...) {
  cat("RateScape Threshold Model\n")
  cat("=========================\n\n")
  cat(sprintf("States: %d | Tips: %d\n", x$k, length(x$states)))
  cat(sprintf("Posterior samples: %d\n\n", x$nsamples))

  cat(sprintf("BM rate (sigma^2): mean = %.4f, 95%% CI = [%.4f, %.4f]\n",
              mean(x$sigma2), quantile(x$sigma2, 0.025), quantile(x$sigma2, 0.975)))

  cat(sprintf("Root liability: mean = %.3f, 95%% CI = [%.3f, %.3f]\n",
              mean(x$root_liability),
              quantile(x$root_liability, 0.025),
              quantile(x$root_liability, 0.975)))

  if (!is.null(x$thresholds)) {
    cat("\nThresholds (t_1 = 0 fixed):\n")
    for (j in 1:ncol(x$thresholds)) {
      cat(sprintf("  t_%d: mean = %.3f, 95%% CI = [%.3f, %.3f]\n",
                  j + 1, mean(x$thresholds[, j]),
                  quantile(x$thresholds[, j], 0.025),
                  quantile(x$thresholds[, j], 0.975)))
    }
  }

  # Ancestral state summary for root
  root_states <- x$ancestral_states[, 1]
  state_probs <- tabulate(root_states + 1, nbins = x$k) / x$nsamples
  cat("\nRoot state posterior probabilities:\n")
  for (s in 0:(x$k-1)) {
    cat(sprintf("  State %d: %.3f\n", s, state_probs[s + 1]))
  }

  invisible(x)
}


#' Plot threshold model results
#'
#' Creates a multi-panel plot showing the posterior distributions of the BM
#' rate, threshold positions, and liability at the root.
#'
#' @param x An object of class "ratescapeThreshold".
#' @param y Unused.
#' @param ... Additional arguments passed to hist().
#'
#' @export
plot.ratescapeThreshold <- function(x, y = NULL, ...) {
  n_panels <- 2 + max(0, x$k - 2)  # sigma2 + root + free thresholds
  n_col <- min(3, n_panels)
  n_row <- ceiling(n_panels / n_col)

  old_par <- par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))
  on.exit(par(old_par))

  hist(x$sigma2, main = expression(sigma^2~"(BM rate)"), xlab = expression(sigma^2),
       col = "lightblue", border = "white", ...)

  hist(x$root_liability, main = "Root liability", xlab = "Liability",
       col = "lightgreen", border = "white", ...)

  if (!is.null(x$thresholds)) {
    for (j in 1:ncol(x$thresholds)) {
      hist(x$thresholds[, j], main = bquote(t[.(j+1)]~"(threshold)"),
           xlab = "Value", col = "lightyellow", border = "white", ...)
    }
  }
}
