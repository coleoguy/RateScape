#' Maximum Likelihood Discretized Gamma Rate Categories
#'
#' Fits a discretized gamma model for rate heterogeneity using maximum likelihood
#' and the EM algorithm. Automatically selects the optimal number of rate categories
#' via BIC model selection.
#'
#' @param tree An object of class "phylo" (ape package).
#' @param data A data frame with rows as species and a single column of discrete states.
#' @param Q A transition rate matrix (k x k).
#' @param k_min Integer. Minimum number of rate categories. Default is 2.
#' @param k_max Integer. Maximum number of rate categories. Default is 8.
#' @param em_tol Numeric. EM convergence tolerance. Default is 1e-6.
#' @param em_maxiter Integer. Maximum EM iterations. Default is 100.
#' @param root_prior Character. Root state prior. Default is "fitzjohn".
#' @param seed Integer or NULL. Random seed.
#'
#' @return An object of class "ratescapeML".
#' @export
rateCategories <- function(
    tree,
    data,
    Q,
    k_min = 2,
    k_max = 8,
    em_tol = 1e-6,
    em_maxiter = 100,
    root_prior = "fitzjohn",
    seed = NULL) {

  if (!inherits(tree, "phylo")) stop("tree must be an object of class 'phylo'")
  if (!is.data.frame(data) && !is.matrix(data)) stop("data must be a data frame or matrix")

  ntips <- length(tree$tip.label)
  if (nrow(data) != ntips) stop("Number of rows in data must match number of tips")

  states <- as.numeric(data[, 1])

  if (!is.null(rownames(data))) {
    match_idx <- match(tree$tip.label, rownames(data))
    if (any(is.na(match_idx))) stop("Some tree tip labels do not match data row names")
    states <- states[match_idx]
  }

  states <- as.integer(states)
  if (min(states) > 0) states <- states - min(states)

  if (!is.matrix(Q) || nrow(Q) != ncol(Q)) stop("Q must be a square matrix")
  if (k_min < 1 || k_max < k_min) stop("k_min must be >= 1 and k_max >= k_min")

  root_prior <- match.arg(root_prior, c("fitzjohn", "equal", "stationary"))
  if (!is.null(seed)) set.seed(seed)

  nedges <- nrow(tree$edge)

  k_values <- k_min:k_max
  bic_scores <- numeric(length(k_values))
  loglik_values <- numeric(length(k_values))
  models <- vector("list", length(k_values))

  message(sprintf("RateScape ML: Testing k = %d to %d categories", k_min, k_max))

  for (i in seq_along(k_values)) {
    k_cat <- k_values[i]
    message(sprintf("  Testing k = %d...", k_cat))

    fit_k <- fit_em_gamma(
      tree = tree, states = states, Q = Q,
      k_categories = k_cat, em_tol = em_tol,
      em_maxiter = em_maxiter, root_prior = root_prior
    )

    models[[i]] <- fit_k
    loglik_values[i] <- fit_k$loglik
    # Parameters: (k-1) free weights + 1 shape parameter = k total
    n_params <- k_cat
    n_obs <- nedges

    bic_scores[i] <- -2 * fit_k$loglik + n_params * log(n_obs)
    message(sprintf("    loglik = %.2f, BIC = %.2f", fit_k$loglik, bic_scores[i]))
  }

  best_idx <- which.min(bic_scores)
  best_k <- k_values[best_idx]
  best_fit <- models[[best_idx]]

  message(sprintf("\nBest fit: k = %d (BIC = %.2f)", best_k, bic_scores[best_idx]))

  result <- list(
    best_fit = best_fit,
    k_tested = k_values,
    bic_scores = bic_scores,
    loglik_values = loglik_values,
    best_k = best_k,
    loglik = best_fit$loglik,
    n_params = best_k,
    tree = tree,
    data = states,
    Q = Q,
    call = match.call(),
    root_prior = root_prior,
    em_settings = list(em_tol = em_tol, em_maxiter = em_maxiter)
  )

  class(result) <- c("ratescapeML", "list")
  return(result)
}


#' Fit Gamma Mixture via EM Algorithm
#'
#' @keywords internal
fit_em_gamma <- function(tree, states, Q, k_categories, em_tol, em_maxiter, root_prior) {

  nedges <- nrow(tree$edge)

  # Initialize: gamma shape=1 gives exponential, discretize into quantiles
  alpha <- 1.0

  # Category midpoints from gamma quantiles
  get_rates <- function(alpha, k) {
    if (k == 1) return(1.0)
    probs <- ((1:k) - 0.5) / k
    rates <- qgamma(probs, shape = alpha, rate = alpha)  # rate=alpha so mean=1
    rates / mean(rates)  # Normalize to mean 1
  }

  rates <- get_rates(alpha, k_categories)
  weights <- rep(1 / k_categories, k_categories)

  # Pre-compute full tree likelihood for each rate category (uniform scalars)
  compute_cat_loglik <- function(rate_val) {
    r_scalars <- rep(rate_val, nedges)
    compute_likelihood(tree = tree, data = states, Q = Q,
                       r_scalars = r_scalars, root_prior = root_prior)
  }

  loglik_old <- -Inf

  for (em_iter in 1:em_maxiter) {

    # E-step: compute log-likelihood under each category for the full tree
    cat_logliks <- sapply(rates, compute_cat_loglik)

    # Convert to log-responsibilities
    # log P(cat j | data) prop to log(w_j) + loglik_j
    log_weighted <- log(weights + 1e-300) + cat_logliks
    log_max <- max(log_weighted)
    weighted_liks <- exp(log_weighted - log_max)
    total <- sum(weighted_liks)

    if (total <= 0 || is.nan(total)) {
      message(sprintf("    EM stopped: numerical issues at iteration %d", em_iter))
      break
    }

    responsibility <- weighted_liks / total
    loglik_new <- log_max + log(total)

    # Check convergence
    if (abs(loglik_new - loglik_old) < em_tol) {
      message(sprintf("    EM converged at iteration %d", em_iter))
      break
    }
    loglik_old <- loglik_new

    # M-step: update weights
    weights <- responsibility
    weights <- pmax(weights, 1e-10)
    weights <- weights / sum(weights)

    # M-step: update alpha via numerical optimization
    # Maximize expected complete-data log-likelihood w.r.t. alpha
    obj_fn <- function(log_alpha) {
      a <- exp(log_alpha)
      r <- get_rates(a, k_categories)
      ll <- sapply(r, compute_cat_loglik)
      -sum(responsibility * ll)
    }

    opt <- tryCatch(
      optimize(obj_fn, interval = c(log(0.01), log(100))),
      error = function(e) NULL
    )

    if (!is.null(opt)) {
      alpha <- exp(opt$minimum)
      rates <- get_rates(alpha, k_categories)
    }
  }

  list(
    rates = rates,
    weights = weights,
    alpha = alpha,
    loglik = loglik_old,
    responsibility = responsibility,
    branch_rates = rates  # Per-category rates
  )
}
