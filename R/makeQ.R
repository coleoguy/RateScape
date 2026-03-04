#' Construct Discrete Character Rate Matrices
#'
#' Creates transition rate matrices (Q matrices) for discrete character evolution
#' models. Supports standard models (Mk, SYM, ARD), constrained models (ordered,
#' Dollo/irreversible), and custom specifications with user-defined zero
#' constraints.
#'
#' @param tree An object of class "phylo" (optional, used for determining number of states).
#' @param model Character. Type of rate matrix model:
#'   \describe{
#'     \item{"mk"}{Mk model. All transitions between different states have equal rate.}
#'     \item{"sym"}{Symmetric model (SYM). Transitions are symmetric: Q_ij = Q_ji.}
#'     \item{"ard"}{All-rates-different (ARD). Each transition has its own rate.}
#'     \item{"ordered"}{Ordered model. Only transitions between adjacent states
#'       (i <-> i+1) are allowed. Appropriate when states have a natural ordering
#'       (e.g., 0 < 1 < 2 < 3). Both forward and backward transitions allowed.}
#'     \item{"dollo"}{Dollo (irreversible) model. Transitions only go in one
#'       direction: 0 -> 1 -> 2 -> ... Reversals are forbidden. Useful for
#'       complex traits assumed to be lost but never regained.}
#'     \item{"irreversible"}{Alias for "dollo".}
#'     \item{"meristic"}{Meristic model. Ordered model where forward and reverse
#'       rates between adjacent states are equal (symmetric ordered).}
#'     \item{"chromosome"}{Specialized model for chromosome number evolution with gains,
#'       losses, and polyploidy events.}
#'     \item{"custom"}{User-provides the matrix directly via the Q parameter.}
#'   }
#' @param k Integer. Number of states (character states or chromosome numbers).
#'   Required unless model = "custom" or tree is provided.
#' @param Q Matrix. User-specified transition matrix (required if model = "custom").
#' @param polyploid Logical. For chromosome model, whether to include polyploidy
#'   transitions (doubling). Default is TRUE.
#' @param constraints A list or matrix specifying zero constraints on the Q matrix.
#'   Can be:
#'   \describe{
#'     \item{matrix}{A k x k logical or 0/1 matrix where TRUE/1 means the
#'       transition is ALLOWED. FALSE/0 entries will be forced to zero. Diagonal
#'       is ignored (always set from row sums).}
#'     \item{character}{A named constraint style:
#'       "no_diagonal_jumps" forbids transitions that skip states (e.g., 0->2);
#'       "upper_triangular" allows only forward (ascending) transitions;
#'       "lower_triangular" allows only backward (descending) transitions.}
#'   }
#' @param rate_classes A k x k integer matrix mapping each transition (i, j) to
#'   a rate class index. Transitions sharing the same class index are constrained
#'   to have the same rate during optimization. Zero entries in rate_classes mean
#'   the transition is forbidden. This provides maximum flexibility for custom
#'   rate structures (e.g., gains faster than losses, transitions between certain
#'   states sharing rates).
#'
#' @details
#'
#' **Standard Models:**
#'
#' The Mk, SYM, and ARD models are standard in phylogenetic comparative biology.
#' They differ in the number of free rate parameters.
#'
#' **Ordered Model:** Only adjacent transitions are allowed (Q_ij = 0 if |i-j| > 1).
#' This is appropriate for characters with a natural ordering, such as the number
#' of petals, digit count, or body size categories.
#'
#' **Dollo/Irreversible Model:** Transitions go 0->1->2->...->k-1 but never
#' reverse. This implements Dollo's law (complex features, once lost, are not
#' regained). The classic use case is for irreversible losses of complex
#' structures (eyes, wings, teeth).
#'
#' **Meristic Model:** Like ordered, but forward and backward rates between the
#' same pair of adjacent states are equal. Useful for characters like chromosome
#' number where gains and losses of one unit are mechanistically similar.
#'
#' **Custom Constraints:** You can apply zero constraints on top of any base
#' model. For example, start with ARD and apply a constraints matrix that forbids
#' specific biologically implausible transitions.
#'
#' **Rate Classes:** The most flexible option. You define groups of transitions
#' that must share the same rate, while other transitions are forbidden (zero).
#' This lets you encode biological hypotheses directly into the rate matrix
#' structure.
#'
#' @return A transition rate matrix (Q matrix) with dimensions k x k. Includes
#'   attributes:
#'   \describe{
#'     \item{param_map}{Integer matrix mapping each (i,j) entry to a rate class.}
#'     \item{n_params}{Number of free rate parameters.}
#'     \item{model}{The model name used.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Standard models
#'   Q_mk <- makeQ(model = "mk", k = 3)
#'   Q_ard <- makeQ(model = "ard", k = 4)
#'
#'   # Ordered: only adjacent transitions
#'   Q_ord <- makeQ(model = "ordered", k = 5)
#'
#'   # Dollo (irreversible): forward only
#'   Q_dollo <- makeQ(model = "dollo", k = 3)
#'
#'   # ARD with custom zero constraints
#'   allowed <- matrix(TRUE, 4, 4)
#'   diag(allowed) <- FALSE
#'   allowed[1, 4] <- FALSE  # forbid 0->3
#'   allowed[4, 1] <- FALSE  # forbid 3->0
#'   Q_const <- makeQ(model = "ard", k = 4, constraints = allowed)
#'
#'   # Rate classes: gains share one rate, losses share another
#'   rc <- matrix(0, 4, 4)
#'   rc[1,2] <- rc[2,3] <- rc[3,4] <- 1  # forward = class 1
#'   rc[2,1] <- rc[3,2] <- rc[4,3] <- 2  # backward = class 2
#'   Q_2rate <- makeQ(model = "custom", k = 4, rate_classes = rc)
#' }
#'
#' @export
makeQ <- function(
    tree = NULL,
    model = "mk",
    k = NULL,
    Q = NULL,
    polyploid = TRUE,
    constraints = NULL,
    rate_classes = NULL) {

  # Validate model
  model <- match.arg(model, c("mk", "sym", "ard", "ordered", "dollo",
                                "irreversible", "meristic", "chromosome", "custom"))

  # Handle aliases

  if (model == "irreversible") model <- "dollo"

  # Custom Q matrix
  if (model == "custom" && !is.null(Q) && is.null(rate_classes)) {
    if (!is.matrix(Q) || nrow(Q) != ncol(Q)) {
      stop("Q must be a square matrix")
    }
    # If constraints provided, apply them
    if (!is.null(constraints)) {
      Q <- apply_constraints(Q, constraints)
    }
    return(Q)
  }

  # Rate classes (maximum flexibility path)
  if (!is.null(rate_classes)) {
    if (is.null(k)) k <- nrow(rate_classes)
    return(build_Q_from_rate_classes(rate_classes, k))
  }

  # Determine k
  if (is.null(k)) {
    if (!is.null(tree) && inherits(tree, "phylo")) {
      warning("k not provided; using default k = 2")
      k <- 2
    } else {
      stop("k must be specified for model = '", model, "'")
    }
  }
  if (k < 2) stop("k must be at least 2")

  # Chromosome model
  if (model == "chromosome") {
    return(make_chromosome_Q(k, polyploid = polyploid))
  }

  # Build base Q matrix and parameter map
  Q <- matrix(0, nrow = k, ncol = k)
  param_map <- matrix(0L, nrow = k, ncol = k)

  if (model == "mk") {
    # All off-diagonals = 1 (single rate class)
    Q[upper.tri(Q)] <- 1
    Q[lower.tri(Q)] <- 1
    param_map[Q > 0] <- 1L
    n_params <- 1

  } else if (model == "sym") {
    # Symmetric: Q_ij = Q_ji, each pair shares a rate
    idx <- 0
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        idx <- idx + 1
        Q[i, j] <- idx
        Q[j, i] <- idx
        param_map[i, j] <- idx
        param_map[j, i] <- idx
      }
    }
    n_params <- idx

  } else if (model == "ard") {
    # All different
    idx <- 0
    for (i in 1:k) {
      for (j in 1:k) {
        if (i != j) {
          idx <- idx + 1
          Q[i, j] <- 1
          param_map[i, j] <- idx
        }
      }
    }
    n_params <- idx

  } else if (model == "ordered") {
    # Only adjacent transitions: i <-> i+1
    idx <- 0
    for (i in 1:(k-1)) {
      idx <- idx + 1
      Q[i, i+1] <- 1   # forward
      param_map[i, i+1] <- idx
      idx <- idx + 1
      Q[i+1, i] <- 1   # backward
      param_map[i+1, i] <- idx
    }
    n_params <- idx

  } else if (model == "dollo") {
    # Irreversible: 0->1, 1->2, ..., (k-2)->(k-1). No reversals.
    idx <- 0
    for (i in 1:(k-1)) {
      idx <- idx + 1
      Q[i, i+1] <- 1
      param_map[i, i+1] <- idx
    }
    n_params <- idx

  } else if (model == "meristic") {
    # Ordered + symmetric between adjacent pairs
    idx <- 0
    for (i in 1:(k-1)) {
      idx <- idx + 1
      Q[i, i+1] <- 1
      Q[i+1, i] <- 1
      param_map[i, i+1] <- idx
      param_map[i+1, i] <- idx
    }
    n_params <- idx
  }

  # Apply zero constraints if provided
  if (!is.null(constraints)) {
    result <- apply_constraints_to_parameterized(Q, param_map, constraints, k)
    Q <- result$Q
    param_map <- result$param_map
    n_params <- max(param_map)
  }

  # Set diagonal so row sums = 0
  diag(Q) <- 0
  Q[Q > 0] <- 1  # normalize non-zero entries to 1 for initial values
  diag(Q) <- -rowSums(Q)

  attr(Q, "param_map") <- param_map
  attr(Q, "n_params") <- n_params
  attr(Q, "model") <- model
  return(Q)
}


#' Build Q matrix from rate class specification
#'
#' @keywords internal
build_Q_from_rate_classes <- function(rate_classes, k) {
  if (!is.matrix(rate_classes) || nrow(rate_classes) != k || ncol(rate_classes) != k)
    stop("rate_classes must be a k x k matrix")

  Q <- matrix(0, nrow = k, ncol = k)
  param_map <- matrix(0L, nrow = k, ncol = k)

  for (i in 1:k) {
    for (j in 1:k) {
      if (i != j && rate_classes[i, j] > 0) {
        Q[i, j] <- 1
        param_map[i, j] <- as.integer(rate_classes[i, j])
      }
    }
  }

  n_params <- max(param_map)
  diag(Q) <- -rowSums(Q)

  attr(Q, "param_map") <- param_map
  attr(Q, "n_params") <- n_params
  attr(Q, "model") <- "rate_classes"
  Q
}


#' Apply constraints matrix to a base Q
#'
#' @keywords internal
apply_constraints <- function(Q, constraints) {
  k <- nrow(Q)
  if (is.character(constraints)) {
    constraints <- match.arg(constraints, c("no_diagonal_jumps",
                                             "upper_triangular",
                                             "lower_triangular"))
    mask <- matrix(TRUE, k, k)
    if (constraints == "no_diagonal_jumps") {
      for (i in 1:k) for (j in 1:k)
        if (abs(i - j) > 1) mask[i, j] <- FALSE
    } else if (constraints == "upper_triangular") {
      mask[lower.tri(mask)] <- FALSE
    } else if (constraints == "lower_triangular") {
      mask[upper.tri(mask)] <- FALSE
    }
    diag(mask) <- FALSE
  } else {
    mask <- as.matrix(constraints)
    if (nrow(mask) != k || ncol(mask) != k)
      stop("constraints matrix must be k x k")
    mask <- mask > 0
    diag(mask) <- FALSE
  }

  Q[!mask] <- 0
  diag(Q) <- -rowSums(Q)
  Q
}


#' Apply constraints to parameterized Q matrix
#' @keywords internal
apply_constraints_to_parameterized <- function(Q, param_map, constraints, k) {
  if (is.character(constraints)) {
    constraints <- match.arg(constraints, c("no_diagonal_jumps",
                                             "upper_triangular",
                                             "lower_triangular"))
    mask <- matrix(TRUE, k, k)
    if (constraints == "no_diagonal_jumps") {
      for (i in 1:k) for (j in 1:k)
        if (abs(i - j) > 1) mask[i, j] <- FALSE
    } else if (constraints == "upper_triangular") {
      mask[lower.tri(mask)] <- FALSE
    } else if (constraints == "lower_triangular") {
      mask[upper.tri(mask)] <- FALSE
    }
    diag(mask) <- FALSE
  } else {
    mask <- as.matrix(constraints)
    if (nrow(mask) != k || ncol(mask) != k)
      stop("constraints matrix must be k x k")
    mask <- mask > 0
    diag(mask) <- FALSE
  }

  # Zero out forbidden transitions
  Q[!mask] <- 0
  param_map[!mask] <- 0L

  # Re-index parameters (remove gaps)
  old_ids <- sort(unique(as.vector(param_map[param_map > 0])))
  new_map <- param_map
  for (new_idx in seq_along(old_ids)) {
    new_map[param_map == old_ids[new_idx]] <- new_idx
  }

  diag(Q) <- -rowSums(Q)
  list(Q = Q, param_map = new_map)
}


#' Create chromosome number evolution Q matrix
#'
#' @keywords internal
make_chromosome_Q <- function(k, polyploid = TRUE) {
  Q <- matrix(0, nrow = k, ncol = k)

  for (i in 1:(k - 1)) {
    # Gain: i -> i+1
    Q[i, i + 1] <- 1.0

    # Loss: i -> i-1 (if not at minimum)
    if (i > 1) {
      Q[i, i - 1] <- 1.0
    }

    # Polyploidy: i -> 2i (if within bounds)
    if (polyploid && 2 * i <= k) {
      Q[i, 2 * i] <- 1.0
    }
  }

  # Handle losses for higher states
  for (i in 2:k) {
    if (Q[i, i - 1] == 0) {
      Q[i, i - 1] <- 1.0
    }
  }

  # Set diagonals
  diag(Q) <- -rowSums(Q)

  # Parameter map: 3 rate classes (gain, loss, polyploidy)
  pm <- matrix(0L, k, k)
  for (i in 1:(k-1)) {
    pm[i, i+1] <- 1L  # gain
  }
  for (i in 2:k) {
    pm[i, i-1] <- 2L  # loss
  }
  if (polyploid) {
    for (i in 1:k) {
      if (2*i <= k) pm[i, 2*i] <- 3L
    }
  }

  attr(Q, "param_map") <- pm
  attr(Q, "n_params") <- if (polyploid) 3L else 2L
  attr(Q, "model") <- "chromosome"
  Q
}
