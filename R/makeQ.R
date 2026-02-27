#' Construct Discrete Character Rate Matrices
#'
#' Creates transition rate matrices (Q matrices) for discrete character evolution
#' models. Supports standard models (Mk, SYM, ARD) and specialized models for
#' chromosome number evolution and custom specifications.
#'
#' @param tree An object of class "phylo" (optional, used for determining number of states).
#' @param model Character. Type of rate matrix model:
#'   - "mk": Mk model. All transitions between different states have equal rate.
#'   - "sym": Symmetric model (SYM). Transitions are symmetric: Q_ij = Q_ji.
#'   - "ard": All-rates-different (ARD). Each transition has its own rate.
#'   - "chromosome": Specialized model for chromosome number evolution with gains,
#'     losses, and polyploidy events.
#'   - "custom": User-provides the matrix directly via the Q parameter.
#' @param k Integer. Number of states (character states or chromosome numbers).
#'   Required unless model = "custom" or tree is provided.
#' @param Q Matrix. User-specified transition matrix (required if model = "custom").
#' @param polyploid Logical. For chromosome model, whether to include polyploidy
#'   transitions (doubling). Default is TRUE.
#'
#' @details
#'
#' **Mk Model:** All off-diagonal entries are equal (Q_ij = μ for i ≠ j).
#' Diagonal entries are set so row sums = 0.
#'
#' **Symmetric (SYM) Model:** Forward and reverse transitions are equal (Q_ij = Q_ji).
#' Commonly used for morphological traits with reversible evolution.
#'
#' **All-Rates-Different (ARD) Model:** Each off-diagonal entry is independent.
#' Maximum complexity; can lead to overparameterization with small datasets.
#'
#' **Chromosome Model:** Specialized for chromosome number evolution:
#'   - Gain: i → i+1 (rate = λ_gain)
#'   - Loss: i → i-1 (rate = λ_loss)
#'   - Polyploid: i → 2i (rate = λ_poly)
#'   - All other transitions have rate 0.
#'
#' **Custom:** User supplies the entire Q matrix.
#'
#' @return A transition rate matrix (Q matrix) with dimensions k × k.
#'   Off-diagonal entries are positive rates. Diagonal entries are negative,
#'   set so that each row sums to 0.
#'
#' @examples
#' \dontrun{
#'   # Mk model with 3 states
#'   Q_mk <- makeQ(model = "mk", k = 3)
#'
#'   # ARD model with 4 states
#'   Q_ard <- makeQ(model = "ard", k = 4)
#'
#'   # Chromosome model for states 4–20
#'   Q_chrom <- makeQ(model = "chromosome", k = 17)
#'
#'   # Custom user matrix
#'   my_Q <- matrix(c(-1, 1, 1, -1), nrow = 2)
#'   Q_custom <- makeQ(model = "custom", Q = my_Q)
#' }
#'
#' @export
makeQ <- function(
    tree = NULL,
    model = "mk",
    k = NULL,
    Q = NULL,
    polyploid = TRUE) {

  # Validate model
  model <- match.arg(model, c("mk", "sym", "ard", "chromosome", "custom"))

  # Determine k
  if (model == "custom") {
    if (is.null(Q)) {
      stop("For model = 'custom', Q matrix must be provided")
    }
    if (!is.matrix(Q) || nrow(Q) != ncol(Q)) {
      stop("Q must be a square matrix")
    }
    return(Q)
  }

  if (is.null(k)) {
    if (!is.null(tree) && inherits(tree, "phylo")) {
      # Try to infer from tree (not typically possible without data)
      # This is a placeholder; usually k must be specified
      warning("k not provided; using default k = 2")
      k <- 2
    } else {
      stop("k must be specified for model = '", model, "'")
    }
  }

  if (k < 2) {
    stop("k must be at least 2")
  }

  if (model == "chromosome") {
    return(make_chromosome_Q(k, polyploid = polyploid))
  }

  # Standard models
  Q <- matrix(0, nrow = k, ncol = k)

  if (model == "mk") {
    # All off-diagonals = 1 (will be rescaled)
    Q[upper.tri(Q)] <- 1
    Q[lower.tri(Q)] <- 1
  } else if (model == "sym") {
    # Symmetric: Q_ij = Q_ji
    # Initialize as Mk, then symmetrize
    Q[upper.tri(Q)] <- 1
    Q[lower.tri(Q)] <- 1
    # Already symmetric for all-ones
  } else if (model == "ard") {
    # All different: assign distinct values for identifiability
    # Use integers 1, 2, 3, ... for off-diagonals
    idx <- 1
    for (i in 1:k) {
      for (j in 1:k) {
        if (i != j) {
          Q[i, j] <- idx
          idx <- idx + 1
        }
      }
    }
  }

  # Set diagonal so row sums = 0
  diag(Q) <- -rowSums(Q)

  return(Q)
}


#' Create chromosome number evolution Q matrix
#'
#' @keywords internal
make_chromosome_Q <- function(k, polyploid = TRUE) {
  Q <- matrix(0, nrow = k, ncol = k)

  for (i in 1:(k - 1)) {
    # Gain: i → i+1
    Q[i, i + 1] <- 1.0

    # Loss: i → i-1 (if not at minimum)
    if (i > 1) {
      Q[i, i - 1] <- 1.0
    }

    # Polyploidy: i → 2i (if within bounds)
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

  return(Q)
}
