#' Construct a Structured Rate Matrix
#'
#' Helper function for building Q matrices from biological parameters.
#' Supports common model structures for chromosome number evolution,
#' general Mk models, and custom parameterizations.
#'
#' @param type Character specifying model type. One of:
#'   \describe{
#'     \item{"mk"}{Equal-rates Mk model (Lewis 2001)}
#'     \item{"sym"}{Symmetric (SYM) model}
#'     \item{"ard"}{All-rates-different (ARD) model}
#'     \item{"chromosome"}{Chromosome number model with gain, loss, polyploidy}
#'     \item{"custom"}{User-specified Q matrix}
#'   }
#' @param nstates Number of character states (required for mk, sym, ard)
#' @param params Named list of model parameters. Required elements depend on type:
#'   \describe{
#'     \item{mk}{\code{rate}: single transition rate}
#'     \item{sym}{\code{rates}: vector of k*(k-1)/2 rates}
#'     \item{ard}{\code{rates}: vector of k*(k-1) rates, filled by column}
#'     \item{chromosome}{\code{gain}: ascending dysploidy rate;
#'       \code{loss}: descending dysploidy rate;
#'       \code{polyploidy}: polyploidization rate (genome doubling);
#'       \code{max_chrom}: maximum chromosome number modeled}
#'   }
#' @param Q For type = "custom", a matrix to validate and return.
#' @return A valid k x k rate matrix
#' @export
#' @examples
#' # 4-state equal rates model
#' Q_mk <- makeQ("mk", nstates = 4, params = list(rate = 0.1))
#'
#' # Chromosome number model
#' Q_chr <- makeQ("chromosome", params = list(
#'   gain = 0.2, loss = 0.3, polyploidy = 0.01, max_chrom = 30
#' ))
makeQ <- function(type = c("mk", "sym", "ard", "chromosome", "custom"),
                  nstates = NULL, params = list(), Q = NULL) {

  type <- match.arg(type)

  if (type == "mk") {
    if (is.null(nstates)) stop("nstates required for mk model")
    rate <- params$rate
    if (is.null(rate)) stop("params$rate required for mk model")
    Q_mat <- matrix(rate, nstates, nstates)
    diag(Q_mat) <- 0
    diag(Q_mat) <- -rowSums(Q_mat)
    return(Q_mat)
  }

  if (type == "sym") {
    if (is.null(nstates)) stop("nstates required for sym model")
    n_rates <- nstates * (nstates - 1) / 2
    if (length(params$rates) != n_rates) {
      stop("params$rates must have length k*(k-1)/2 = ", n_rates)
    }
    Q_mat <- matrix(0, nstates, nstates)
    Q_mat[upper.tri(Q_mat)] <- params$rates
    Q_mat <- Q_mat + t(Q_mat)
    diag(Q_mat) <- -rowSums(Q_mat)
    return(Q_mat)
  }

  if (type == "ard") {
    if (is.null(nstates)) stop("nstates required for ard model")
    n_rates <- nstates * (nstates - 1)
    if (length(params$rates) != n_rates) {
      stop("params$rates must have length k*(k-1) = ", n_rates)
    }
    Q_mat <- matrix(0, nstates, nstates)
    Q_mat[row(Q_mat) != col(Q_mat)] <- params$rates
    diag(Q_mat) <- -rowSums(Q_mat)
    return(Q_mat)
  }

  if (type == "chromosome") {
    gain <- params$gain
    loss <- params$loss
    polyploidy <- params$polyploidy
    max_chrom <- params$max_chrom
    if (is.null(gain) || is.null(loss) || is.null(max_chrom)) {
      stop("params must include gain, loss, and max_chrom")
    }
    if (is.null(polyploidy)) polyploidy <- 0

    k <- max_chrom
    Q_mat <- matrix(0, k, k)

    for (i in 1:k) {
      # Gain: i -> i+1
      if (i < k) Q_mat[i, i + 1] <- gain
      # Loss: i -> i-1
      if (i > 1) Q_mat[i, i - 1] <- loss
      # Polyploidy: i -> 2i (if 2i <= max_chrom)
      if (polyploidy > 0 && 2 * i <= k) {
        Q_mat[i, 2 * i] <- polyploidy
      }
    }

    diag(Q_mat) <- -rowSums(Q_mat)
    return(Q_mat)
  }

  if (type == "custom") {
    if (is.null(Q)) stop("Q matrix required for custom type")
    if (!is.matrix(Q)) stop("Q must be a matrix")
    # Validate
    if (any(Q[row(Q) != col(Q)] < 0)) {
      stop("Off-diagonal elements of Q must be non-negative")
    }
    row_sums <- rowSums(Q)
    if (any(abs(row_sums) > 1e-10)) {
      warning("Rows of Q do not sum to zero; correcting diagonals")
      diag(Q) <- 0
      diag(Q) <- -rowSums(Q)
    }
    return(Q)
  }
}


#' Create a Q-matrix function for joint estimation
#'
#' Returns a function that maps a parameter vector to a Q matrix,
#' suitable for use with \code{fitRateScape(Q = ...)} for joint estimation.
#'
#' @param type Model type (see \code{\link{makeQ}})
#' @param nstates Number of states (if applicable)
#' @param max_chrom Maximum chromosome number (for chromosome model)
#' @return A function that takes a numeric parameter vector and returns a Q matrix
#' @export
makeQ_func <- function(type = c("mk", "chromosome"),
                       nstates = NULL, max_chrom = NULL) {

  type <- match.arg(type)

  if (type == "mk") {
    if (is.null(nstates)) stop("nstates required")
    f <- function(par) {
      # par[1] = rate
      makeQ("mk", nstates = nstates, params = list(rate = exp(par[1])))
    }
    attr(f, "n_params") <- 1
    attr(f, "param_names") <- "log_rate"
    return(f)
  }

  if (type == "chromosome") {
    if (is.null(max_chrom)) stop("max_chrom required")
    f <- function(par) {
      # par = c(log_gain, log_loss, log_polyploidy)
      makeQ("chromosome", params = list(
        gain = exp(par[1]),
        loss = exp(par[2]),
        polyploidy = exp(par[3]),
        max_chrom = max_chrom
      ))
    }
    attr(f, "n_params") <- 3
    attr(f, "param_names") <- c("log_gain", "log_loss", "log_polyploidy")
    return(f)
  }
}
