// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <unordered_map>

using namespace Rcpp;
using namespace arma;


// ============================================================
// Matrix exponential via Padé approximation (scaling & squaring)
// ============================================================
// We implement our own to avoid R-level calls to expm::expm
// inside tight C++ loops. Uses the degree-13 Padé approximant
// with scaling and squaring (Higham 2005).

// Padé coefficients for degree 13
static const double pade_coeff[] = {
  1.0,
  0.5,
  0.12,
  1.833333333333333333e-02,
  1.992753623188405797e-03,
  1.630434782608695652e-04,
  1.035196687370600414e-05,
  5.175983436853002070e-07,
  2.043151356652500817e-08,
  6.306022705717595115e-10,
  1.483770048404140027e-11,
  2.529153491597965955e-13,
  2.810170546219962172e-15,
  1.544049750670308885e-17
};

mat matrix_expm(const mat& A) {
  int n = A.n_rows;

  // Compute norms to decide scaling
  double normA = norm(A, "inf");

  if (normA == 0.0) {
    return eye<mat>(n, n);
  }

  // Scaling: find s such that ||A/2^s|| < theta_13 = 5.37
  int s = std::max(0, (int)std::ceil(std::log2(normA / 5.371920351148152)));
  mat As = A / std::pow(2.0, s);

  // Compute Padé approximant of degree 13
  mat As2 = As * As;
  mat As4 = As2 * As2;
  mat As6 = As4 * As2;

  mat U = As * (As6 * (pade_coeff[13] * As6 + pade_coeff[11] * As4 + pade_coeff[9] * As2) +
                pade_coeff[7] * As6 + pade_coeff[5] * As4 + pade_coeff[3] * As2 +
                pade_coeff[1] * eye<mat>(n, n));

  mat V = As6 * (pade_coeff[12] * As6 + pade_coeff[10] * As4 + pade_coeff[8] * As2) +
          pade_coeff[6] * As6 + pade_coeff[4] * As4 + pade_coeff[2] * As2 +
          pade_coeff[0] * eye<mat>(n, n);

  mat P = V + U;   // numerator
  mat Qm = V - U;  // denominator

  // Solve Q * result = P  =>  result = Q^{-1} P
  mat result = solve(Qm, P);

  // Squaring phase
  for (int i = 0; i < s; i++) {
    result = result * result;
  }

  return result;
}


// ============================================================
// Cache for transition probability matrices
// ============================================================
// Key: discretized (rate * branch_length) to avoid recomputing
// identical P matrices. Essential when many branches share the
// same rate category (discretized gamma) or are in the spike.

struct PmatCache {
  std::unordered_map<long long, mat> cache;
  const mat& Q;
  int nstates;
  double precision;

  PmatCache(const mat& Q_, double prec = 1e-8)
    : Q(Q_), nstates(Q_.n_rows), precision(prec) {}

  // Discretize t to a cache key
  long long make_key(double t) {
    return (long long)(t / precision + 0.5);
  }

  const mat& get(double t) {
    long long key = make_key(t);
    auto it = cache.find(key);
    if (it != cache.end()) {
      return it->second;
    }
    mat P = matrix_expm(Q * t);
    // Clamp small negatives from numerical error
    P.clamp(0.0, datum::inf);
    // Renormalize rows
    for (int i = 0; i < nstates; i++) {
      double rs = accu(P.row(i));
      if (rs > 0) P.row(i) /= rs;
    }
    cache[key] = P;
    return cache[key];
  }

  void clear() { cache.clear(); }
};


// ============================================================
// Pruning likelihood (scaled) — C++ implementation
// ============================================================

// [[Rcpp::export]]
double pruning_llik_cpp(IntegerVector parent_vec,
                        IntegerVector child_vec,
                        NumericVector edge_lengths,
                        IntegerVector tip_states,  // 1-indexed
                        NumericMatrix Q_mat,
                        NumericVector rates,
                        NumericVector root_prior,
                        int ntip,
                        int nnode) {

  int nedge = parent_vec.size();
  int nstates = Q_mat.nrow();
  int nnodes_total = ntip + nnode;

  // Convert Q to arma::mat
  mat Q(Q_mat.begin(), nstates, nstates, false);

  // Initialize P-matrix cache
  PmatCache pcache(Q);

  // Partial likelihoods: nstates x nnodes_total
  mat partial_lk(nstates, nnodes_total, fill::zeros);

  // Log scaling factors per node
  vec log_scale(nnodes_total, fill::zeros);

  // Set tip likelihoods
  for (int i = 0; i < ntip; i++) {
    int state = tip_states[i] - 1;  // Convert to 0-indexed
    partial_lk(state, i) = 1.0;
  }

  // Track accumulated children contributions
  // We need to multiply conditional likelihoods from all children
  std::vector<vec> node_accum(nnodes_total);
  std::vector<double> node_logscale(nnodes_total, 0.0);
  std::vector<bool> node_initialized(nnodes_total, false);

  // Traverse edges in postorder (assumed to be in postorder already)
  for (int e = 0; e < nedge; e++) {
    int parent = parent_vec[e] - 1;  // 0-indexed
    int child = child_vec[e] - 1;
    double eff_t = rates[e] * edge_lengths[e];

    // Get transition probability matrix (cached)
    const mat& P = pcache.get(eff_t);

    // Conditional likelihood: P * partial_lk[child]
    vec child_lk = partial_lk.col(child);
    vec cond_lk = P * child_lk;

    // Multiply into parent's accumulator
    if (!node_initialized[parent]) {
      node_accum[parent] = cond_lk;
      node_logscale[parent] = log_scale[child];
      node_initialized[parent] = true;
    } else {
      node_accum[parent] %= cond_lk;  // element-wise multiply
      node_logscale[parent] += log_scale[child];
    }

    // Check if this is the last child of parent
    // (by checking remaining edges with same parent after this one)
    bool is_last_child = true;
    for (int f = e + 1; f < nedge; f++) {
      if (parent_vec[f] - 1 == parent) {
        is_last_child = false;
        break;
      }
    }

    if (is_last_child) {
      // Scale to avoid underflow
      double max_lk = node_accum[parent].max();
      if (max_lk > 0 && std::isfinite(max_lk)) {
        partial_lk.col(parent) = node_accum[parent] / max_lk;
        log_scale[parent] = std::log(max_lk) + node_logscale[parent];
      } else {
        partial_lk.col(parent) = node_accum[parent];
        log_scale[parent] = node_logscale[parent];
      }
    }
  }

  // Root likelihood
  int root = ntip;  // 0-indexed root node
  vec root_lk = partial_lk.col(root);

  // Apply root prior
  double llik = 0.0;
  for (int s = 0; s < nstates; s++) {
    llik += root_lk[s] * root_prior[s];
  }

  if (llik > 0) {
    return std::log(llik) + log_scale[root];
  } else {
    return -1e300;  // -Inf substitute
  }
}


// ============================================================
// Full Gibbs sweep for spike-and-slab indicators + rate updates
// ============================================================
// This does ONE full pass over all branches, updating z_i and r_i
// for each branch. Returns the updated rates, indicators, and
// log-likelihood. This is the critical inner loop.

// [[Rcpp::export]]
List gibbs_sweep_cpp(IntegerVector parent_vec,
                     IntegerVector child_vec,
                     NumericVector edge_lengths,
                     IntegerVector tip_states,
                     NumericMatrix Q_mat,
                     NumericVector rates,          // current rates
                     IntegerVector z_indicators,   // current indicators
                     double pi_val,                // current pi
                     double sigma2,                // current slab variance
                     NumericVector root_prior,
                     int ntip,
                     int nnode,
                     double rate_proposal_sd) {

  int nedge = rates.size();
  int nstates = Q_mat.nrow();

  // Copy inputs to modify
  NumericVector new_rates = clone(rates);
  IntegerVector new_z = clone(z_indicators);

  // Current full log-likelihood
  double llik_current = pruning_llik_cpp(parent_vec, child_vec, edge_lengths,
                                          tip_states, Q_mat, new_rates,
                                          root_prior, ntip, nnode);

  int n_accepted = 0;
  int n_proposed = 0;
  double sd_log = std::sqrt(sigma2);

  for (int e = 0; e < nedge; e++) {

    // === Update z_i (Gibbs) ===
    double saved_rate = new_rates[e];

    // Log-likelihood under spike (r = 1)
    new_rates[e] = 1.0;
    double llik_spike = pruning_llik_cpp(parent_vec, child_vec, edge_lengths,
                                          tip_states, Q_mat, new_rates,
                                          root_prior, ntip, nnode);

    // For slab: if currently in slab, use current rate; else draw one
    double slab_rate;
    if (new_z[e] == 1) {
      slab_rate = saved_rate;
    } else {
      // Draw from slab prior
      slab_rate = std::exp(R::rnorm(0.0, sd_log));
    }

    new_rates[e] = slab_rate;
    double llik_slab = pruning_llik_cpp(parent_vec, child_vec, edge_lengths,
                                         tip_states, Q_mat, new_rates,
                                         root_prior, ntip, nnode);

    // Log prior for slab rate
    double log_prior_slab = R::dlnorm(slab_rate, 0.0, sd_log, 1);

    // Log posterior probabilities
    double log_prob_spike = std::log(pi_val) + llik_spike;
    double log_prob_slab = std::log(1.0 - pi_val) + llik_slab + log_prior_slab;

    // Normalize
    double max_log = std::max(log_prob_spike, log_prob_slab);
    double prob_spike = std::exp(log_prob_spike - max_log) /
      (std::exp(log_prob_spike - max_log) + std::exp(log_prob_slab - max_log));

    if (R::runif(0.0, 1.0) < prob_spike) {
      new_z[e] = 0;
      new_rates[e] = 1.0;
      llik_current = llik_spike;
    } else {
      new_z[e] = 1;
      new_rates[e] = slab_rate;
      llik_current = llik_slab;
    }

    // === Update r_i if in slab (MH) ===
    if (new_z[e] == 1) {
      n_proposed++;
      double log_r_current = std::log(new_rates[e]);
      double log_r_proposed = R::rnorm(log_r_current, rate_proposal_sd);
      double r_proposed = std::exp(log_r_proposed);

      double saved_rate2 = new_rates[e];
      new_rates[e] = r_proposed;

      double llik_proposed = pruning_llik_cpp(parent_vec, child_vec, edge_lengths,
                                               tip_states, Q_mat, new_rates,
                                               root_prior, ntip, nnode);

      double log_prior_curr = R::dlnorm(saved_rate2, 0.0, sd_log, 1);
      double log_prior_prop = R::dlnorm(r_proposed, 0.0, sd_log, 1);

      double log_alpha = (llik_proposed - llik_current) +
        (log_prior_prop - log_prior_curr);

      if (std::log(R::runif(0.0, 1.0)) < log_alpha) {
        // Accept
        llik_current = llik_proposed;
        n_accepted++;
      } else {
        // Reject: revert
        new_rates[e] = saved_rate2;
      }
    }
  }

  return List::create(
    Named("rates") = new_rates,
    Named("z") = new_z,
    Named("loglik") = llik_current,
    Named("n_accepted") = n_accepted,
    Named("n_proposed") = n_proposed
  );
}


// ============================================================
// Batch likelihood computation for EM rate categories
// ============================================================
// For each branch, compute the likelihood under each rate category.
// Returns a matrix (nedge x ncat) of log-likelihoods.

// [[Rcpp::export]]
NumericMatrix batch_branch_llik_cpp(IntegerVector parent_vec,
                                     IntegerVector child_vec,
                                     NumericVector edge_lengths,
                                     IntegerVector tip_states,
                                     NumericMatrix Q_mat,
                                     NumericVector base_rates,
                                     NumericVector rate_categories,
                                     NumericVector root_prior,
                                     int ntip,
                                     int nnode) {

  int nedge = parent_vec.size();
  int ncat = rate_categories.size();

  NumericMatrix result(nedge, ncat);

  for (int e = 0; e < nedge; e++) {
    for (int c = 0; c < ncat; c++) {
      // Set branch e to rate category c, all others at base_rates
      NumericVector test_rates = clone(base_rates);
      test_rates[e] = rate_categories[c];

      result(e, c) = pruning_llik_cpp(parent_vec, child_vec, edge_lengths,
                                       tip_states, Q_mat, test_rates,
                                       root_prior, ntip, nnode);
    }
  }

  return result;
}
