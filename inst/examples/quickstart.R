# RateScape Quick Start Example
# Demonstrates workflow for detecting rate heterogeneity in discrete character evolution

library(RateScape)
library(ape)

# ============================================================================
# STEP 1: Simulate a tree and character data with rate heterogeneity
# ============================================================================

cat("Step 1: Simulating phylogenetic tree and character data\n")

# Create a random phylogenetic tree with 50 species
set.seed(42)
tree <- ape::rtree(n = 50, rooted = TRUE)
tree$edge.length <- tree$edge.length / max(tree$edge.length)  # Rescale to [0,1]

# Create a transition rate matrix for a binary trait
Q <- makeQ(model = "mk", k = 2)
cat("Transition rate matrix (Mk model, k=2):\n")
print(Q)

# Simulate character data under spike-and-slab prior
# 70% of branches at background rate, 30% with heterogeneous rates
sim_data <- simRateScape(
  tree = tree,
  Q = Q,
  lambda_sigma = 1.0,           # E[σ²] = 1.0 (standard heterogeneity)
  pi = 0.7,                      # 70% background rate
  nrep = 1,                      # Single replicate
  return_rates = TRUE,           # Also return branch rates
  seed = 42
)

data <- data.frame(state = sim_data$data[1, ])
rownames(data) <- tree$tip.label

cat("Simulated character states at tips (first 10):\n")
print(head(data, 10))


# ============================================================================
# STEP 2: Run prior diagnostic check
# ============================================================================

cat("\n\nStep 2: Prior diagnostic check (checkPrior)\n")

prior_check <- checkPrior(
  tree = tree,
  Q = Q,
  lambda_sigma = 1.0,            # Same prior as simulation
  expected_background = 0.7,     # Expected 70% background
  nsim = 50                       # 50 prior-predictive simulations
)

print(prior_check)

cat("\nPrior-posterior overlap: ", sprintf("%.1f%%", prior_check$overlap_measure * 100), "\n")
if (prior_check$overlap_measure < 0.5) {
  cat("  ✓ Prior-posterior overlap is low; data are informative\n")
} else {
  cat("  ⚠ Prior-posterior overlap is high; consider adjusting lambda_sigma\n")
}


# ============================================================================
# STEP 3: Fit Bayesian model via MCMC
# ============================================================================

cat("\n\nStep 3: Bayesian MCMC fitting\n")

fit_bayes <- fitRateScape(
  tree = tree,
  data = data,
  Q = Q,
  estimate_Q = FALSE,            # Q is fixed (known from external data)
  lambda_sigma = 1.0,            # Prior on heterogeneity variance
  expected_background = 0.7,     # Prior: expect 70% background rate
  root_prior = "fitzjohn",       # Default FitzJohn root prior
  ngen = 5000,                   # MCMC generations
  burn_in = 1000,                # Burn-in period
  thin = 5,                      # Thin every 5 iterations
  tau_init = 0.3,                # Proposal SD for MH
  target_acceptance = 0.30,      # Target 30% MH acceptance
  seed = 42
)

cat("MCMC fitting complete.\n")
cat("Acceptance rate:", sprintf("%.1f%%", 100 * fit_bayes$acceptance_rate), "\n")


# ============================================================================
# STEP 4: Summarize results
# ============================================================================

cat("\n\nStep 4: Posterior summary\n")

summary_bayes <- summarizeRates(fit_bayes)
print(summary_bayes)

cat("\nBayes factor for heterogeneity vs. homogeneous rates:\n")
cat("  BF =", sprintf("%.2f", summary_bayes$bayes_factor), "\n")

if (summary_bayes$bayes_factor > 10) {
  cat("  Interpretation: STRONG evidence for rate heterogeneity\n")
} else if (summary_bayes$bayes_factor > 3) {
  cat("  Interpretation: WEAK evidence for rate heterogeneity\n")
} else {
  cat("  Interpretation: NO clear evidence for rate heterogeneity\n")
}


# ============================================================================
# STEP 5: Visualize rate-painted tree
# ============================================================================

cat("\n\nStep 5: Visualize rate-painted phylogeny\n")

png(filename = "ratescaped_tree.png", width = 800, height = 600)
plot(fit_bayes, main = "Rate-Painted Phylogeny (Bayesian Fit)")
dev.off()

cat("Rate-painted tree saved to 'ratescaped_tree.png'\n")


# ============================================================================
# STEP 6: Alternative: Fit ML discretized gamma model
# ============================================================================

cat("\n\nStep 6: Maximum likelihood (discretized gamma approach)\n")

fit_ml <- rateCategories(
  tree = tree,
  data = data,
  Q = Q,
  k_min = 2,
  k_max = 8,                     # Test 2-8 rate categories
  em_tol = 1e-6,                 # EM convergence tolerance
  em_maxiter = 100,              # Max EM iterations
  root_prior = "fitzjohn"
)

summary_ml <- summarizeRates(fit_ml)
print(summary_ml)

cat("\nOptimal number of rate categories (by BIC):", fit_ml$best_k, "\n")


# ============================================================================
# STEP 7: Compare models
# ============================================================================

cat("\n\nStep 7: Model comparison\n")

comparison <- compareModels(fit_bayes, fit_ml, names = c("Bayesian", "ML"))
print(comparison)


cat("\n\n=== Quick start example complete ===\n")
cat("For more details, see the RateScape documentation and vignettes.\n")
