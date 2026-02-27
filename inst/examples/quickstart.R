##########################################################
## RateScape Quick Start Example
##########################################################
## This script demonstrates both the Bayesian (spike-and-slab)
## and ML (discretized-gamma) approaches using a simple
## simulated dataset with known rate heterogeneity.

library(RateScape)
library(ape)

set.seed(2026)

## --- 1. Set up tree and model ---
tree <- rtree(50)

# 4-state equal-rates model
Q <- makeQ("mk", nstates = 4, params = list(rate = 0.3))

## --- 2. Simulate with known rate heterogeneity ---
# Create rate scalars: most branches at background,
# a few clade-specific accelerations
true_rates <- rep(1.0, Nedge(tree))

# Make edges 1-8 fast (3x background)
true_rates[1:8] <- 3.0

# Make edges 20-25 slow (0.2x background)
true_rates[20:25] <- 0.2

# Simulate data
data <- simRateScape(tree, Q, rates = true_rates)

## --- 3. Fit Bayesian model (spike-and-slab) ---
fit_bayes <- fitRateScape(
  tree, data, Q,
  ngen = 5000,       # Use more (e.g., 50000) for real analyses
  sample_freq = 5,
  burnin = 0.25,
  verbose = TRUE
)

print(fit_bayes)

## --- 4. Visualize ---
# Painted phylogeny
par(mfrow = c(1, 2))
plotRateTree(fit_bayes, main = "Posterior Mean Rates")
plotRateTree(fit_bayes, color_by = "shift_prob", main = "P(Rate Shift)")

## --- 5. Summarize branch-specific rates ---
rates_df <- summarizeRates(fit_bayes)
head(rates_df, 10)  # Top 10 most shifted branches

## --- 6. Model comparison ---
comp <- compareModels(fit_bayes)
print(comp)

## --- 7. ML approach (discretized gamma) ---
fit_ml <- rateCategories(tree, data, Q, ncat = 4)
cat("\nML Results:\n")
cat("  Alpha:", fit_ml$alpha, "\n")
cat("  Rate categories:", round(fit_ml$rate_categories, 3), "\n")
cat("  Log-likelihood:", fit_ml$loglik, "\n")
cat("  AIC:", fit_ml$AIC, "\n")

# Plot ML result
par(mfrow = c(1, 1))
plotRateTree(fit_ml, main = "ML Rate Estimates")

## --- 8. Chromosome number example ---
# Build a chromosome number Q matrix
Q_chr <- makeQ("chromosome", params = list(
  gain = 0.2,        # ascending dysploidy
  loss = 0.3,        # descending dysploidy
  polyploidy = 0.01, # genome doubling
  max_chrom = 20
))

# Simulate chromosome data with rate heterogeneity
chr_rates <- simRates(Nedge(tree), pi = 0.7, sigma2 = 0.5)
chr_data <- simRateScape(tree, Q_chr, rates = chr_rates, root_state = 8)

# Fit (short run for demo)
fit_chr <- fitRateScape(
  tree, chr_data, Q_chr,
  ngen = 5000,
  sample_freq = 5,
  verbose = TRUE
)

plotRateTree(fit_chr, main = "Chromosome Number Rate Variation")
